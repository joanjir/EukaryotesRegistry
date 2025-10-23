#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import csv
import requests
import pandas as pd
from itertools import islice
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from threading import Semaphore

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
USER_AGENT = "EukaryotesRegistry/1.0 (contact: your_email@domain)"
API_KEY = None  # si tienes una API key de NCBI, colócala aquí (string)

# ===== Concurrencia y sesión robusta =====
MAX_WORKERS = 8           # hilos totales
MAX_PARALLEL_CALLS = 8    # semáforo para limitar llamadas simultáneas
NET_TIMEOUT = 45          # segundos por llamada
sema = Semaphore(MAX_PARALLEL_CALLS)

def make_session():
    s = requests.Session()
    retries = Retry(
        total=5, connect=5, read=5, backoff_factor=0.8,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"]
    )
    s.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    if API_KEY:
        s.headers.update({"X-API-Key": API_KEY})
    adapter = HTTPAdapter(max_retries=retries, pool_connections=100, pool_maxsize=100)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s

SESSION = make_session()

def _call(method, path, *, params=None, json=None, timeout=NET_TIMEOUT):
    url = f"{BASE}{path}"
    with sema:
        r = SESSION.request(method, url, params=params, json=json, timeout=timeout)
    r.raise_for_status()
    return r

def _get(path, *, params=None, timeout=NET_TIMEOUT):
    return _call("GET", path, params=params, timeout=timeout)

def _post(path, *, json=None, timeout=NET_TIMEOUT):
    return _call("POST", path, json=json, timeout=timeout)

# ===== TAXONOMY =====
def taxid_by_name(name: str) -> str:
    r = _post("/taxonomy/name_report", json={"taxons": [name]})
    rep = (r.json().get("reports") or [])
    if not rep or "taxonomy" not in rep[0]:
        raise ValueError(f"No encontrado: {name}")
    return str(rep[0]["taxonomy"]["tax_id"])

def related_ids(root_tax_id: str, rank_upper: str, page_size: int = 1000) -> list[str]:
    out, token = [], None
    while True:
        params = {"ranks": rank_upper, "page_size": page_size}
        if token:
            params["page_token"] = token
        r = _get(f"/taxonomy/taxon/{root_tax_id}/related_ids", params=params)
        js = r.json()
        out.extend([str(t) for t in js.get("tax_ids", [])])
        token = js.get("next_page_token")
        if not token:
            break
    return out

def names_for_ids(ids: list[str]) -> list[dict]:
    ids = [str(i) for i in ids]
    out = []
    CHUNK = 200
    it = iter(ids)
    while True:
        chunk = list(islice(it, CHUNK))
        if not chunk:
            break
        joined = ",".join(chunk)
        r = _get(f"/taxonomy/taxon/{joined}/name_report")
        reps = r.json().get("reports") or []
        for rep in reps:
            tax = rep.get("taxonomy", {})
            out.append({
                "tax_id": str(tax.get("tax_id")),
                "rank": (tax.get("rank") or "").upper(),
                "name": (tax.get("current_scientific_name") or {}).get("name")
            })
    return out

# ===== GENOMES =====
REFSEQ_SCORE = {"REFERENCE GENOME": 3, "REPRESENTATIVE GENOME": 2}
LEVEL_SCORE  = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

def genome_dataset_report_for_taxid(tax_id: str) -> list[dict]:
    # 1) pedir por tax_id (evita ambigüedades)
    r = _get(f"/genome/taxon/{tax_id}/dataset_report")
    reps = r.json().get("reports") or []
    # 2) filtrar por pertenencia exacta (opcionalmente: o por set de descendientes)
    out = []
    for rep in reps:
        org_tid = str((rep.get("organism") or {}).get("tax_id") or "")
        if org_tid == str(tax_id):
            out.append(rep)
    return out


def pick_best_assembly(reports: list[dict]) -> dict | None:
    if not reports:
        return None
    def score(rep):
        info  = rep.get("assembly_info", {}) or {}
        stats = rep.get("assembly_stats", {}) or {}
        refcat = (info.get("refseq_category") or "").upper()
        level  = info.get("assembly_level") or ""
        scaffold_n50 = float(stats.get("scaffold_n50") or 0.0) / 1000.0
        contig_n50   = float(stats.get("contig_n50") or 0.0) / 1000.0
        try:
            n_scaff = int(stats.get("number_of_scaffolds") or 0)
            scaff_pen = -n_scaff / 1e5
        except Exception:
            scaff_pen = 0.0
        recent = 0.0
        date = info.get("release_date") or info.get("submission_date")
        if date:
            try: recent = (int(str(date)[:4]) - 2000) / 100.0
            except Exception: pass
        return (
            REFSEQ_SCORE.get(refcat, 0),
            LEVEL_SCORE.get(level, 0),
            scaffold_n50,
            contig_n50,
            scaff_pen,
            recent
        )
    return sorted(reports, key=score, reverse=True)[0]

def extract_metrics(rep: dict) -> dict:
    info  = rep.get("assembly_info", {}) or {}
    stats = rep.get("assembly_stats", {}) or {}
    org   = rep.get("organism", {}) or {}
    def f(x): return None if x is None else float(x) / 1000.0
    return {
        "Accession": rep.get("current_accession") or rep.get("accession"),
        "RefSeq category": info.get("refseq_category"),
        "Genome level": info.get("assembly_level"),
        "Genome coverage": stats.get("genome_coverage"),
        "Contig N50 (kb)": f(stats.get("contig_n50")),
        "Scaffold N50 (kb)": f(stats.get("scaffold_n50")),
        "Organism": org.get("organism_name"),
        "TaxID": org.get("tax_id")
    }

# ===== PIPELINE =====
def best_species_rows_for_class(class_id: str, class_name: str, species_limit: int) -> list[dict]:
    # 1) listar species taxIDs de la clase
    sp_ids = related_ids(class_id, rank_upper="SPECIES")
    if species_limit and species_limit > 0:
        sp_ids = sp_ids[:species_limit]

    # 2) resolver nombres (opcional, para mostrar Species)
    sp_meta = {x["tax_id"]: x["name"]
               for x in names_for_ids(sp_ids)
               if x["rank"] == "SPECIES" and x["name"]}

    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {}
        for sp_tid, sp_name in sp_meta.items():
            futs[ex.submit(genome_dataset_report_for_taxid, sp_tid)] = (sp_tid, sp_name)

        for fut in as_completed(futs):
            sp_tid, sp_name = futs[fut]
            try:
                reps = fut.result()
            except Exception as e:
                # timeout u otro error: registra fila vacía
                results.append({
                    "Class": class_name,
                    "Species": sp_name,
                    "Accession": None, "RefSeq category": None, "Genome level": None,
                    "Genome coverage": None, "Contig N50 (kb)": None, "Scaffold N50 (kb)": None
                })
                continue

            best = pick_best_assembly(reps)
            if not best:
                results.append({
                    "Class": class_name,
                    "Species": sp_name,
                    "Accession": None, "RefSeq category": None, "Genome level": None,
                    "Genome coverage": None, "Contig N50 (kb)": None, "Scaffold N50 (kb)": None
                })
            else:
                m = extract_metrics(best)
                results.append({
                    "Class": class_name,
                    "Species": sp_name,
                    "Accession": m["Accession"],
                    "RefSeq category": m["RefSeq category"],
                    "Genome level": m["Genome level"],
                    "Genome coverage": m["Genome coverage"],
                    "Contig N50 (kb)": m["Contig N50 (kb)"],
                    "Scaffold N50 (kb)": m["Scaffold N50 (kb)"]
                })
    return results

def best_two_species_per_class(phylum_name: str, species_per_class: int = 2) -> list[dict]:
    phylum_id = taxid_by_name(phylum_name)

    # IDs de clases bajo el filo/clado
    class_ids = related_ids(phylum_id, rank_upper="CLASS")
    class_info = {x["tax_id"]: x["name"]
                  for x in names_for_ids(class_ids)
                  if x["rank"] == "CLASS" and x["name"]}

    all_rows = []
    # paralelizar por clase también
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {}
        for cid, cname in class_info.items():
            futs[ex.submit(best_species_rows_for_class, cid, cname, species_per_class)] = (cid, cname)
        for fut in as_completed(futs):
            cid, cname = futs[fut]
            try:
                rows = fut.result()
            except Exception as e:
                rows = []
            for r in rows:
                r["Phylum"] = phylum_name
                r["Class"]  = r.get("Class") or cname
                all_rows.append(r)

    # ordenar un poco
    all_rows.sort(key=lambda r: (r.get("Phylum",""), r.get("Class",""), r.get("Species","")))
    return all_rows

# ===== MAIN =====
if __name__ == "__main__":
    PHYLA = [
       
        "Rhodophyta",
        "Tracheophyta",
    ]
    N_SPECIES_PER_CLASS = 20  # ajusta si quieres más

    header = [
        "Phylum","Class","Species",
        "Accession","RefSeq category","Genome level",
        "Genome coverage","Contig N50 (kb)","Scaffold N50 (kb)"
    ]

    all_rows = []
    for ph in PHYLA:
        try:
            rows = best_two_species_per_class(ph, species_per_class=N_SPECIES_PER_CLASS)
            all_rows.extend(rows)
            print(f"[OK] {ph}: {len(rows)} filas")
        except Exception as e:
            print(f"[WARN] {ph}: error -> {e}")

    # CSV
    with open("best_genomes_by_class.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(all_rows)
    print("[OK] CSV -> best_genomes_by_class.csv")

    # Excel
    df = pd.DataFrame(all_rows, columns=header)
    df.sort_values(by=["Phylum","Class","Species"], inplace=True, kind="stable")
    df.to_excel("best_genomes_by_class.xlsx", index=False)
    print("[OK] Excel -> best_genomes_by_class.xlsx")
