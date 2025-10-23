# -*- coding: utf-8 -*-

import time
import csv
import requests
import pandas as pd
from itertools import islice

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
HEADERS = {
    "User-Agent": "EukaryotesRegistry/1.0 (contact: your_email@domain)",
    "Accept": "application/json"
}

# -------------------------- HTTP helpers --------------------------
def _get(path, params=None, timeout=60, tries=4):
    url = f"{BASE}{path}"
    for i in range(tries):
        r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
        if r.status_code not in (429,) and r.status_code < 500:
            r.raise_for_status()
            return r
        time.sleep(min(2 ** i, 8))
    r.raise_for_status()
    return r


def _post(path, json=None, timeout=60, tries=4):
    url = f"{BASE}{path}"
    for i in range(tries):
        r = requests.post(url, json=json, headers=HEADERS, timeout=timeout)
        if r.status_code not in (429,) and r.status_code < 500:
            r.raise_for_status()
            return r
        time.sleep(min(2 ** i, 8))
    r.raise_for_status()
    return r


# -------------------------- Taxonomy steps --------------------------
def taxid_by_name(name: str) -> str:
    """Obtiene tax_id del nombre científico."""
    rep = _post("/taxonomy/name_report", json={"taxons": [name]})
    reports = rep.json().get("reports") or []
    if not reports or "taxonomy" not in reports[0]:
        raise ValueError(f"No encontrado: {name}")
    return str(reports[0]["taxonomy"]["tax_id"])


def related_ids(root_tax_id: str, rank_upper: str, page_size: int = 1000) -> list[str]:
    """Obtiene taxIDs de un rango (clase, especie, etc.)."""
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
    """Obtiene nombres y rangos para una lista de taxIDs."""
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
        reports = r.json().get("reports") or []
        for rep in reports:
            tax = rep.get("taxonomy", {})
            out.append({
                "tax_id": str(tax.get("tax_id")),
                "rank": (tax.get("rank") or "").upper(),
                "name": (tax.get("current_scientific_name") or {}).get("name")
            })
    return out


# -------------------------- Genome helpers --------------------------
REFSEQ_SCORE = {"REFERENCE GENOME": 3, "REPRESENTATIVE GENOME": 2}
LEVEL_SCORE = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}


def genome_dataset_report_for_taxid(taxid: str) -> list[dict]:
    """Obtiene dataset_report por taxID (más seguro que por nombre)."""
    r = _get(f"/genome/taxon/{taxid}/dataset_report")
    return r.json().get("reports") or []


def pick_best_assembly(reports: list[dict]) -> dict | None:
    """Selecciona el mejor ensamblaje según nivel, categoría y N50."""
    if not reports:
        return None

    def score(rep):
        info = rep.get("assembly_info", {}) or {}
        stats = rep.get("assembly_stats", {}) or {}
        refcat = (info.get("refseq_category") or "").upper()
        level = info.get("assembly_level") or ""
        scaffold_n50 = float(stats.get("scaffold_n50") or 0.0) / 1000.0
        contig_n50 = float(stats.get("contig_n50") or 0.0) / 1000.0
        try:
            n_scaff = int(stats.get("number_of_scaffolds") or 0)
            scaff_pen = -n_scaff / 1e5
        except Exception:
            scaff_pen = 0.0
        recent = 0.0
        date = info.get("release_date") or info.get("submission_date")
        if date:
            try:
                recent = (int(str(date)[:4]) - 2000) / 100.0
            except Exception:
                pass
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
    """Extrae las métricas relevantes."""
    info = rep.get("assembly_info", {}) or {}
    stats = rep.get("assembly_stats", {}) or {}
    org = rep.get("organism", {}) or {}
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


# -------------------------- Pipeline principal --------------------------
def best_two_species_per_class(phylum_name: str, species_per_class: int = 2):
    phylum_id = taxid_by_name(phylum_name)
    print(f"Phylum '{phylum_name}' -> taxid {phylum_id}")

    class_ids = related_ids(phylum_id, rank_upper="CLASS")
    class_info = {x["tax_id"]: x["name"]
                  for x in names_for_ids(class_ids)
                  if x["rank"] == "CLASS" and x["name"]}
    print(f"  Clases encontradas: {len(class_info)}")

    rows = []
    for class_id, class_name in sorted(class_info.items(), key=lambda kv: kv[1].lower()):
        sp_ids = related_ids(class_id, rank_upper="SPECIES")
        if species_per_class and species_per_class > 0:
            sp_ids = sp_ids[:species_per_class]

        sp_meta = {x["tax_id"]: x["name"]
                   for x in names_for_ids(sp_ids)
                   if x["rank"] == "SPECIES" and x["name"]}

        # Filtrar "sp." o "uncultured"
        sp_meta = {tid: name for tid, name in sp_meta.items()
                   if all(bad not in name.lower() for bad in [" sp.", "uncultured", "environmental sample"])}

        print(f"  Clase '{class_name}': {len(sp_meta)} especies válidas")

        for sp_taxid, sp_name in sp_meta.items():
            try:
                reps = genome_dataset_report_for_taxid(sp_taxid)
                best = pick_best_assembly(reps)
                if not best:
                    print(f"    Especie '{sp_name}': sin ensamblajes válidos.")
                    continue
                m = extract_metrics(best)
                print(f"    ✓ {sp_name}: mejor ensamblaje {m['Accession']} ({m['Genome level']}, cobertura={m['Genome coverage']})")
                rows.append({
                    "Phylum": phylum_name,
                    "Class": class_name,
                    "Species": sp_name,
                    **{k: m[k] for k in ["Accession", "RefSeq category", "Genome level", "Genome coverage", "Contig N50 (kb)", "Scaffold N50 (kb)"]}
                })
            except Exception as e:
                print(f"    ⚠ Error en especie '{sp_name}': {e}")

    return rows


# -------------------------- Ejecución --------------------------
if __name__ == "__main__":
    PHYLA = ["Tracheophyta"]  # ejemplo
    N_SPECIES_PER_CLASS = 50

    header = [
        "Phylum", "Class", "Species",
        "Accession", "RefSeq category", "Genome level",
        "Genome coverage", "Contig N50 (kb)", "Scaffold N50 (kb)"
    ]

    all_rows = []
    for ph in PHYLA:
        try:
            rows = best_two_species_per_class(ph, species_per_class=N_SPECIES_PER_CLASS)
            all_rows.extend(rows)
            print(f"[OK] {ph}: {len(rows)} filas registradas")
        except Exception as e:
            print(f"[WARN] {ph}: error -> {e}")

    if all_rows:
        df = pd.DataFrame(all_rows, columns=header)
        df.sort_values(by=["Phylum", "Class", "Species"], inplace=True)
        df.to_excel("best_genomes_by_class.xlsx", index=False)
        print("[OK] Resultados guardados en best_genomes_by_class.xlsx")
    else:
        print("[WARN] No se encontraron genomas válidos para guardar.")
