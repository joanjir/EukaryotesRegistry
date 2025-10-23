#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Consulta NCBI por Filo -> Clases y, para cada clase, selecciona hasta 2 organismos
con mejores ensamblajes (priorizando RefSeq y nivel cromosómico/completo).
Devuelve métricas: Genome level, Genome coverage, Contig N50 (kb), Scaffold N50 (kb).

Uso:
  python plant_classes_best_genomes.py --phylum "Bryophyta" --per-class 2 --email "tu_correo@dominio"
  python plant_classes_best_genomes.py --phylum "Tracheophyta" --per-class 2 --csv salida.csv
"""

import argparse, time, sys
import requests
from urllib.parse import quote
from collections import defaultdict

# ---------- Utilidades ----------
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# Orden de preferencia
REFSEQ_RANK = {"reference genome": 3, "representative genome": 2}
ASSEMBLY_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

def sleep_backoff(i):
    time.sleep(min(2**i, 8))

def esearch_taxonomy(term, email=None, retmax=10000):
    """Devuelve lista de taxids desde NCBI Taxonomy (E-utilities)."""
    params = {
        "db": "taxonomy",
        "term": term,
        "retmode": "json",
        "retmax": retmax
    }
    if email:
        params["email"] = email
    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])

def esummary_taxonomy(taxids, email=None):
    """Devuelve {taxid: scientific_name, rank}."""
    out = {}
    if not taxids:
        return out
    chunks = [taxids[i:i+500] for i in range(0, len(taxids), 500)]
    for chunk in chunks:
        params = {
            "db": "taxonomy",
            "id": ",".join(chunk),
            "retmode": "json"
        }
        if email:
            params["email"] = email
        r = requests.get(f"{EUTILS_BASE}/esummary.fcgi", params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        for tid, rec in data.get("result", {}).items():
            if tid == "uids":
                continue
            sci = rec.get("scientificname")
            rank = rec.get("rank")
            out[tid] = {"name": sci, "rank": rank}
    return out

def get_classes_under_phylum(phylum, email=None):
    """
    Busca todas las clases dentro del filo dado (subárbol) en NCBI Taxonomy.
    Term de búsqueda: "<phylum>[Subtree] AND rank:class".
    """
    term = f'"{phylum}"[Subtree] AND rank:class'
    taxids = esearch_taxonomy(term, email=email, retmax=5000)
    meta = esummary_taxonomy(taxids, email=email)
    # Ordenar por nombre científico
    classes = sorted([(tid, meta[tid]["name"]) for tid in meta if meta[tid]["rank"] == "class"],
                     key=lambda x: x[1].lower())
    return classes  # lista de (taxid, class_name)

def datasets_list_genomes_for_taxon(taxid, page_token=None, only_refseq=True):
    """
    Lista ensamblados para un taxón usando NCBI Datasets v2.
    Filtra por RefSeq si only_refseq=True.
    """
    url = f"{DATASETS_BASE}/genome/taxon/{taxid}"
    params = {
        "returned_content": "COMPLETE",   # incluir assembly_info + assembly_stats
        "page_size": 100
    }
    if only_refseq:
        # Preferir RefSeq; el filtro 'refseq_only=true' limita a ensamblados RefSeq
        params["refseq_only"] = "true"
    if page_token:
        params["page_token"] = page_token

    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.json()

def _safe(n, factor=1.0):
    if n is None:
        return None
    try:
        return float(n) * factor
    except Exception:
        return None

def extract_metrics(asm):
    """Extrae métricas de un registro de datasets."""
    org = asm.get("organism", {}).get("organism_name")
    acc = asm.get("accession")
    info = asm.get("assembly", {}).get("assembly_info", {})
    stats = asm.get("assembly", {}).get("assembly_stats", {})
    ann  = asm.get("annotation_info", {}) or {}

    assembly_level = info.get("assembly_level")  # 'Complete Genome', 'Chromosome', etc.
    refseq_cat = info.get("refseq_category")  # 'reference genome', 'representative genome', None
    # Cobertura: distintos campos posibles según dataset
    coverage = info.get("sequencing_coverage") or info.get("coverage") or stats.get("coverage")
    # N50: en bases; convertir a kb
    contig_n50  = _safe(stats.get("contig_n50"), factor=1.0/1000.0)
    scaffold_n50 = _safe(stats.get("scaffold_n50"), factor=1.0/1000.0)
    num_scaff = stats.get("number_of_scaffolds") or stats.get("scaffold_count")
    num_contig = stats.get("number_of_contigs") or stats.get("contig_count")
    date = info.get("submission_date") or info.get("release_date")

    return {
        "organism": org,
        "accession": acc,
        "assembly_level": assembly_level,
        "refseq_category": refseq_cat,
        "coverage": coverage,
        "contig_n50_kb": contig_n50,
        "scaffold_n50_kb": scaffold_n50,
        "num_scaffolds": num_scaff,
        "num_contigs": num_contig,
        "date": date
    }

def quality_score(m):
    """Puntaje para ordenar: RefSeq > Assembly level > N50 > menor #scaffolds > fecha."""
    ref = REFSEQ_RANK.get((m.get("refseq_category") or "").lower(), 0)
    lvl = ASSEMBLY_RANK.get(m.get("assembly_level") or "", 0)
    n50 = m.get("scaffold_n50_kb") or 0.0
    # penalizar muchos scaffolds
    sc_penalty = 0.0
    try:
        sc = float(m.get("num_scaffolds") or 0)
        sc_penalty = -sc / 1e5  # penalización pequeña
    except Exception:
        pass
    # Prioriza registros más recientes levemente
    recent = 0.0
    if m.get("date"):
        try:
            y = int(str(m["date"])[:4])
            recent = (y - 2000) / 100.0
        except Exception:
            pass
    return (ref, lvl, n50, sc_penalty, recent)

def best_genomes_for_class(class_taxid, per_class=2, include_genbank_fallback=True):
    """
    Devuelve hasta 'per_class' genomas de mejor calidad para una clase (preferencia RefSeq).
    Si no hay RefSeq, puede caer a GenBank (include_genbank_fallback=True).
    """
    # 1) Intento RefSeq only
    recs = []
    page = None
    for i in range(4):
        data = datasets_list_genomes_for_taxon(class_taxid, page_token=page, only_refseq=True)
        for asm in data.get("assemblies", []):
            recs.append(extract_metrics(asm))
        page = data.get("next_page_token")
        if not page:
            break

    # 2) Si no hay suficientes y se permite fallback, traer GenBank también
    if include_genbank_fallback and len(recs) < per_class:
        page = None
        for i in range(3):
            data = datasets_list_genomes_for_taxon(class_taxid, page_token=page, only_refseq=False)
            for asm in data.get("assemblies", []):
                recs.append(extract_metrics(asm))
            page = data.get("next_page_token")
            if not page or len(recs) > 500:  # límite de seguridad
                break

    # Filtrar por niveles aceptables
    valid_levels = set(ASSEMBLY_RANK.keys())
    recs = [r for r in recs if (r.get("assembly_level") in valid_levels)]

    # Ordenar por calidad y tomar top-N
    recs.sort(key=quality_score, reverse=True)
    return recs[:per_class]

def main():
    ap = argparse.ArgumentParser(description="NCBI: Clases por Filo y mejores genomas por clase (RefSeq prior).")
    ap.add_argument("--phylum", required=True, help="Nombre del filo (e.g., 'Bryophyta', 'Tracheophyta').")
    ap.add_argument("--per-class", type=int, default=2, help="Máximo de organismos por clase.")
    ap.add_argument("--email", default=None, help="Correo para E-utilities (recomendado).")
    ap.add_argument("--csv", default=None, help="Ruta para exportar CSV (opcional).")
    args = ap.parse_args()

    # 1) Obtener clases
    classes = get_classes_under_phylum(args.phylum, email=args.email)
    if not classes:
        print(f"[WARN] No se encontraron clases bajo el filo '{args.phylum}'.")
        sys.exit(0)

    # 2) Para cada clase, seleccionar mejores genomas
    rows = []
    for tid, cname in classes:
        try:
            best = best_genomes_for_class(tid, per_class=args.per_class)
        except requests.HTTPError as e:
            print(f"[ERROR] Datasets API fallo para clase {cname} ({tid}): {e}", file=sys.stderr)
            continue
        if not best:
            rows.append({
                "Phylum": args.phylum, "Class": cname, "Organism": None,
                "Accession": None, "RefSeq category": None,
                "Genome level": None, "Genome coverage": None,
                "Contig N50 (kb)": None, "Scaffold N50 (kb)": None
            })
        for r in best:
            rows.append({
                "Phylum": args.phylum,
                "Class": cname,
                "Organism": r.get("organism"),
                "Accession": r.get("accession"),
                "RefSeq category": r.get("refseq_category"),
                "Genome level": r.get("assembly_level"),
                "Genome coverage": r.get("coverage"),
                "Contig N50 (kb)": r.get("contig_n50_kb"),
                "Scaffold N50 (kb)": r.get("scaffold_n50_kb"),
            })

    # 3) Salida en tabla
    if args.csv:
        import csv
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "Phylum", "Class", "Organism", "Accession",
                "RefSeq category", "Genome level", "Genome coverage",
                "Contig N50 (kb)", "Scaffold N50 (kb)"
            ])
            w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"[OK] Exportado CSV -> {args.csv}")
    else:
        # Imprimir tabla simple
        from textwrap import shorten
        hdr = ["Phylum","Class","Organism","Accession","RefSeq category","Genome level","Genome coverage","Contig N50 (kb)","Scaffold N50 (kb)"]
        print("\t".join(hdr))
        for row in rows:
            row["Organism"] = shorten(str(row["Organism"]), width=40, placeholder="…")
            print("\t".join([str(row.get(k,"")) for k in hdr]))

if __name__ == "__main__":
    main()
