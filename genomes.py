# -*- coding: utf-8 -*-
from typing import Tuple
from config import log_kv
from http_client import http_get

# Reglas y pesos
REFSEQ_SCORE = {"REFERENCE GENOME": 3, "REPRESENTATIVE GENOME": 2}
LEVEL_SCORE  = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
MIN_COV_REFSEQ = 40.0
MIN_COV_OTHER  = 30.0

# Boost para especies modelo (si aparecen)
MODEL_BOOST = {
    # planta modelo, gramíneas, solanáceas, musgos, hepáticas…
    "Arabidopsis thaliana": 2.0,
    "Oryza sativa": 2.0,
    "Nicotiana tabacum": 1.5,
    "Zea mays": 1.5,
    "Physcomitrium patens": 1.5,
    "Marchantia polymorpha": 1.5,
    "Sorghum bicolor": 1.2,
    "Glycine max": 1.2,
    "Triticum aestivum": 1.2,
}

def genome_dataset_report_for_taxid(tax_id: str) -> list[dict]:
    r = http_get(f"/genome/taxon/{tax_id}/dataset_report")
    reps = r.json().get("reports") or []
    # filtrar exactamente por tax_id
    out = []
    for rep in reps:
        org_tid = str((rep.get("organism") or {}).get("tax_id") or "")
        if org_tid == str(tax_id):
            out.append(rep)
    log_kv("info", "genome_reports", tax_id=tax_id, count=len(out))
    return out

def _f_or_none(x):
    try: return float(x)
    except: return None

def assess_assembly(rep: dict) -> Tuple[bool, str]:
    info  = rep.get("assembly_info", {}) or {}
    stats = rep.get("assembly_stats", {}) or {}

    refcat = (info.get("refseq_category") or "").upper()
    level  = info.get("assembly_level") or ""

    cov = _f_or_none(stats.get("genome_coverage"))
    c50 = _f_or_none(stats.get("contig_n50"))
    s50 = _f_or_none(stats.get("scaffold_n50"))
    nsc = _f_or_none(stats.get("number_of_scaffolds")) or 0.0

    # Cobertura mínima
    min_cov = MIN_COV_REFSEQ if refcat in {"REFERENCE GENOME","REPRESENTATIVE GENOME"} else MIN_COV_OTHER
    if cov is not None and cov < min_cov:
        return (False, f"low_coverage:{cov}x<{min_cov}x")

    # Si no hay coverage: acepta si el nivel/estructura son fuertes
    if cov is None:
        if level in {"Complete Genome","Chromosome"} and ((s50 and s50 >= 1_000_000) or (c50 and c50 >= 100_000)) and nsc <= 1000:
            return (True, "fallback_ok_no_coverage")
        return (False, "no_coverage_and_low_structure")

    # Con cobertura: pide mínimos razonables
    if level in {"Complete Genome","Chromosome"} or (s50 and s50 >= 500_000) or (c50 and c50 >= 50_000):
        return (True, "ok")

    return (False, "n50_low")

def score_assembly(rep: dict) -> float:
    info  = rep.get("assembly_info", {}) or {}
    stats = rep.get("assembly_stats", {}) or {}
    org   = rep.get("organism", {}) or {}

    refcat = (info.get("refseq_category") or "").upper()
    level  = info.get("assembly_level") or ""
    s_ref   = REFSEQ_SCORE.get(refcat, 0)
    s_level = LEVEL_SCORE.get(level, 0)
    s_scaf  = _f_or_none(stats.get("scaffold_n50")) or 0.0
    s_cont  = _f_or_none(stats.get("contig_n50")) or 0.0
    n_scaff = _f_or_none(stats.get("number_of_scaffolds")) or 0.0
    recent  = 0.0
    date = info.get("release_date") or info.get("submission_date")
    if date:
        try: recent = (int(str(date)[:4]) - 2000) / 100.0
        except: pass

    boost = MODEL_BOOST.get(org.get("organism_name",""), 0.0)

    return (10*s_ref) + (5*s_level) + (s_scaf/1e6) + (s_cont/1e6) - (n_scaff/1e6) + recent + boost

def pick_best_accepted(reports: list[dict]) -> tuple[dict|None, dict]:
    audited = {"accepted": [], "rejected": []}
    for rep in reports:
        ok, why = assess_assembly(rep)
        rec = {
            "accession": rep.get("current_accession") or rep.get("accession"),
            "refseq": (rep.get("assembly_info",{}) or {}).get("refseq_category"),
            "level":  (rep.get("assembly_info",{}) or {}).get("assembly_level"),
            "cov":    (rep.get("assembly_stats",{}) or {}).get("genome_coverage"),
            "contig_n50": (rep.get("assembly_stats",{}) or {}).get("contig_n50"),
            "scaffold_n50": (rep.get("assembly_stats",{}) or {}).get("scaffold_n50"),
            "why": why
        }
        (audited["accepted"] if ok else audited["rejected"]).append(rec)

    if not audited["accepted"]:
        return (None, audited)

    acc_reps = []
    acc_ids  = {a["accession"] for a in audited["accepted"]}
    for rp in reports:
        acc = rp.get("current_accession") or rp.get("accession")
        if acc in acc_ids:
            acc_reps.append(rp)

    acc_reps.sort(key=score_assembly, reverse=True)
    best = acc_reps[0]
    return (best, audited)

def extract_metrics(rep: dict) -> dict:
    info  = rep.get("assembly_info", {}) or {}
    stats = rep.get("assembly_stats", {}) or {}
    org   = rep.get("organism", {}) or {}
    def kb(x): 
        try: return float(x)/1000.0
        except: return None
    return {
        "Accession": rep.get("current_accession") or rep.get("accession"),
        "RefSeq category": info.get("refseq_category"),
        "Genome level": info.get("assembly_level"),
        "Genome coverage": stats.get("genome_coverage"),
        "Contig N50 (kb)": kb(stats.get("contig_n50")),
        "Scaffold N50 (kb)": kb(stats.get("scaffold_n50")),
        "Organism": org.get("organism_name"),
        "TaxID": org.get("tax_id")
    }
