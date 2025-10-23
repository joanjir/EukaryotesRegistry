#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genomes.py
Consulta y evaluación de ensamblajes genómicos mediante NCBI Datasets API:
- Recupera reportes de genomas asociados a un taxID
- Evalúa métricas estructurales y de calidad
- Selecciona el mejor ensamblaje basado en un puntaje ponderado
"""

from http_client import _get
from config import log_kv

# ===================== PONDERACIONES =====================

REFSEQ_SCORE = {"REFERENCE GENOME": 3, "REPRESENTATIVE GENOME": 2}
LEVEL_SCORE = {"COMPLETE GENOME": 4, "CHROMOSOME": 3, "SCAFFOLD": 2, "CONTIG": 1}

MIN_COVERAGE_REFSEQ = 40.0  # Cobertura mínima para RefSeq
MIN_COVERAGE_GENBANK = 30.0 # Cobertura mínima para GenBank

# ===================== FUNCIONES PRINCIPALES =====================

def genome_dataset_report_for_taxid(tax_id: str) -> list[dict]:
    """
    Recupera los ensamblajes genómicos asociados a un tax_id específico.
    Endpoint: /genome/taxon/{tax_id}/dataset_report
    """
    log_kv("INFO", "Consultando genomas", TaxID=tax_id)
    r = _get(f"/genome/taxon/{tax_id}/dataset_report")
    reports = r.json().get("reports") or []
    if not reports:
        log_kv("WARN", "Sin genomas disponibles", TaxID=tax_id)
    return reports


def extract_metrics(rep: dict) -> dict:
    """
    Extrae las métricas relevantes del reporte de ensamblaje.
    """
    info = rep.get("assembly_info", {}) or {}
    stats = rep.get("assembly_stats", {}) or {}
    org = rep.get("organism", {}) or {}

    def f(x):  # convierte a kilobases
        try:
            return None if x is None else float(x) / 1000.0
        except Exception:
            return None

    return {
        "Accession": rep.get("current_accession") or rep.get("accession"),
        "RefSeq category": info.get("refseq_category"),
        "Genome level": info.get("assembly_level"),
        "Genome coverage": float(stats.get("genome_coverage") or 0.0),
        "Contig N50 (kb)": f(stats.get("contig_n50")),
        "Scaffold N50 (kb)": f(stats.get("scaffold_n50")),
        "Number of scaffolds": int(stats.get("number_of_scaffolds") or 0),
        "Organism": org.get("organism_name"),
        "TaxID": org.get("tax_id")
    }


def passes_quality_filter(metrics: dict) -> bool:
    """
    Evalúa si un ensamblaje cumple los requisitos mínimos de cobertura.
    """
    refcat = (metrics.get("RefSeq category") or "").upper()
    cov = metrics.get("Genome coverage") or 0.0

    if refcat == "REFERENCE GENOME" and cov < MIN_COVERAGE_REFSEQ:
        return False
    if refcat != "REFERENCE GENOME" and cov < MIN_COVERAGE_GENBANK:
        return False
    return True


def compute_score(metrics: dict) -> float:
    """
    Calcula el puntaje compuesto (ponderado) de un ensamblaje.
    """
    refcat = (metrics.get("RefSeq category") or "").upper()
    level = (metrics.get("Genome level") or "").upper()
    scaffold_n50 = metrics.get("Scaffold N50 (kb)") or 0.0
    contig_n50 = metrics.get("Contig N50 (kb)") or 0.0
    n_scaff = metrics.get("Number of scaffolds") or 0
    coverage = metrics.get("Genome coverage") or 0.0

    score = (
        REFSEQ_SCORE.get(refcat, 0) * 1.0 +
        LEVEL_SCORE.get(level, 0) * 1.2 +
        (coverage / 50.0) +                 # mayor cobertura = mejor
        (scaffold_n50 + contig_n50) / 100.0 -
        (n_scaff / 1e5)                     # penalización por fragmentación
    )

    return round(score, 3)


def pick_best_assembly(reports: list[dict]) -> dict | None:
    """
    Selecciona el mejor ensamblaje de una lista basándose en su puntaje total.
    """
    if not reports:
        return None

    evaluated = []
    for rep in reports:
        metrics = extract_metrics(rep)
        if not passes_quality_filter(metrics):
            log_kv("WARN", "Descartado por baja cobertura",
                   Accession=metrics["Accession"],
                   Coverage=metrics["Genome coverage"])
            continue

        metrics["Score"] = compute_score(metrics)
        evaluated.append((metrics["Score"], rep, metrics))

    if not evaluated:
        return None

    evaluated.sort(key=lambda x: x[0], reverse=True)
    best_score, best_rep, best_metrics = evaluated[0]
    log_kv("INFO", "Mejor ensamblaje seleccionado",
           Accession=best_metrics["Accession"],
           Score=best_score,
           Coverage=best_metrics["Genome coverage"])
    return best_rep
