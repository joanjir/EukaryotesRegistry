#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genomes.py
Consulta y evaluación de ensamblajes genómicos mediante NCBI Datasets API:
- Recupera reportes de genomas asociados a un taxID
- Evalúa métricas estructurales y de calidad
- Selecciona el mejor ensamblaje basado en un puntaje ponderado
"""

from http_client import _get, _get_binary
from config import log_kv
import io, zipfile, json


# ===================== PONDERACIONES =====================

REFSEQ_SCORE = {"REFERENCE GENOME": 3, "REPRESENTATIVE GENOME": 2}
LEVEL_SCORE = {"COMPLETE GENOME": 4, "CHROMOSOME": 3, "SCAFFOLD": 2, "CONTIG": 1}

MIN_COVERAGE_REFSEQ = 40.0  # Cobertura mínima para RefSeq
MIN_COVERAGE_GENBANK = 30.0 # Cobertura mínima para GenBank

# ===================== POLÍTICA DE FILTRADO DURO =====================

# Aceptamos solo estas categorías RefSeq (maximizan probabilidad de proteoma):
REQUIRED_REFSEQ_CATEGORIES = {"REFERENCE GENOME", "REPRESENTATIVE GENOME"}

# Niveles de ensamblaje permitidos; "SCAFFOLD" solo si supera N50 fuertes.
ALLOWED_LEVELS = {"COMPLETE GENOME", "CHROMOSOME"}  # "SCAFFOLD" condicional

# Umbrales mínimos (ajusta por grupo si lo necesitas):
MIN_SCAFFOLD_N50_KB = 1000.0   # ≥ 1 Mb
MIN_CONTIG_N50_KB   = 100.0    # ≥ 100 kb

# Si True, además de filtros duros, verificaremos en el catálogo del ZIP que exista .faa
REQUIRE_PROTEOME = True


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
    Filtro duro para maximizar la probabilidad de que exista PROT_FASTA (protein.faa):
    - RefSeq category ∈ {Reference, Representative}
    - Genome level ∈ {Complete, Chromosome} (o Scaffold con N50 fuertes)
    - Cobertura mínima (reutiliza umbrales existentes)
    """
    refcat = (metrics.get("RefSeq category") or "").upper()
    level  = (metrics.get("Genome level") or "").upper()
    cov    = float(metrics.get("Genome coverage") or 0.0)
    sN50   = float(metrics.get("Scaffold N50 (kb)") or 0.0)
    cN50   = float(metrics.get("Contig N50 (kb)") or 0.0)

    # 1) Categoría RefSeq (clave para anotación)
    if refcat not in REQUIRED_REFSEQ_CATEGORIES:
        return False

    # 2) Nivel de ensamblaje (alto), o "Scaffold" con N50 fuertes
    if level not in ALLOWED_LEVELS:
        if level == "SCAFFOLD":
            if not (sN50 >= MIN_SCAFFOLD_N50_KB or cN50 >= MIN_CONTIG_N50_KB):
                return False
        else:
            return False

    # 3) Cobertura mínima (reutilizamos los umbrales que ya definiste)
    if refcat == "REFERENCE GENOME" and cov < MIN_COVERAGE_REFSEQ:
        return False
    if refcat == "REPRESENTATIVE GENOME" and cov < MIN_COVERAGE_GENBANK:
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

def has_protein_faa_in_catalog(accession: str) -> bool:
    """
    Descarga el ZIP con archivos .faa y verifica:
    1. Que exista al menos un archivo protein.faa
    2. Que la mayoría de las proteínas tengan descripciones válidas (no "unnamed", "hypothetical", etc.)
    """
    UNNAMED_KEYS = (
        "unnamed protein product",
        "hypothetical protein",
        "uncharacterized protein",
    )

    path = f"/genome/accession/{accession}/download"
    params = {
        "include_annotation_type": "PROT_FASTA",
        "hydrated": "FULLY_HYDRATED"
    }
    try:
        blob = _get_binary(path, params=params, accept="application/zip")
        with zipfile.ZipFile(io.BytesIO(blob), "r") as z:
            catalog_name = next((n for n in z.namelist() if n.endswith("dataset_catalog.json")), None)
            if not catalog_name:
                return False

            # --- verificar que haya al menos un archivo protein.faa ---
            faa_files = [n for n in z.namelist() if n.endswith(".faa")]
            if not faa_files:
                return False

            # --- abrir un .faa y revisar encabezados ---
            unnamed_count, total = 0, 0
            for faa_name in faa_files:
                with z.open(faa_name, "r") as fh:
                    for line in io.TextIOWrapper(fh, encoding="utf-8"):
                        if line.startswith(">"):
                            total += 1
                            header = line.lower()
                            if any(k in header for k in UNNAMED_KEYS):
                                unnamed_count += 1

                # Si ya hay suficientes para estimar proporción, no necesitamos leer todo
                if total >= 1000:
                    break

            if total == 0:
                return False

            ratio_unnamed = unnamed_count / total
            # Filtramos si más del 80% son "unnamed"/"hypothetical"/"uncharacterized"
            if ratio_unnamed > 0.8:
                log_kv("WARN", "Proteoma descartado por pobre anotación",
                       Accession=accession, UnnamedRatio=f"{ratio_unnamed:.2f}")
                return False

        return True

    except Exception as e:
        log_kv("ERROR", "Fallo verificando proteoma", Accession=accession, Error=str(e))
        return False

def pick_best_assembly(reports: list[dict]) -> dict | None:
    """
    Selecciona el mejor ensamblaje de una lista basándose en su puntaje total,
    aplicando primero filtros duros y (opcionalmente) la verificación de proteoma.
    """
    if not reports:
        return None

    evaluated = []
    for rep in reports:
        metrics = extract_metrics(rep)

        # 1) Filtros duros (categoría, nivel, N50, cobertura)
        if not passes_quality_filter(metrics):
            log_kv("WARN", "Descartado por filtros duros",
                   Accession=metrics.get("Accession"),
                   RefSeq=metrics.get("RefSeq category"),
                   Level=metrics.get("Genome level"),
                   Coverage=metrics.get("Genome coverage"),
                   ScN50=metrics.get("Scaffold N50 (kb)"),
                   CtN50=metrics.get("Contig N50 (kb)"))
            continue

        # 2) Verificación opcional de proteoma real en el catálogo
        if REQUIRE_PROTEOME:
            acc = metrics.get("Accession")
            if not has_protein_faa_in_catalog(acc):
                log_kv("WARN", "Descartado: sin protein.faa en catálogo", Accession=acc)
                continue

        # 3) Puntuar y acumular
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
