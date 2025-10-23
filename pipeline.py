#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pipeline.py
Flujo principal del sistema EukaryotesRegistry:
- Para cada filo (phylum), obtiene clases y especies
- Descarga sus ensamblajes genómicos
- Evalúa la calidad y selecciona los mejores por clase
- Ejecuta en paralelo con control de carga y registro
"""

import csv
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from taxonomy import taxid_by_name, related_ids, names_for_ids
from genomes import genome_dataset_report_for_taxid, pick_best_assembly, extract_metrics, compute_score, passes_quality_filter
from config import MAX_WORKERS, log_kv, log_header, log_line


# ===================== FUNCIONES INTERNAS =====================

def best_species_rows_for_class(class_id: str, class_name: str, species_limit: int = 20) -> list[dict]:
    """
    Obtiene las mejores especies para una clase, verificando métricas de calidad.
    """
    rows = []
    try:
        species_ids = related_ids(class_id, rank_upper="SPECIES")
        if not species_ids:
            log_kv("WARN", "Clase sin especies detectadas", Class=class_name)
            return []

        species_ids = species_ids[:species_limit]
        species_meta = {
            x["tax_id"]: x["name"]
            for x in names_for_ids(species_ids)
            if x["rank"] == "SPECIES" and x["name"]
        }

        log_kv("INFO", "Analizando clase", Class=class_name, SpeciesCount=len(species_meta))

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {
                executor.submit(genome_dataset_report_for_taxid, sp_id): (sp_id, sp_name)
                for sp_id, sp_name in species_meta.items()
            }

            for future in as_completed(futures):
                sp_id, sp_name = futures[future]
                try:
                    reports = future.result()
                    best = pick_best_assembly(reports)
                    if not best:
                        continue

                    metrics = extract_metrics(best)

                    # Validar cobertura y métricas completas
                    if not passes_quality_filter(metrics):
                        log_kv("WARN", "Descartado por baja cobertura", Class=class_name, Species=sp_name)
                        continue

                    for key in ["Genome coverage", "Contig N50 (kb)", "Scaffold N50 (kb)"]:
                        if metrics.get(key) is None:
                            metrics[key] = 0.0

                    metrics["Score"] = compute_score(metrics)
                    rows.append({
                        "Class": class_name,
                        "Species": sp_name,
                        **{k: metrics.get(k) for k in [
                            "Accession", "RefSeq category", "Genome level",
                            "Genome coverage", "Contig N50 (kb)", "Scaffold N50 (kb)", "Score"
                        ]}
                    })
                except Exception as e:
                    log_kv("ERROR", "Error analizando especie", Class=class_name, Species=sp_name, Error=str(e))

    except Exception as e:
        log_kv("ERROR", "Error procesando clase", Class=class_name, Error=str(e))

    # Ordenar dentro de la clase por score descendente
    rows.sort(key=lambda r: r.get("Score", 0), reverse=True)
    return rows

def best_two_species_per_class(phylum_name: str, species_per_class: int = 20) -> list[dict]:
    """
    Para un filo dado:
    1. Obtiene sus clases.
    2. Evalúa sus especies.
    3. Conserva las 2 mejores especies por clase según score.
    """
    try:
        phylum_id = taxid_by_name(phylum_name)
    except Exception as e:
        log_kv("ERROR", "Filo no encontrado", Phylum=phylum_name, Error=str(e))
        return []

    log_header(f"Procesando filo {phylum_name}")
    class_ids = related_ids(phylum_id, rank_upper="CLASS")
    class_info = {
        x["tax_id"]: x["name"]
        for x in names_for_ids(class_ids)
        if x["rank"] == "CLASS" and x["name"]
    }

    all_rows = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(best_species_rows_for_class, cid, cname, species_per_class): (cid, cname)
            for cid, cname in class_info.items()
        }

        for future in as_completed(futures):
            cid, cname = futures[future]
            try:
                rows = future.result()
                if not rows:
                    continue

                # Ordenar por calidad si existe cobertura o N50
                def sort_key(r):
                    cov = r.get("Genome coverage") or 0
                    n50 = r.get("Scaffold N50 (kb)") or 0
                    return (cov, n50)

                top2 = sorted(rows, key=lambda r: r.get("Score", 0), reverse=True)[:2]
                for r in top2:
                    r["Phylum"] = phylum_name
                    r["Class"] = cname
                all_rows.extend(top2)

                log_kv("INFO", "Clase procesada", Phylum=phylum_name, Class=cname, Selected=len(top2))
            except Exception as e:
                log_kv("ERROR", "Error en clase", Class=cname, Error=str(e))

    log_kv("INFO", "Filo completado", Phylum=phylum_name, Rows=len(all_rows))
    return all_rows


def export_results(all_rows: list[dict], output_csv="best_genomes_by_class.csv", output_xlsx="best_genomes_by_class.xlsx"):
    """
    Exporta los resultados consolidados a CSV y Excel.
    """
    header = [
        "Phylum", "Class", "Species",
        "Accession", "RefSeq category", "Genome level",
        "Genome coverage", "Contig N50 (kb)", "Scaffold N50 (kb)"
    ]

    if not all_rows:
        log_line("⚠️ No se generaron resultados para exportar.")
        return

    # CSV
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(all_rows)
    log_kv("INFO", "Archivo CSV generado", File=output_csv, Rows=len(all_rows))

    # Excel
    df = pd.DataFrame(all_rows, columns=header)
    df.sort_values(by=["Phylum", "Class", "Species"], inplace=True, kind="stable")
    df.to_excel(output_xlsx, index=False)
    log_kv("INFO", "Archivo Excel generado", File=output_xlsx, Rows=len(df))
