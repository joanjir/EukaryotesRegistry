#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py
Punto de entrada principal del sistema EukaryotesRegistry.
Orquesta la ejecución completa del pipeline taxonómico y genómico.
"""

from pipeline import best_two_species_per_class, export_results
from config import log_header, log_line, log_kv
import time

# ===================== CONFIGURACIÓN GLOBAL =====================
#Aca iria la lista de filos a analizar se pudiera cargar desde un archivo externo si se desea
PHYLA = [
    "Bryophyta",
    "Chlorophyta",
    "Charophyta",
    "Anthocerotophyta",
    "Glaucophyta",
    "Marchantiophyta",
    "Rhodophyta",
    "Tracheophyta"
    ]

SPECIES_PER_CLASS = 500  # número máximo de especies analizadas por clase
TOP_PER_CLASS = 3        # número de especies a conservar por clase


# ===================== EJECUCIÓN PRINCIPAL =====================

def run_pipeline(phylum_list: list[str], species_per_class: int = 20, top_per_class: int = 2):
    """
    Ejecuta el flujo completo para una lista de filos y exporta los resultados.
    """
    start_time = time.time()
    log_header("INICIO DEL PIPELINE EukaryotesRegistry")

    all_rows = []
    for phylum in phylum_list:
        log_kv("INFO", "Analizando filo", Phylum=phylum)
        try:
            rows = best_two_species_per_class(phylum, species_per_class=species_per_class)
            all_rows.extend(rows)
            log_kv("INFO", "Filo completado", Phylum=phylum, Resultados=len(rows))
        except Exception as e:
            log_kv("ERROR", "Error procesando filo", Phylum=phylum, Error=str(e))

    if all_rows:
        export_results(all_rows)
    else:
        log_line("No se generaron resultados. Verifica los registros de pipeline.log.")

    total_time = round(time.time() - start_time, 2)
    log_kv("INFO", "Pipeline finalizado", TotalRows=len(all_rows), Tiempo=f"{total_time}s")


# ===================== EJECUCIÓN DIRECTA =====================

if __name__ == "__main__":
    try:
        run_pipeline(PHYLA, species_per_class=SPECIES_PER_CLASS, top_per_class=TOP_PER_CLASS)
    except KeyboardInterrupt:
        log_line(" Ejecución interrumpida manualmente por el usuario.")
    except Exception as e:
        log_kv("ERROR", "Fallo crítico del sistema", Error=str(e))
