#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py
Punto de entrada principal del sistema EukaryotesRegistry.
Ejecuta el pipeline completo por filo, seleccionando las mejores especies por clase.
"""

from pipeline import best_species_per_class, export_results
from config import log_header, log_line, log_kv
import time

# ===================== CONFIGURACIÓN GLOBAL =====================

PHYLA = [
    
    "Arthropoda",
    
]

SPECIES_PER_CLASS = 10000  # revisar todas las especies posibles
TOP_PER_CLASS = 3          # conservar las 3 mejores por clase

# ===================== EJECUCIÓN PRINCIPAL =====================

def run_pipeline(phylum_list: list[str], species_per_class: int | None = None, top_per_class: int | None = None):
    """
    Ejecuta el flujo completo para una lista de filos y exporta los resultados.
    """
    species_per_class = species_per_class or SPECIES_PER_CLASS
    top_per_class = top_per_class or TOP_PER_CLASS

    start_time = time.time()
    log_header("INICIO DEL PIPELINE EukaryotesRegistry")

    all_rows = []
    for phylum in phylum_list:
        log_kv("INFO", "Analizando filo", Phylum=phylum)
        try:
            rows = best_species_per_class(phylum, species_per_class=species_per_class, top_per_class=top_per_class)
            all_rows.extend(rows)
        except Exception as e:
            log_kv("ERROR", "Error procesando filo", Phylum=phylum, Error=str(e))

    if all_rows:
        export_results(all_rows)
    else:
        log_line("No se generaron resultados. Verifica pipeline.log.")

    total_time = round(time.time() - start_time, 2)
    log_kv("INFO", "PIPELINE FINALIZADO", TotalRows=len(all_rows), Tiempo=f"{total_time}s")

# ===================== EJECUCIÓN DIRECTA =====================

if __name__ == "__main__":
    try:
        run_pipeline(PHYLA)
    except KeyboardInterrupt:
        log_line("Ejecución interrumpida manualmente.")
    except Exception as e:
        log_kv("ERROR", "Fallo crítico del sistema", Error=str(e))
