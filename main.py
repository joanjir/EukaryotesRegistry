#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py
Pipeline principal del sistema EukaryotesRegistry.
Ejecuta el análisis por reinos y guarda una hoja Excel por filo
inmediatamente después de completarlo.
"""

from pipeline import best_species_per_class, export_results
from config import log_header, log_line, log_kv, set_log_for_phylum
from phylos import (
    ANIMALIA_PHYLA, PLANT_PHYLA, FUNGAL_PHYLA,
    CHROMISTA_PHYLA, PROTOZOA_PHYLA
)
import os, time

# ===================== CONFIGURACIÓN GLOBAL =====================

REINOS = {
    #"Animalia": ANIMALIA_PHYLA,
    "Plantae": PLANT_PHYLA,
    #"Fungi": FUNGAL_PHYLA,
    #"Chromista": CHROMISTA_PHYLA,
    #"Protozoa": PROTOZOA_PHYLA,
}

SPECIES_PER_CLASS = 10000
TOP_PER_CLASS = 3

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ===================== FUNCIÓN PRINCIPAL =====================

def run_pipeline():
    start_time = time.time()
    log_header("INICIO DEL PIPELINE EukaryotesRegistry")

    # --- Recorre cada reino ---
    for reino, phyla in REINOS.items():
        log_kv("INFO", "Procesando reino", Reino=reino)
        output_path = os.path.join(RESULTS_DIR, f"{reino}.xlsx")

        for phylum in phyla:
            set_log_for_phylum(phylum, reino)
            log_kv("INFO", "Analizando filo", Reino=reino, Phylum=phylum)

            try:
                rows = best_species_per_class(phylum, SPECIES_PER_CLASS, TOP_PER_CLASS)
                if not rows:
                    log_kv("WARN", "Sin resultados válidos para filo", Reino=reino, Phylum=phylum)
                    continue

                # ✅ Exportar inmediatamente la hoja del filo
                export_results(rows, output_excel=output_path)
                log_kv("INFO", "Filo exportado al Excel del reino",
                       Reino=reino, Phylum=phylum, Archivo=output_path, Filas=len(rows))

            except Exception as e:
                log_kv("ERROR", "Error procesando filo", Reino=reino, Phylum=phylum, Error=str(e))

        log_kv("INFO", "Reino completado", Reino=reino, Archivo=output_path)

    total_time = round(time.time() - start_time, 2)
    log_kv("INFO", "PIPELINE FINALIZADO", Tiempo=f"{total_time}s")
    log_line(f"Pipeline completado en {total_time}s. Archivos en '{RESULTS_DIR}/'.")

# ===================== EJECUCIÓN DIRECTA =====================

if __name__ == "__main__":
    try:
        run_pipeline()
    except KeyboardInterrupt:
        log_line("Ejecución interrumpida manualmente.")
    except Exception as e:
        log_kv("ERROR", "Fallo crítico del sistema", Error=str(e))
