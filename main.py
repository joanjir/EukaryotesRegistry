# -*- coding: utf-8 -*-
import csv, pandas as pd, time
from config import logger
from pipeline import best_two_species_per_class

def run_pipeline(phylum_list, species_per_class=20,
                 csv_path="best_genomes_by_class.csv",
                 xlsx_path="best_genomes_by_class.xlsx"):
    header = [
        "Phylum","Class","Species",
        "Accession","RefSeq category","Genome level",
        "Genome coverage","Contig N50 (kb)","Scaffold N50 (kb)"
    ]

    all_rows = []
    t0 = time.time()
    logger.info(f"[MAIN] Starting with {len(phylum_list)} phyla")
    for ph in phylum_list:
        rows = best_two_species_per_class(ph, species_per_class=species_per_class)
        all_rows.extend(rows)
        logger.info(f"[MAIN] {ph}: {len(rows)} filas")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(all_rows)
    logger.info(f"[MAIN] CSV -> {csv_path}")

    df = pd.DataFrame(all_rows, columns=header)
    df.sort_values(by=["Phylum","Class","Species"], inplace=True, kind="stable")
    df.to_excel(xlsx_path, index=False)
    logger.info(f"[MAIN] Excel -> {xlsx_path}")
    logger.info(f"[MAIN] Finished. Total rows={len(all_rows)} in {round(time.time()-t0,1)}s")

if __name__ == "__main__":
    PHYLA = [
        
        "Glaucophyta",
        "Rhodophyta",
    ]
    # Explora varias especies por clase y luego quedas con 2 mejores
    run_pipeline(PHYLA, species_per_class=5)
