#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import zipfile
import pandas as pd
import requests
from tqdm import tqdm

# ---- Configuración ----
BASE_DIR   = "proteins_download"
EXCEL_FILE = "best_genomes_by_class.xlsx"
LOG_FILE   = "protein_download_api.log"
TIMEOUT    = 90
CHUNK      = 1024 * 1024
DELETE_ZIP_AFTER_EXTRACT = True  # ponlo en False si quieres conservar los .zip

INCLUDE_PARAMS = (
    "include_annotation_type=PROT_FASTA&"
    "hydrated=FULLY_HYDRATED"
)

def build_url(accession: str) -> str:
    return (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{accession}/download?{INCLUDE_PARAMS}")

def safe_species_name(name: str) -> str:
    return (name.strip()
            .replace(" ", "_").replace("/", "_").replace("\\", "_")
            .replace("(", "").replace(")", ""))

def download_zip(accession: str, folder: str) -> tuple[str, str]:
    """Descarga el ZIP fully-hydrated para el accession indicado."""
    url = build_url(accession)
    zip_path = os.path.join(folder, f"{accession}.zip")

    if os.path.exists(zip_path) and os.path.getsize(zip_path) > 50_000:
        return accession, "already_downloaded"

    try:
        with requests.get(url, stream=True, timeout=TIMEOUT,
                          headers={"Accept": "application/zip", "User-Agent": "datasets-script/1.0"}) as r:
            if r.status_code != 200:
                return accession, f"HTTP {r.status_code}"
            with open(zip_path, "wb") as f:
                for chunk in r.iter_content(CHUNK):
                    if chunk:
                        f.write(chunk)
        return accession, "downloaded"
    except Exception as e:
        return accession, f"error {e}"

def extract_only_proteins(zip_path: str, species: str) -> tuple[str, str]:
    """
    Busca en cualquier subcarpeta del ZIP archivos *.faa (p. ej. ncbi_dataset/data/<ACC>/protein.faa),
    concatena todos y escribe: <Especie>.proteins.fa
    """
    out_dir = os.path.dirname(zip_path)
    sp = safe_species_name(species)
    buf = io.BytesIO()
    found = 0

    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            for name in z.namelist():
                if name.endswith("/"):
                    continue
                if name.lower().endswith(".faa"):
                    with z.open(name, "r") as fh:
                        buf.write(fh.read())
                    found += 1

        if found == 0:
            # (opcional) borrar el zip si no sirve
            if DELETE_ZIP_AFTER_EXTRACT and os.path.exists(zip_path):
                try: os.remove(zip_path)
                except Exception: pass
            return (f"{sp}.fa", "no_faa_found")

        out_path = os.path.join(out_dir, f"{sp}.fa")
        with open(out_path, "wb") as out:
            out.write(buf.getvalue())

        if DELETE_ZIP_AFTER_EXTRACT:
            try:
                os.remove(zip_path)
            except Exception:
                pass

        return (out_path, f"ok_concat_{found}")

    except Exception as e:
        return (zip_path, f"error {e}")


def process_sheet(sheet_name: str, df: pd.DataFrame) -> list[tuple]:
    """Secuencial: descarga todos los ZIP y luego extrae SOLO protein.faa."""
    dst = os.path.join(BASE_DIR, sheet_name)
    os.makedirs(dst, exist_ok=True)
    logs = []

    print(f"\n🧫 Hoja: {sheet_name} → {len(df)} accesiones")

    # 1) Descarga
    for _, row in tqdm(df.iterrows(), total=len(df), desc=f"{sheet_name} download"):
        acc = str(row["Accession"]).strip()
        sp  = str(row["Species"]).strip()
        acc, st = download_zip(acc, dst)
        logs.append((sp, acc, f"download:{st}"))

    # 2) Extracción SOLO protein.faa
    for _, row in tqdm(df.iterrows(), total=len(df), desc=f"{sheet_name} extract"):
        acc = str(row["Accession"]).strip()
        sp  = str(row["Species"]).strip()
        zip_path = os.path.join(dst, f"{acc}.zip")
        if not os.path.exists(zip_path) and not DELETE_ZIP_AFTER_EXTRACT:
            logs.append((sp, acc, "extract:missing_zip"))
            continue
        out, st = extract_only_proteins(zip_path, sp)
        logs.append((sp, acc, f"extract:{st}:{os.path.basename(out)}"))

    return logs

def main():
    if not os.path.exists(EXCEL_FILE):
        print(f" No se encontró el archivo {EXCEL_FILE}")
        return

    xls = pd.ExcelFile(EXCEL_FILE)
    all_logs = []

    for sheet_name in xls.sheet_names:
        df = pd.read_excel(EXCEL_FILE, sheet_name=sheet_name)
        df = df.dropna(subset=["Accession", "Species"])
        if df.empty:
            continue
        all_logs.extend(process_sheet(sheet_name, df))

    with open(LOG_FILE, "w", encoding="utf-8") as fh:
        for sp, acc, st in all_logs:
            fh.write(f"{sp}\t{acc}\t{st}\n")

    print("\n Proceso finalizado.")
    print(f" Log: {LOG_FILE}")

if __name__ == "__main__":
    main()
