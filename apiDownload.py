#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import json
import re
import ast
import zipfile
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

# ================= Configuración =================
BASE_DIR   = Path("proteins_download")
EXCEL_FILE = "best_genomes_by_class.xlsx"
PHYLOS_MD  = "phylos.md"                       # definiciones de reinos y sus listas
LOG_FILE   = "protein_download_api.log"
STATE_FILE = "protein_download_state.json"     # estado persistente
TIMEOUT    = 90
CHUNK      = 1024 * 1024
DELETE_ZIP_AFTER_EXTRACT = True                # borra el ZIP tras extraer

# Solo pedimos PROT_FASTA (ZIP mínimo)
INCLUDE_PARAMS = (
    "include_annotation_type=PROT_FASTA&"
    "hydrated=FULLY_HYDRATED"
)

SECTION_KEYS = {
    "ANIMALIA_PHYLA": "animalia",
    "PLANT_PHYLA": "plantae",
    "FUNGAL_PHYLA": "fungi",
    "CHROMISTA_PHYLA": "chromista",
    "PROTOZOA_PHYLA": "protozoa",
}

# ================ Utilidades generales ================

def build_url(accession: str) -> str:
    return (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{accession}/download?{INCLUDE_PARAMS}")

def safe_species_name(name: str) -> str:
    return (str(name).strip()
            .replace(" ", "_").replace("/", "_").replace("\\", "_")
            .replace("(", "").replace(")", ""))

def load_state(path: str = str(STATE_FILE)) -> dict:
    p = Path(path)
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}

def save_state(state: dict, path: str = str(STATE_FILE)):
    tmp = f"{path}.tmp"
    Path(tmp).write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")
    Path(path).write_text(Path(tmp).read_text(encoding="utf-8"), encoding="utf-8")
    Path(tmp).unlink(missing_ok=True)

def append_log(lines: list[str]):
    if not lines:
        return
    with open(LOG_FILE, "a", encoding="utf-8") as fh:
        for ln in lines:
            fh.write(ln.rstrip() + "\n")

# ================ Parser de phylos.md ================

def parse_phylos_md(md_path: str) -> dict[str, list[str]]:
    """
    Lee phylos.md y extrae listas tipo:
      ANIMALIA_PHYLA = [ "A", "B", ... ]
    Devuelve: { 'animalia': [...], 'plantae': [...], ... }
    """
    result = {v: [] for v in SECTION_KEYS.values()}
    p = Path(md_path)
    if not p.exists():
        return result

    text = p.read_text(encoding="utf-8")
    pattern = re.compile(
        r'^\s*([A-Z_]+)\s*=\s*\[(.*?)\]\s*$',
        flags=re.MULTILINE | re.DOTALL
    )

    for name, payload in pattern.findall(text):
        if name not in SECTION_KEYS:
            continue
        kingdom = SECTION_KEYS[name]
        literal = f'[{payload}]'
        try:
            items = ast.literal_eval(literal)
        except Exception:
            sanitized = literal.replace("“", '"').replace("”", '"').replace("’", "'")
            items = ast.literal_eval(sanitized)
        clean = [str(x).strip() for x in items if str(x).strip()]
        result[kingdom] = clean

    return result

# ================ Red / descarga / extracción ================

def download_zip(accession: str, folder: Path) -> tuple[str, str]:
    url = build_url(accession)
    zip_path = folder / f"{accession}.zip"

    if zip_path.exists() and zip_path.stat().st_size > 50_000:
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

def extract_only_proteins(zip_path: Path, species: str) -> tuple[str, str]:
    out_dir = zip_path.parent
    sp = safe_species_name(species)
    out_path = out_dir / f"{sp}.fa"

    if out_path.exists() and out_path.stat().st_size > 0:
        if DELETE_ZIP_AFTER_EXTRACT and zip_path.exists():
            try: zip_path.unlink()
            except Exception: pass
        return (str(out_path), "already_extracted")

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
            if DELETE_ZIP_AFTER_EXTRACT and zip_path.exists():
                try: zip_path.unlink()
                except Exception: pass
            return (str(out_path), "no_faa_found")

        with open(out_path, "wb") as out:
            out.write(buf.getvalue())

        if DELETE_ZIP_AFTER_EXTRACT and zip_path.exists():
            try: zip_path.unlink()
            except Exception: pass

        return (str(out_path), f"ok_concat_{found}")

    except Exception as e:
        return (str(zip_path), f"error {e}")

# ================ Proceso por clase (hoja) con estado ================

def process_class_sheet(kingdom: str, phylum: str, class_name: str, df: pd.DataFrame, state: dict) -> list[tuple]:
    """
    Filtra la hoja por 'Phylum' == phylum y procesa accesiones.
    Guarda en: proteins_download/<kingdom>/<phylum>/<class>/
    """
    logs = []

    if "Phylum" not in df.columns:
        logs.append((class_name, "-", f"skip:no_Phylum_column_in_sheet"))
        return logs

    sub = df[df["Phylum"].astype(str).str.strip() == phylum].dropna(subset=["Accession", "Species"])
    if sub.empty:
        return logs

    # *** IMPORTANTE: incluye la CLASE en el path ***
    class_dir = BASE_DIR / kingdom / phylum / class_name
    class_dir.mkdir(parents=True, exist_ok=True)

    # Estado por kingdom::phylum::class
    key = f"{kingdom}::{phylum}::{class_name}"
    class_state = state.setdefault(key, {})  # Accession -> status

    print(f"  • Clase: {class_name} ({len(sub)} accesiones)")

    # 1) Descarga
    for _, row in tqdm(sub.iterrows(), total=len(sub), desc=f"{class_name} download"):
        acc = str(row["Accession"]).strip()
        sp  = str(row["Species"]).strip()

        if class_state.get(acc) == "done":
            out_file = class_dir / f"{safe_species_name(sp)}.fa"
            if out_file.exists() and out_file.stat().st_size > 0:
                logs.append((sp, acc, "skip:already_done"))
                continue

        out_file = class_dir / f"{safe_species_name(sp)}.fa"
        if out_file.exists() and out_file.stat().st_size > 0:
            class_state[acc] = "done"
            logs.append((sp, acc, "skip:output_exists"))
            continue

        acc, st = download_zip(acc, class_dir)
        logs.append((sp, acc, f"download:{st}"))

    # 2) Extracción
    for _, row in tqdm(sub.iterrows(), total=len(sub), desc=f"{class_name} extract"):
        acc = str(row["Accession"]).strip()
        sp  = str(row["Species"]).strip()

        if class_state.get(acc) == "done":
            out_file = class_dir / f"{safe_species_name(sp)}.fa"
            if out_file.exists() and out_file.stat().st_size > 0:
                logs.append((sp, acc, "skip:already_done"))
                continue

        zip_path = class_dir / f"{acc}.zip"
        if not zip_path.exists() and not DELETE_ZIP_AFTER_EXTRACT:
            logs.append((sp, acc, "extract:missing_zip"))
            continue

        out, st = extract_only_proteins(zip_path, sp)
        logs.append((sp, acc, f"extract:{st}:{Path(out).name}"))

        if st.startswith("ok_") or st == "already_extracted" or st == "no_faa_found":
            class_state[acc] = "done"
        else:
            class_state[acc] = f"error:{st}"

        save_state(state)  # progresivo

    return logs

# ================ Orquestación principal ================

def main():
    if not Path(EXCEL_FILE).exists():
        print(f"❌ No se encontró el archivo {EXCEL_FILE}")
        return

    kingdom_map = parse_phylos_md(PHYLOS_MD)
    if not any(kingdom_map.values()):
        print("⚠️  No se pudo leer 'phylos.md' o no contiene listas válidas. No se crearán carpetas a ciegas.")

    state = load_state()
    xls = pd.ExcelFile(EXCEL_FILE)
    all_logs = []

    for kingdom, phyla in kingdom_map.items():
        if not phyla:
            continue
        print(f"\n🧭 Reino: {kingdom} ({len(phyla)} phyla)")

        for phylum in phyla:
            print(f"\n🧬 Phylum: {phylum}")
            processed_any = False

            for sheet_name in xls.sheet_names:
                df = pd.read_excel(EXCEL_FILE, sheet_name=sheet_name)
                if "Phylum" not in df.columns:
                    print(f"  - Hoja '{sheet_name}' sin columna 'Phylum' → se salta.")
                    continue

                logs = process_class_sheet(kingdom, phylum, sheet_name, df, state)
                if logs:
                    processed_any = True
                    all_logs.extend(logs)

            if not processed_any:
                print(f"  (sin filas para '{phylum}' en ninguna hoja; no se crean carpetas)")

    out_lines = [f"{sp}\t{acc}\t{st}" for (sp, acc, st) in all_logs]
    append_log(out_lines)

    print("\n✅ Proceso finalizado.")
    print(f"📘 Log: {LOG_FILE}")
    print(f"🗂️  Estado: {STATE_FILE}")

if __name__ == "__main__":
    main()
