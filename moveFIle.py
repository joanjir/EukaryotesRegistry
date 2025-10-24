#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, json, shutil
from pathlib import Path

BASE_DIR = Path("proteins_download")
DEST_ROOT = Path("proteins_unnamed")           # a dónde mover los “sin nombre”
STATE_FILE = Path(".unnamed_state.json")       # registra carpetas ya procesadas
SUMMARY_LOG = Path("unnamed_summary.log")      # resumen de lo hecho

# Frases que marcan “sin nombre” en headers FASTA (minúsculas)
UNNAMED_KEYS = (
    "unnamed protein product",
    "hypothetical protein",
    "uncharacterized protein",
)

# -------------------- utilidades --------------------

def is_fasta_header(line: str) -> bool:
    return line.startswith(">")

def header_is_unnamed(header_line: str) -> bool:
    h = header_line.lower()
    return any(k in h for k in UNNAMED_KEYS)

def split_named_vs_unnamed_file(src: Path, named_out: Path, unnamed_out: Path):
    """
    Separa un FASTA en dos: nombradas vs sin nombre.
    Devuelve dict con contadores.
    """
    stats = {"total": 0, "named": 0, "unnamed": 0}
    w_named = named_out.open("w", encoding="utf-8")
    w_unnamed = unnamed_out.open("w", encoding="utf-8")

    write_to = None
    keep_current_in = None

    with src.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if is_fasta_header(line):
                stats["total"] += 1
                if header_is_unnamed(line):
                    stats["unnamed"] += 1
                    write_to = w_unnamed
                else:
                    stats["named"] += 1
                    write_to = w_named
                write_to.write(line)
            else:
                if write_to is not None:
                    write_to.write(line)

    w_named.close()
    w_unnamed.close()
    return stats

def load_state() -> dict:
    if STATE_FILE.exists():
        try:
            return json.loads(STATE_FILE.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}

def save_state(state: dict):
    STATE_FILE.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")

def append_summary(msg: str):
    with SUMMARY_LOG.open("a", encoding="utf-8") as log:
        log.write(msg.rstrip() + "\n")

# -------------------- procesamiento --------------------

def process_phylum_folder(folder: Path) -> dict:
    """
    Procesa una carpeta de phylum: separa secuencias y mueve las sin nombre.
    Regresa un dict resumen.
    """
    phylum = folder.name
    dest_phylum = DEST_ROOT / phylum
    dest_phylum.mkdir(parents=True, exist_ok=True)

    total_files = 0
    moved_whole = 0
    rewrote_named = 0
    created_unnamed = 0
    total_named_seqs = 0
    total_unnamed_seqs = 0

    fasta_exts = (".fa", ".fasta", ".faa")
    for fa in folder.glob("*"):
        if not fa.is_file() or fa.suffix.lower() not in fasta_exts:
            continue
        total_files += 1

        # archivos de salida
        base = fa.stem
        named_tmp = folder / f"{base}.named.tmp"
        unnamed_tmp = folder / f"{base}.unnamed.tmp"

        stats = split_named_vs_unnamed_file(fa, named_tmp, unnamed_tmp)
        total_named_seqs += stats["named"]
        total_unnamed_seqs += stats["unnamed"]

        if stats["unnamed"] == 0:
            # no hay “sin nombre”, dejar el original intacto
            named_tmp.unlink(missing_ok=True)
            unnamed_tmp.unlink(missing_ok=True)
            append_summary(f"[{phylum}] {fa.name}: 0 unnamed; keep original")
            continue

        if stats["named"] == 0:
            # todas sin nombre → mover archivo completo al destino y borrar el src
            dest = dest_phylum / fa.name
            shutil.move(str(fa), str(dest))
            moved_whole += 1
            # limpiar temporales
            named_tmp.unlink(missing_ok=True)
            # renombrar el tmp de unnamed como respaldo adicional (opcional)
            unnamed_tmp_dest = dest_phylum / f"{base}.unnamed.fa"
            shutil.move(str(unnamed_tmp), str(unnamed_tmp_dest))
            append_summary(f"[{phylum}] {fa.name}: ALL unnamed → moved whole file")
            continue

        # Mezcla: sobrescribimos el original con solo “named”
        shutil.move(str(named_tmp), str(fa))
        rewrote_named += 1

        # y movemos las sin nombre a la carpeta de destino
        unnamed_dest = dest_phylum / f"{base}.unnamed.fa"
        shutil.move(str(unnamed_tmp), str(unnamed_dest))
        created_unnamed += 1

        append_summary(
            f"[{phylum}] {fa.name}: named={stats['named']}, unnamed={stats['unnamed']} → kept named; moved unnamed"
        )

    return {
        "phylum": phylum,
        "total_files": total_files,
        "moved_whole": moved_whole,
        "rewrote_named": rewrote_named,
        "created_unnamed": created_unnamed,
        "total_named_seqs": total_named_seqs,
        "total_unnamed_seqs": total_unnamed_seqs,
    }

def run_filter(base_dir: Path = BASE_DIR, force: bool = False):
    """
    Revisa cada carpeta de phylum bajo base_dir.
    - Si ya fue procesada (según STATE_FILE) y no se usa force=True, la salta.
    - Escribe resumen en SUMMARY_LOG y estado en STATE_FILE.
    """
    state = load_state()
    overall = []

    for phylum_dir in sorted(p for p in base_dir.iterdir() if p.is_dir()):
        key = str(phylum_dir.resolve())
        if (not force) and state.get(key) == "DONE":
            append_summary(f"[SKIP] {phylum_dir.name} (already processed)")
            continue

        stats = process_phylum_folder(phylum_dir)
        overall.append(stats)
        state[key] = "DONE"
        save_state(state)

    # resumen general
    append_summary("=== RESUMEN GLOBAL ===")
    for s in overall:
        append_summary(
            f"{s['phylum']}: files={s['total_files']}, moved_whole={s['moved_whole']}, "
            f"rewrote_named={s['rewrote_named']}, created_unnamed={s['created_unnamed']}, "
            f"named_seqs={s['total_named_seqs']}, unnamed_seqs={s['total_unnamed_seqs']}"
        )
    print("Filtro de 'unnamed' completado. Revisa", SUMMARY_LOG)

# -------------------- uso directo --------------------
if __name__ == "__main__":
    run_filter(BASE_DIR, force=False)  # cambia a True si quieres re-procesar todo
