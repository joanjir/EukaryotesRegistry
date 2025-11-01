#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, shutil
from pathlib import Path

BASE_DIR   = Path("proteins_download")
DEST_ROOT  = Path("proteins_unnamed")    # espejo de BASE_DIR para las "sin nombre"
STATE_FILE = Path(".unnamed_state.json") # estado por carpeta relativa a BASE_DIR
SUMMARY_LOG= Path("unnamed_summary.log")

FA_EXTS = {".fa", ".fasta", ".faa"}

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
    Separa un FASTA en dos: nombradas vs sin nombre. Devuelve contadores.
    """
    stats = {"total": 0, "named": 0, "unnamed": 0}
    with named_out.open("w", encoding="utf-8") as w_named, \
         unnamed_out.open("w", encoding="utf-8") as w_unnamed, \
         src.open("r", encoding="utf-8", errors="replace") as fh:
        write_to = None
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

def process_one_dir(dir_abs: Path, state: dict):
    """
    Procesa una carpeta con posibles FASTA.
    Refleja la ruta relativa a BASE_DIR en DEST_ROOT.
    """
    rel = dir_abs.relative_to(BASE_DIR)  # p.ej. animalia/Chordata/Mammalia
    key = str(rel)
    if state.get(key) == "DONE":
        append_summary(f"[SKIP] {rel} (already processed)")
        return

    fasta_files = [p for p in dir_abs.iterdir() if p.is_file() and p.suffix.lower() in FA_EXTS]
    if not fasta_files:
        state[key] = "DONE"
        save_state(state)
        return

    dest_dir = DEST_ROOT / rel
    dest_dir.mkdir(parents=True, exist_ok=True)

    total_files = moved_whole = rewrote_named = created_unnamed = 0
    total_named_seqs = total_unnamed_seqs = 0

    for fa in fasta_files:
        total_files += 1
        base = fa.stem
        named_tmp   = dir_abs / f"{base}.named.tmp"
        unnamed_tmp = dir_abs / f"{base}.unnamed.tmp"

        stats = split_named_vs_unnamed_file(fa, named_tmp, unnamed_tmp)
        total_named_seqs   += stats["named"]
        total_unnamed_seqs += stats["unnamed"]

        if stats["unnamed"] == 0:
            named_tmp.unlink(missing_ok=True)
            unnamed_tmp.unlink(missing_ok=True)
            append_summary(f"[{rel}] {fa.name}: 0 unnamed; keep original")
            continue

        if stats["named"] == 0:
            shutil.move(str(fa), str(dest_dir / fa.name))
            moved_whole += 1
            named_tmp.unlink(missing_ok=True)
            shutil.move(str(unnamed_tmp), str(dest_dir / f"{base}.unnamed.fa"))
            append_summary(f"[{rel}] {fa.name}: ALL unnamed → moved whole file")
            continue

        shutil.move(str(named_tmp), str(fa))
        rewrote_named += 1
        shutil.move(str(unnamed_tmp), str(dest_dir / f"{base}.unnamed.fa"))
        created_unnamed += 1

        append_summary(
            f"[{rel}] {fa.name}: named={stats['named']}, unnamed={stats['unnamed']} → kept named; moved unnamed"
        )

    state[key] = "DONE"
    save_state(state)

    append_summary(
        f"[SUMMARY {rel}] files={total_files}, moved_whole={moved_whole}, rewrote_named={rewrote_named}, "
        f"created_unnamed={created_unnamed}, named_seqs={total_named_seqs}, unnamed_seqs={total_unnamed_seqs}"
    )

def run_filter_recursive(base_dir: Path = BASE_DIR, force: bool = False):
    """
    Recorre recursivamente BASE_DIR y procesa toda carpeta que contenga FASTA.
    Usa estado por carpeta relativa (p.ej., animalia/Chordata/ClaseX).
    """
    state = load_state()
    candidate_dirs = set()
    for ext in FA_EXTS:
        for f in base_dir.rglob(f"*{ext}"):
            candidate_dirs.add(f.parent)

    if not candidate_dirs:
        print("No se encontraron FASTA bajo", base_dir)
        return

    for d in sorted(candidate_dirs):
        rel = d.relative_to(base_dir)
        if (not force) and state.get(str(rel)) == "DONE":
            append_summary(f"[SKIP] {rel} (already processed)")
            continue
        process_one_dir(d, state)

    print("Filtro de 'unnamed' completado. Revisa", SUMMARY_LOG)

# -------------------- uso directo --------------------
if __name__ == "__main__":
    run_filter_recursive(BASE_DIR, force=False)  # force=True para reprocesar
