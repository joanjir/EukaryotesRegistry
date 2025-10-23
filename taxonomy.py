#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
taxonomy.py
Módulo para la resolución taxonómica dentro del sistema EukaryotesRegistry:
- Obtiene tax_id por nombre científico
- Recupera IDs relacionados (clases, especies)
- Resuelve nombres científicos asociados a IDs
"""

from itertools import islice
from http_client import _get, _post
from config import log_kv

# ===================== FUNCIONES TAXONÓMICAS =====================

def taxid_by_name(name: str) -> str:
    """
    Obtiene el tax_id de un organismo a partir de su nombre científico.
    Utiliza /taxonomy/name_report.
    """
    log_kv("INFO", "Resolviendo tax_id", Name=name)
    r = _post("/taxonomy/name_report", json={"taxons": [name]})
    reports = r.json().get("reports") or []
    if not reports or "taxonomy" not in reports[0]:
        log_kv("WARN", "TaxID no encontrado", Name=name)
        raise ValueError(f"No se encontró tax_id para '{name}'")
    tax_id = str(reports[0]["taxonomy"]["tax_id"])
    log_kv("INFO", "TaxID encontrado", Name=name, TaxID=tax_id)
    return tax_id


def related_ids(root_tax_id: str, rank_upper: str, page_size: int = 1000) -> list[str]:
    """
    Obtiene los taxIDs descendientes con un rango dado (CLASS, SPECIES, etc.).
    Endpoint: /taxonomy/taxon/{id}/related_ids?ranks=RANK
    """
    log_kv("INFO", "Buscando descendientes", TaxID=root_tax_id, Rank=rank_upper)
    out, token = [], None
    while True:
        params = {"ranks": rank_upper, "page_size": page_size}
        if token:
            params["page_token"] = token
        r = _get(f"/taxonomy/taxon/{root_tax_id}/related_ids", params=params)
        js = r.json()
        out.extend([str(t) for t in js.get("tax_ids", [])])
        token = js.get("next_page_token")
        if not token:
            break
    log_kv("INFO", "Descendientes encontrados", Count=len(out), Rank=rank_upper)
    return out


def names_for_ids(ids: list[str]) -> list[dict]:
    """
    Resuelve nombres científicos para una lista de taxIDs.
    Endpoint: /taxonomy/taxon/{ids}/name_report
    """
    log_kv("INFO", "Resolviendo nombres de taxIDs", Count=len(ids))
    ids = [str(i) for i in ids]
    out = []
    CHUNK = 200  # límite para evitar requests demasiado grandes
    it = iter(ids)
    while True:
        chunk = list(islice(it, CHUNK))
        if not chunk:
            break
        joined = ",".join(chunk)
        r = _get(f"/taxonomy/taxon/{joined}/name_report")
        reports = r.json().get("reports") or []
        for rep in reports:
            tax = rep.get("taxonomy", {})
            out.append({
                "tax_id": str(tax.get("tax_id")),
                "rank": (tax.get("rank") or "").upper(),
                "name": (tax.get("current_scientific_name") or {}).get("name")
            })
    log_kv("INFO", "Nombres resueltos", Resueltos=len(out))
    return out
