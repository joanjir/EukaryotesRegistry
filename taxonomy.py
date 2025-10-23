# -*- coding: utf-8 -*-
from itertools import islice
from config import logger, log_kv
from http_client import http_get, http_post

def taxid_by_name(name: str) -> str:
    r = http_post("/taxonomy/name_report", json={"taxons": [name]})
    reps = r.json().get("reports") or []
    if not reps or "taxonomy" not in reps[0]:
        raise ValueError(f"No encontrado: {name}")
    tid = str(reps[0]["taxonomy"]["tax_id"])
    log_kv("info", "name_report", query=name, tax_id=tid, rank=reps[0]["taxonomy"].get("rank"))
    return tid

def related_ids(root_tax_id: str, rank_upper: str, page_size: int = 1000) -> list[str]:
    out, token = [], None
    pages = 0
    while True:
        params = {"ranks": rank_upper, "page_size": page_size}
        if token: params["page_token"] = token
        r = http_get(f"/taxonomy/taxon/{root_tax_id}/related_ids", params=params)
        js = r.json()
        got = [str(t) for t in js.get("tax_ids", [])]
        out.extend(got)
        token = js.get("next_page_token")
        pages += 1
        log_kv("debug", "related_ids_page", root=root_tax_id, rank=rank_upper, page=pages, got=len(got), total=len(out))
        if not token: break
    log_kv("info", "related_ids_done", root=root_tax_id, rank=rank_upper, total=len(out))
    return out

def names_for_ids(ids: list[str]) -> list[dict]:
    ids = [str(i) for i in ids]
    out, CHUNK = [], 50  # chunks pequeños para evitar 429
    it = iter(ids)
    idx = 0
    while True:
        chunk = list(islice(it, CHUNK))
        if not chunk: break
        idx += 1
        joined = ",".join(chunk)
        r = http_get(f"/taxonomy/taxon/{joined}/name_report")
        reps = r.json().get("reports") or []
        for rep in reps:
            tax = rep.get("taxonomy", {})
            out.append({
                "tax_id": str(tax.get("tax_id")),
                "rank": (tax.get("rank") or "").upper(),
                "name": (tax.get("current_scientific_name") or {}).get("name")
            })
        log_kv("debug", "name_report_chunk", chunk=idx, count=len(reps), acc=len(out))
    log_kv("info", "name_report_done", count=len(out))
    return out
