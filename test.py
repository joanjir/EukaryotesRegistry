#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import time

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
HEADERS = {
    "User-Agent": "EukaryotesRegistry/1.0 (contact: your_email@domain)",
    "Accept": "application/json"
}

def query_name(name):
    try:
        r = requests.post(f"{BASE}/taxonomy/name_report", headers=HEADERS,
                          json={"taxons": [name]}, timeout=30)
        r.raise_for_status()
        js = r.json()
        reps = js.get("reports", [])
        if not reps:
            return None
        rep = reps[0]
        return rep.get("taxonomy", None)
    except Exception:
        return None


def check_name_in_ncbi(name: str):
    """
    Consulta el nombre en NCBI Taxonomy y devuelve estado, tax_id y rank.
    Incluye fallback por sinónimos.
    """
    tax = query_name(name)
    if tax:
        return {
            "name": name,
            "tax_id": tax.get("tax_id"),
            "rank": tax.get("rank"),
            "scientific_name": tax.get("current_scientific_name", {}).get("name"),
            "status": "found"
        }

    # Intentar búsqueda por sinónimos
    try:
        s = requests.get(f"{BASE}/taxonomy/search", params={"q": name}, headers=HEADERS, timeout=30)
        s.raise_for_status()
        js = s.json()
        hits = js.get("hits", [])
        if hits:
            h = hits[0].get("taxonomy", {})
            return {
                "name": name,
                "tax_id": h.get("tax_id"),
                "rank": h.get("rank"),
                "scientific_name": h.get("current_scientific_name", {}).get("name"),
                "status": "found (via search)"
            }
    except Exception as e:
        return {"name": name, "status": f"error: {e}"}

    return {"name": name, "status": "not_found"}


def check_list(names: list[str], delay=0.5):
    results = []
    for n in names:
        res = check_name_in_ncbi(n)
        print(f"{n:<25} -> {res['status']} ({res.get('tax_id', '-')})")
        results.append(res)
        time.sleep(delay)
    return results


if __name__ == "__main__":
    CHROMISTA_PHYLA = [
       
    "Orthonectida",  # verifique posible falta de ortografía
    
    ]

    print("🔍 Verificando nombres de filos en NCBI Taxonomy (modo extendido)...\n")
    results = check_list(CHROMISTA_PHYLA)

    found = [r for r in results if r["status"].startswith("found")]
    not_found = [r for r in results if not r["status"].startswith("found")]

    print("\n✅ Encontrados:", len(found))
    for f in found:
        print(f"  - {f['name']} (tax_id={f['tax_id']}, rank={f['rank']}, status={f['status']})")

    print("\n⚠️ No encontrados:", len(not_found))
    for nf in not_found:
        print(f"  - {nf['name']}: {nf['status']}")
