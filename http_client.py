#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
http_client.py
Cliente HTTP centralizado para el sistema EukaryotesRegistry:
- Manejo de sesiones persistentes
- Reintentos automáticos con backoff exponencial
- Control de concurrencia mediante semáforo
- Registro estructurado de solicitudes y errores
"""

import time
import requests
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from config import BASE_URL, USER_AGENT, API_KEY, NET_TIMEOUT, sema, log_kv

# ===================== SESIÓN HTTP GLOBAL =====================

def make_session() -> requests.Session:
    """Crea una sesión HTTP persistente con reintentos y cabeceras comunes."""
    session = requests.Session()

    retries = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=0.8,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"]
    )

    adapter = HTTPAdapter(max_retries=retries, pool_connections=100, pool_maxsize=100)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    session.headers.update({
        "User-Agent": USER_AGENT,
        "Accept": "application/json"
    })

    if API_KEY:
        session.headers.update({"X-API-Key": API_KEY})

    return session

SESSION = make_session()

# ===================== FUNCIONES DE PETICIÓN =====================

def _call(method: str, path: str, *, params=None, json=None, timeout=NET_TIMEOUT):
    """
    Realiza una llamada HTTP con control de concurrencia y registro.
    Se usa internamente por _get() y _post().
    """
    url = f"{BASE_URL}{path}"
    with sema:
        try:
            response = SESSION.request(method, url, params=params, json=json, timeout=timeout)
            response.raise_for_status()
            log_kv("INFO", f"{method} {path}", status=response.status_code)
            return response
        except requests.exceptions.RequestException as e:
            log_kv("ERROR", f"HTTP {method} error", path=path, error=str(e))
            raise

def _get(path: str, *, params=None, timeout=NET_TIMEOUT):
    """Realiza una petición GET segura al API."""
    return _call("GET", path, params=params, timeout=timeout)

def _post(path: str, *, json=None, timeout=NET_TIMEOUT):
    """Realiza una petición POST segura al API."""
    return _call("POST", path, json=json, timeout=timeout)

# ===================== PRUEBA RÁPIDA =====================

if __name__ == "__main__":
    # Ejemplo de prueba rápida del cliente
    try:
        print("[TEST] Realizando GET /taxonomy/taxon/2759/name_report")
        r = _get("/taxonomy/taxon/2759/name_report")
        print(f"Estado: {r.status_code} | Bytes recibidos: {len(r.content)}")
    except Exception as e:
        print("[ERROR]", e)
