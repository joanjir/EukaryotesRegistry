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
        
def _get_binary(path: str, *, params=None, timeout=NET_TIMEOUT, accept: str = "application/octet-stream") -> bytes:
    """
    Realiza una petición GET y devuelve el cuerpo binario (bytes).
    - Usa SESSION y semáforo (sema) para limitar concurrencia.
    - Permite ajustar el header 'Accept' (útil para ZIPs).
    - Descarga en chunks para no saturar memoria innecesariamente.
    Lanza excepción si la respuesta no es 2xx.
    """
    url = f"{BASE_URL}{path}"
    with sema:
        try:
            # Clonar headers de la sesión y forzar 'Accept' para este request
            headers = dict(SESSION.headers)
            headers["Accept"] = accept

            resp = SESSION.get(
                url,
                params=params,
                timeout=timeout,
                headers=headers,
                stream=True,
            )
            status = resp.status_code
            log_kv("INFO", "HTTP GET(bin)", url=url, status=status)

            # Levanta excepción para códigos no exitosos
            resp.raise_for_status()

            # Descargar en chunks a memoria
            chunks = []
            for chunk in resp.iter_content(chunk_size=1024 * 1024):  # 1 MiB
                if chunk:
                    chunks.append(chunk)

            data = b"".join(chunks)
            log_kv("INFO", "HTTP GET(bin) ok", url=url, bytes=len(data))
            return data

        except requests.exceptions.RequestException as e:
            # Log detallado y re-lanzar para que el caller decida
            try:
                # Si hay respuesta, incluir código/status si existe
                sc = getattr(e.response, "status_code", None)
                log_kv("ERROR", "HTTP GET(bin) failed", url=url, status=sc, error=str(e))
            except Exception:
                log_kv("ERROR", "HTTP GET(bin) failed", url=url, error=str(e))
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
