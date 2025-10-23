#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
config.py
Configuración global del sistema EukaryotesRegistry:
- Variables de entorno (timeouts, concurrencia, API base)
- Inicialización del logger estructurado
- Constantes globales de control
"""

import logging
from threading import Semaphore
from datetime import datetime

# ===================== CONFIGURACIÓN GLOBAL =====================

BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
USER_AGENT = "EukaryotesRegistry/1.0 (contact: your_email@domain)"
API_KEY = None  # Si posees una clave de NCBI, colócala aquí

# ===================== CONCURRENCIA =====================

MAX_WORKERS = 8             # Número máximo de hilos paralelos bajar de ser necesario para evitar un exceso de carga
MAX_PARALLEL_CALLS = 8      # Límite de llamadas simultáneas al API
NET_TIMEOUT = 45            # Tiempo máximo por solicitud en segundos
sema = Semaphore(MAX_PARALLEL_CALLS)

# ===================== LOGGING ESTRUCTURADO =====================

LOG_FILE = "pipeline.log"

logging.basicConfig(
    filename=LOG_FILE,
    filemode="a",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def log_kv(level: str, msg: str, **kwargs):
    """
    Registra una entrada en el log con pares clave-valor.
    Ejemplo:
        log_kv("INFO", "Procesando clase", Class="Bryopsida", SpeciesCount=25)
    """
    context = " ".join([f"{k}={v}" for k, v in kwargs.items()])
    text = f"{msg} {context}".strip()
    if level.upper() == "INFO":
        logging.info(text)
    elif level.upper() == "WARN" or level.upper() == "WARNING":
        logging.warning(text)
    elif level.upper() == "ERROR":
        logging.error(text)
    else:
        logging.debug(text)

def log_header(title: str):
    """Imprime y registra un encabezado visual en el log."""
    banner = f"\n{'='*80}\n{title.upper()} - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n{'='*80}\n"
    print(banner)
    logging.info(banner)

def log_line(msg: str):
    """Imprime y registra una línea de información simple."""
    print(msg)
    logging.info(msg)
