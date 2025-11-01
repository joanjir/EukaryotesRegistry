#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
config.py
Configuración global del sistema EukaryotesRegistry:
- Variables de entorno (timeouts, concurrencia, API base)
- Inicialización del logger estructurado
- Soporte para logs individuales por filo, con fecha y hora exacta
"""

import logging
import os
from threading import Semaphore
from datetime import datetime

# ===================== CONFIGURACIÓN GLOBAL =====================

BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
USER_AGENT = "EukaryotesRegistry/1.0 (contact: your_email@domain)"
API_KEY = None  

# ===================== CONCURRENCIA =====================

MAX_WORKERS = 8
MAX_PARALLEL_CALLS = 8
NET_TIMEOUT = 45
sema = Semaphore(MAX_PARALLEL_CALLS)

# ===================== DIRECTORIO DE LOGS =====================

LOGS_DIR = "logsPhylas"
os.makedirs(LOGS_DIR, exist_ok=True)

# Logger principal (se redefinirá dinámicamente por filo)
logger = logging.getLogger("EukaryotesRegistry")
logger.setLevel(logging.INFO)
logger.propagate = False

# Remueve handlers previos para evitar duplicados
for h in list(logger.handlers):
    logger.removeHandler(h)

# ===================== FUNCIONES DE LOG =====================

def set_log_for_phylum(phylum_name: str, reino_name: str | None = None):
    """
    Configura un nuevo archivo de log exclusivo para un filo.
    Crea rutas jerárquicas del tipo:
        logsPhylas/<Reino>/<Phylum>_<fecha>.log
    Si no se especifica reino, usa 'Unclassified'.
    """
    global logger

    # === Definir reino y carpeta ===
    reino_name = reino_name or "Unclassified"
    reino_dir = os.path.join(LOGS_DIR, reino_name)
    os.makedirs(reino_dir, exist_ok=True)

    # === Crear nombre de archivo ===
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f")[:-3]
    safe_name = phylum_name.replace(" ", "_")
    log_filename = os.path.join(reino_dir, f"{safe_name}_{timestamp}.log")

    # === Limpiar handlers previos ===
    for h in list(logger.handlers):
        logger.removeHandler(h)

    # === Crear nuevo handler ===
    handler = logging.FileHandler(log_filename, mode="a", encoding="utf-8")
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # === Encabezado visual ===
    banner = (
        f"\n{'='*80}\n"
        f"PROCESANDO REINO {reino_name.upper()}  |  FILO {phylum_name.upper()}  -  "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"{'='*80}\n"
    )
    print(banner)
    logger.info(banner)

    return log_filename

def log_kv(level: str, msg: str, **kwargs):
    """
    Registra una entrada en el log con pares clave=valor.
    """
    context = " ".join([f"{k}={v}" for k, v in kwargs.items()])
    text = f"{msg} {context}".strip()

    if level.upper() == "INFO":
        logger.info(text)
    elif level.upper() in ("WARN", "WARNING"):
        logger.warning(text)
    elif level.upper() == "ERROR":
        logger.error(text)
    else:
        logger.debug(text)


def log_header(title: str):
    """Imprime y registra un encabezado visual."""
    banner = f"\n{'='*80}\n{title.upper()} - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n{'='*80}\n"
    print(banner)
    logger.info(banner)


def log_line(msg: str):
    """Imprime y registra una línea de información simple."""
    print(msg)
    logger.info(msg)
