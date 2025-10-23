# -*- coding: utf-8 -*-
import os, sys, json, logging
from logging.handlers import RotatingFileHandler
from datetime import datetime, timezone

BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
USER_AGENT = "EukaryotesRegistry/1.0 (contact: your_email@domain)"
API_KEY = os.getenv("NCBI_DATASETS_API_KEY")  # opcional

# Concurrencia / red (conservador para evitar 429)
MAX_WORKERS = int(os.getenv("MAX_WORKERS", "4"))
MAX_PARALLEL_CALLS = int(os.getenv("MAX_PARALLEL_CALLS", "4"))
NET_TIMEOUT = int(os.getenv("NET_TIMEOUT", "45"))

# Logging
LOG_DIR = os.path.abspath(".")
LOG_FILE = os.path.join(LOG_DIR, "pipeline.log")
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()  # DEBUG/INFO/WARNING/ERROR

def setup_logging():
    os.makedirs(LOG_DIR, exist_ok=True)
    logger = logging.getLogger("pipeline")
    logger.setLevel(getattr(logging, LOG_LEVEL, logging.INFO))
    logger.handlers[:] = []

    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] [%(threadName)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(getattr(logging, LOG_LEVEL, logging.INFO))
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = RotatingFileHandler(LOG_FILE, maxBytes=5_000_000, backupCount=3, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger

logger = setup_logging()

def log_kv(level, msg, **kv):
    line = f"{msg} | " + " ".join(f"{k}={v}" for k,v in kv.items())
    getattr(logger, level.lower() if hasattr(logger, level.lower()) else "info")(line)
