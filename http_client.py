# -*- coding: utf-8 -*-
import time, random, requests
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from threading import Semaphore

from config import BASE_URL, USER_AGENT, API_KEY, MAX_PARALLEL_CALLS, NET_TIMEOUT, logger, log_kv

sema = Semaphore(MAX_PARALLEL_CALLS)

def make_session():
    s = requests.Session()
    retries = Retry(
        total=5, connect=5, read=5, backoff_factor=0.8,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET","POST"]
    )
    s.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    if API_KEY:
        s.headers.update({"X-API-Key": API_KEY})
    adapter = HTTPAdapter(max_retries=retries, pool_connections=100, pool_maxsize=100)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s

SESSION = make_session()

def _call(method, path, *, params=None, json=None, timeout=NET_TIMEOUT):
    url = f"{BASE_URL}{path}"
    attempts = 6
    for i in range(attempts):
        with sema:
            r = SESSION.request(method, url, params=params, json=json, timeout=timeout)
        if r.status_code != 429:
            r.raise_for_status()
            return r
        ra = r.headers.get("Retry-After")
        wait = float(ra) if ra else (1.5 * (i + 1))
        jitter = random.random()
        log_kv("warning", "HTTP 429", url=url, wait_s=round(wait+jitter,2))
        time.sleep(wait + jitter)
    r.raise_for_status()

def http_get(path, *, params=None, timeout=NET_TIMEOUT):
    return _call("GET", path, params=params, timeout=timeout)

def http_post(path, *, json=None, timeout=NET_TIMEOUT):
    return _call("POST", path, json=json, timeout=timeout)
