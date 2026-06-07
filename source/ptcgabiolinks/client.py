"""Small GDC API client with retries.

Requests honors http_proxy, https_proxy, and all_proxy from the environment,
which is useful for users behind restricted networks.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

GDC_API_ROOT = "https://api.gdc.cancer.gov"
DEFAULT_TIMEOUT = (10, 600)


@lru_cache(maxsize=1)
def get_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        status=5,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def api_url(path: str) -> str:
    if path.startswith("http://") or path.startswith("https://"):
        return path
    return f"{GDC_API_ROOT}/{path.lstrip('/')}"


def get_json(path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    response = get_session().get(api_url(path), params=params, timeout=DEFAULT_TIMEOUT)
    response.raise_for_status()
    return response.json()


def post_data(ids: list[str]) -> requests.Response:
    response = get_session().post(
        api_url("/data/"),
        json={"ids": ids},
        stream=True,
        timeout=DEFAULT_TIMEOUT,
    )
    response.raise_for_status()
    return response


def get_gdc_info() -> dict[str, Any]:
    return get_json("/status")


def is_serve_ok() -> bool:
    info = get_gdc_info()
    if info.get("status") != "OK":
        raise RuntimeError("GDC server down, try to use this package later")
    return True
