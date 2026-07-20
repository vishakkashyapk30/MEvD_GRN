"""I/O utilities: config loading (with single-level `inherit`), JSON/torch helpers."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict

import yaml


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge `override` into `base` (override wins). Returns a new dict."""
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | os.PathLike) -> Dict[str, Any]:
    """Load a YAML config. If it has an `inherit:` key, merge onto the parent
    config (resolved relative to this file's directory). One level of inheritance.
    """
    path = Path(path)
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    parent_name = cfg.pop("inherit", None)
    if parent_name:
        parent_path = (path.parent / parent_name).resolve()
        parent = load_config(parent_path)
        cfg = _deep_merge(parent, cfg)
    return cfg


def save_json(obj: Any, path: str | os.PathLike) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)


def load_json(path: str | os.PathLike) -> Any:
    with open(path) as f:
        return json.load(f)


def ensure_dir(path: str | os.PathLike) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p
