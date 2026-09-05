#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validated loader for installable production figure parameters.

The file is config/plot_params.json (JSONC: // line-comments + trailing commas
allowed); relocate it with env var SCMA_PLOT_PARAMS.

Contract:
  * A key that is ABSENT falls back to the built-in default (backward compatible).
  * A key that is PRESENT but malformed (wrong type / out of range / bad choice)
    fails LOUDLY with a ValueError naming the offending "group.key", never silently.

  API: num / boolean / text / num_list / str_list (validated accessors);
       load() (parsed dict) ; path() (resolved file path).
"""

import os
import re
import json

SELF = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PATH = os.path.join(SELF, "config", "plot_params.json")
_CACHE = {"path": None, "data": None}


def path():
    """Resolved params file path."""
    return os.environ.get("SCMA_PLOT_PARAMS") or DEFAULT_PATH


def _strip_jsonc(text):
    """Remove // line-comments (not inside strings) and trailing commas -> strict JSON.
    Deterministic and byte-for-byte equivalent to rflow/plot_params.R's stripper."""
    out = []
    for line in text.splitlines():
        in_str = esc = False
        buf = []
        i, n = 0, len(line)
        while i < n:
            ch = line[i]
            if in_str:
                buf.append(ch)
                if esc:
                    esc = False
                elif ch == "\\":
                    esc = True
                elif ch == '"':
                    in_str = False
            else:
                if ch == '"':
                    in_str = True
                    buf.append(ch)
                elif ch == "/" and i + 1 < n and line[i + 1] == "/":
                    break  # rest of line is a comment
                else:
                    buf.append(ch)
            i += 1
        out.append("".join(buf))
    txt = "\n".join(out)
    return re.sub(r",(\s*[}\]])", r"\1", txt)  # drop trailing commas


def load():
    """Parse the params file once (cached). Returns {} if the file is absent (all
    defaults). Raises ValueError if the file exists but is not valid JSONC."""
    p = path()
    if _CACHE["path"] == p and _CACHE["data"] is not None:
        return _CACHE["data"]
    data = {}
    if os.path.exists(p):
        with open(p, "r", encoding="utf-8") as fh:
            raw = fh.read()
        try:
            data = json.loads(_strip_jsonc(raw))
        except json.JSONDecodeError as e:
            raise ValueError(f"plot_params: {p} is not valid JSON/JSONC: {e}")
        if not isinstance(data, dict):
            raise ValueError(f"plot_params: {p} top level must be an object")
    _CACHE["path"], _CACHE["data"] = p, data
    return data


_MISSING = object()


def _raw(group, key):
    return load().get(group, {}).get(key, _MISSING)


def num(group, key, default, lo=None, hi=None):
    v = _raw(group, key)
    if v is _MISSING:
        return default
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise ValueError(f"plot_params: '{group}.{key}' must be a number, got {v!r}")
    if (lo is not None and v < lo) or (hi is not None and v > hi):
        raise ValueError(f"plot_params: '{group}.{key}'={v} out of range [{lo}, {hi}]")
    return v


def boolean(group, key, default):
    v = _raw(group, key)
    if v is _MISSING:
        return default
    if not isinstance(v, bool):
        raise ValueError(f"plot_params: '{group}.{key}' must be true/false, got {v!r}")
    return v


def text(group, key, default, choices=None):
    v = _raw(group, key)
    if v is _MISSING:
        return default
    if not isinstance(v, str):
        raise ValueError(f"plot_params: '{group}.{key}' must be a string, got {v!r}")
    if choices is not None and v not in choices:
        raise ValueError(f"plot_params: '{group}.{key}'={v!r} not in {sorted(choices)}")
    return v


def num_list(group, key, default, n=None, lo=None, hi=None):
    v = _raw(group, key)
    if v is _MISSING:
        return default
    if not isinstance(v, (list, tuple)) or any(
        isinstance(x, bool) or not isinstance(x, (int, float)) for x in v
    ):
        raise ValueError(
            f"plot_params: '{group}.{key}' must be a list of numbers, got {v!r}"
        )
    if n is not None and len(v) != n:
        raise ValueError(
            f"plot_params: '{group}.{key}' must have {n} numbers, got {len(v)}"
        )
    for x in v:
        if (lo is not None and x < lo) or (hi is not None and x > hi):
            raise ValueError(
                f"plot_params: '{group}.{key}' element {x} out of range [{lo}, {hi}]"
            )
    return list(v)


def str_list(group, key, default):
    v = _raw(group, key)
    if v is _MISSING:
        return default
    if not isinstance(v, (list, tuple)) or any(not isinstance(x, str) for x in v):
        raise ValueError(
            f"plot_params: '{group}.{key}' must be a list of strings, got {v!r}"
        )
    return list(v)


if __name__ == "__main__":
    print(f"[plot_params] {path()}")
    print(json.dumps(load(), indent=2, ensure_ascii=False))
