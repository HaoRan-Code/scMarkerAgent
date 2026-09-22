#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from .config import (
    CACHE_DIR as CACHE,
    CONFIG_PATH,
    LLM_SETTINGS,
)

MODEL = os.environ.get("SCMA_LLM_MODEL", LLM_SETTINGS.model)
MODEL_OVERRIDE = os.environ.get("SCMA_LLM_MODEL", "").strip()
CONFIG_SHA256 = hashlib.sha256(Path(CONFIG_PATH).read_bytes()).hexdigest()
_RAW_LOG_LOCK = threading.Lock()
DEFAULT_OPENAI_BASE_URL = "https://api.openai.com/v1"


class LlmError(RuntimeError):
    pass


class LlmTransportError(LlmError):
    pass


class LlmResponseError(LlmError):
    pass


def _responses_endpoint(base_url: str) -> str:
    base = str(base_url or DEFAULT_OPENAI_BASE_URL).strip().rstrip("/")
    return base if base.endswith("/responses") else base + "/responses"


def resolve_model(model=None) -> str:
    return MODEL_OVERRIDE or str(model or LLM_SETTINGS.model)


def resolve_api_with_source() -> tuple[str, str, str]:
    if os.environ.get("SCMA_OFFLINE", "0").lower() in {"1", "true", "yes", "on"}:
        return "", "", "offline"
    key = os.environ.get("OPENAI_API_KEY", "").strip()
    base_url = os.environ.get("OPENAI_BASE_URL", DEFAULT_OPENAI_BASE_URL).strip()
    return (
        key,
        _responses_endpoint(base_url),
        "environment" if key else "none",
    )


def resolve_api() -> tuple[str, str]:
    key, url, _ = resolve_api_with_source()
    return key, url


def credential_source() -> str:
    return resolve_api_with_source()[2]


def provider_mode() -> str:
    configured = os.environ.get("OPENAI_BASE_URL", "").strip().rstrip("/")
    default = DEFAULT_OPENAI_BASE_URL.rstrip("/")
    return (
        "official_openai"
        if not configured or configured == default
        else "custom_openai_base"
    )


def _safe_error(value, api_key="", api_url="") -> str:
    text = str(value or "")
    for secret in (api_key, api_url):
        if secret:
            text = text.replace(secret, "<redacted>")
    text = re.sub(
        r"(?i)authorization\s*[:=]\s*bearer\s+\S+",
        "Authorization: Bearer <redacted>",
        text,
    )
    text = re.sub(r"https?://[^\s\"']+", "<redacted_url>", text)
    return text[:2000]


def _raw_log_path() -> Path:
    value = os.environ.get("SCMA_LLM_RAW_LOG", "")
    return Path(value) if value else Path(CACHE) / "llm_cold_calls.jsonl"


def _write_raw_record(record):
    path = _raw_log_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with _RAW_LOG_LOCK:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")


def parse_json(text):
    value = re.sub(
        r"```$",
        "",
        re.sub(r"^```(json)?", "", text.strip()).strip(),
    ).strip()
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        pass
    start = value.find("{")
    if start < 0:
        return None
    try:
        return json.JSONDecoder().raw_decode(value[start:])[0]
    except json.JSONDecodeError:
        return None


def _response_text(payload):
    if not isinstance(payload, dict):
        return None
    if isinstance(payload.get("output_text"), str):
        return payload["output_text"]
    for item in payload.get("output") or []:
        if not isinstance(item, dict):
            continue
        for content in item.get("content") or []:
            if not isinstance(content, dict):
                continue
            if content.get("type") in {"output_text", "text"}:
                text = content.get("text")
                if isinstance(text, str) and text.strip():
                    return text
    return None


def call_llm(
    prompt,
    api_url,
    api_key,
    timeout=None,
    retries=None,
    reasoning_effort=None,
    model=None,
    trace_id=None,
    turn_index=None,
):
    timeout = int(timeout or LLM_SETTINGS.timeout_seconds)
    retries = int(retries or LLM_SETTINGS.retries)
    resolved_model = resolve_model(model)
    per_model_key = os.environ.get(f"SCMA_API_KEY__{resolved_model}", "").strip()
    cred_source = credential_source()
    if per_model_key:
        api_key = per_model_key
        cred_source = "environment(per-model)"
    body = {
        "model": resolved_model,
        "input": prompt,
    }
    if reasoning_effort:
        body["reasoning"] = {"effort": reasoning_effort}
    headers = {
        "Authorization": "Bearer " + api_key,
        "Content-Type": "application/json",
    }
    last_error = "no attempt"
    prompt_sha256 = hashlib.sha256(prompt.encode("utf-8")).hexdigest()
    for attempt in range(retries):
        started = datetime.now(timezone.utc)
        response = None
        provider_payload = None
        content = None
        try:
            response = requests.post(
                api_url,
                headers=headers,
                json=body,
                timeout=timeout,
            )
            try:
                provider_payload = response.json()
            except requests.exceptions.JSONDecodeError:
                provider_payload = {"raw_text": response.text}
            if response.status_code == 200:
                content = _response_text(provider_payload)
                if content and content.strip():
                    _write_raw_record(
                        {
                            "schema_version": "llm-cold-call-v1",
                            "timestamp_utc": started.isoformat(),
                            "prompt_sha256": prompt_sha256,
                            "prompt": prompt,
                            "model": resolved_model,
                            "reasoning_effort": reasoning_effort,
                            "trace_id": trace_id,
                            "turn_index": turn_index,
                            "config_sha256": CONFIG_SHA256,
                            "credential_source": cred_source,
                            "provider_mode": provider_mode(),
                            "request": body,
                            "attempt": attempt + 1,
                            "http_status": response.status_code,
                            "provider_response": provider_payload,
                            "content": content,
                            "parsed_json": parse_json(content),
                            "usage": provider_payload.get("usage"),
                            "cost": {
                                "value": None,
                                "currency": None,
                                "reason": "pricing_not_configured",
                            },
                            "error": None,
                        }
                    )
                    return content, None
                last_error = "empty content"
            else:
                last_error = f"HTTP {response.status_code}"
        except requests.RequestException as exc:
            last_error = _safe_error(repr(exc), api_key, api_url)
        last_error = _safe_error(last_error, api_key, api_url)
        _write_raw_record(
            {
                "schema_version": "llm-cold-call-v1",
                "timestamp_utc": started.isoformat(),
                "prompt_sha256": prompt_sha256,
                "prompt": prompt,
                "model": resolved_model,
                "reasoning_effort": reasoning_effort,
                "trace_id": trace_id,
                "turn_index": turn_index,
                "config_sha256": CONFIG_SHA256,
                "credential_source": cred_source,
                "provider_mode": provider_mode(),
                "request": body,
                "attempt": attempt + 1,
                "http_status": response.status_code if response is not None else None,
                "provider_response": provider_payload,
                "content": content,
                "parsed_json": parse_json(content) if content else None,
                "usage": (
                    provider_payload.get("usage")
                    if isinstance(provider_payload, dict)
                    else None
                ),
                "cost": {
                    "value": None,
                    "currency": None,
                    "reason": "pricing_not_configured",
                },
                "error": last_error,
            }
        )
        if attempt + 1 < retries:
            time.sleep(min(2**attempt, 30))
    return None, _safe_error(last_error, api_key, api_url)


@dataclass
class _ResponseCacheState:
    path: str | None = None
    data: dict[str, str] | None = None
    dirty: bool = False
    hits: int = 0
    misses: int = 0
    credential_misses: int = 0
    lock: Any = field(default_factory=threading.Lock)


_RESPONSE_CACHE = _ResponseCacheState()


def response_cache_path():
    path = os.environ.get("SCMA_LLM_CACHE", "")
    return path if path else os.path.join(CACHE, "llm_response_cache.json")


def _load_vcache() -> dict[str, str]:
    path = response_cache_path()
    with _RESPONSE_CACHE.lock:
        if _RESPONSE_CACHE.path == path and _RESPONSE_CACHE.data is not None:
            return _RESPONSE_CACHE.data
        data: dict[str, str] = {}
        if os.path.exists(path):
            try:
                with open(path, encoding="utf-8") as handle:
                    loaded = json.load(handle)
            except (OSError, json.JSONDecodeError) as error:
                raise ValueError(f"invalid LLM response cache: {path}") from error
            if not isinstance(loaded, dict):
                raise ValueError(f"LLM response cache is not an object: {path}")
            data = {str(key): str(value) for key, value in loaded.items()}
        _RESPONSE_CACHE.path = path
        _RESPONSE_CACHE.data = data
        _RESPONSE_CACHE.dirty = False
        return data


def _vkey(prompt, reasoning_effort, model=None):
    digest = hashlib.sha256()
    digest.update(
        (
            resolve_model(model)
            + "\x1f"
            + str(reasoning_effort or "")
            + "\x1f"
            + prompt
        ).encode("utf-8")
    )
    return digest.hexdigest()


def cached_call_llm(
    prompt,
    api_url,
    api_key,
    timeout=None,
    retries=None,
    reasoning_effort=None,
    model=None,
    trace_id=None,
    turn_index=None,
):
    cache = _load_vcache()
    key = _vkey(prompt, reasoning_effort, model)
    with _RESPONSE_CACHE.lock:
        if key in cache:
            _RESPONSE_CACHE.hits += 1
            return cache[key], None, True
    with _RESPONSE_CACHE.lock:
        _RESPONSE_CACHE.misses += 1
        if not api_url or not api_key:
            _RESPONSE_CACHE.credential_misses += 1
            return None, "cache_miss_no_credentials", False
    response, error = call_llm(
        prompt,
        api_url,
        api_key,
        timeout=timeout,
        retries=retries,
        reasoning_effort=reasoning_effort,
        model=model,
        trace_id=trace_id,
        turn_index=turn_index,
    )
    persisted = False
    with _RESPONSE_CACHE.lock:
        if response is not None:
            cache[key] = response
            _RESPONSE_CACHE.dirty = True
            persisted = True
    if persisted:
        flush_response_cache()
    return response, error, False


def cached_call_llm_many(
    prompts,
    api_url,
    api_key,
    timeout=None,
    retries=None,
    reasoning_effort=None,
    model=None,
    max_workers=None,
    trace_ids=None,
    turn_indexes=None,
):
    if not prompts:
        return []
    workers = max(1, int(max_workers or LLM_SETTINGS.threads))
    traces = list(trace_ids or [None] * len(prompts))
    turns = list(turn_indexes or [None] * len(prompts))

    def _one(item):
        prompt, trace_id, turn_index = item
        content, error, _replayed = cached_call_llm(
            prompt,
            api_url,
            api_key,
            timeout=timeout,
            retries=retries,
            reasoning_effort=reasoning_effort,
            model=model,
            trace_id=trace_id,
            turn_index=turn_index,
        )
        return content, error

    with ThreadPoolExecutor(max_workers=min(workers, len(prompts))) as pool:
        return list(pool.map(_one, zip(prompts, traces, turns)))


def flush_response_cache():
    with _RESPONSE_CACHE.lock:
        if not _RESPONSE_CACHE.dirty or _RESPONSE_CACHE.data is None:
            return
        path = _RESPONSE_CACHE.path or response_cache_path()
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        temporary = path + ".tmp"
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(_RESPONSE_CACHE.data, handle, ensure_ascii=False)
        os.replace(temporary, path)
        _RESPONSE_CACHE.dirty = False


def response_cache_stats():
    with _RESPONSE_CACHE.lock:
        return {
            "hits": _RESPONSE_CACHE.hits,
            "misses": _RESPONSE_CACHE.misses,
            "credential_misses": _RESPONSE_CACHE.credential_misses,
            "size": (len(_RESPONSE_CACHE.data) if _RESPONSE_CACHE.data is not None else 0),
            "path": _RESPONSE_CACHE.path or response_cache_path(),
        }
