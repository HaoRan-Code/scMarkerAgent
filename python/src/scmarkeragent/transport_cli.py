"""Language-neutral JSON adapter for the shared LLM transport."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import llm_client

SCHEMA_VERSION = "scmarkeragent-transport-request-v1"
ALLOWED_EFFORT = {"low", "medium", "high", "xhigh"}


@dataclass(frozen=True)
class TransportRequest:
    prompt: str
    reasoning_effort: str

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "TransportRequest":
        if set(value) != {
            "schema_version",
            "prompt",
            "reasoning_effort",
        }:
            raise ValueError("invalid transport request fields")
        if value["schema_version"] != SCHEMA_VERSION:
            raise ValueError("unsupported transport request schema")
        prompt = value["prompt"]
        effort = value["reasoning_effort"]
        if not isinstance(prompt, str) or not prompt.strip():
            raise ValueError("transport prompt must be non-empty")
        if effort not in ALLOWED_EFFORT:
            raise ValueError("invalid reasoning effort")
        return cls(prompt=prompt, reasoning_effort=effort)


def execute(request: TransportRequest) -> dict[str, Any]:
    api_key, api_url = llm_client.resolve_api()
    content, error, cached = llm_client.cached_call_llm(
        request.prompt,
        api_url,
        api_key,
        reasoning_effort=request.reasoning_effort,
    )
    if content is None and not api_key:
        return {
            "schema_version": SCHEMA_VERSION,
            "status": "skipped_no_credentials",
            "content": None,
            "error": "cache_miss_no_credentials",
            "cached": False,
            "response_sha256": "",
        }
    if content is None:
        raise llm_client.LlmTransportError(error or "LLM transport failed")
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "replayed" if cached else "cold_call",
        "content": content,
        "error": None,
        "cached": cached,
        "response_sha256": hashlib.sha256(content.encode("utf-8")).hexdigest(),
    }


def write_atomic(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="scmarkeragent-transport")
    parser.add_argument("--request", required=True)
    parser.add_argument("--response", required=True)
    args = parser.parse_args(argv)
    request_path = Path(args.request)
    response_path = Path(args.response)
    request = TransportRequest.from_dict(
        json.loads(request_path.read_text(encoding="utf-8"))
    )
    write_atomic(response_path, execute(request))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
