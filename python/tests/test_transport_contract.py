from __future__ import annotations


import pytest

from scmarkeragent.transport_cli import (
    SCHEMA_VERSION,
    TransportRequest,
    execute,
)


def test_transport_request_schema():
    request = TransportRequest.from_dict(
        {
            "schema_version": SCHEMA_VERSION,
            "prompt": "prompt",
            "reasoning_effort": "xhigh",
        }
    )
    assert request.prompt == "prompt"


def test_transport_request_rejects_extra_fields():
    with pytest.raises(ValueError):
        TransportRequest.from_dict(
            {
                "schema_version": SCHEMA_VERSION,
                "prompt": "prompt",
                "reasoning_effort": "xhigh",
                "unsupported": True,
            }
        )


def test_transport_no_key_status(monkeypatch):
    monkeypatch.setenv("SCMA_OFFLINE", "1")
    result = execute(TransportRequest(prompt="prompt", reasoning_effort="xhigh"))
    assert result["status"] == "skipped_no_credentials"
    assert result["content"] is None
