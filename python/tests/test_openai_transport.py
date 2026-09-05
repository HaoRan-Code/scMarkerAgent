from __future__ import annotations

import json

from scmarkeragent import llm_client
from scmarkeragent.config import DEFAULTS

FAKE_KEY = "unit-test-openai-key-not-secret"
FAKE_BASE = "https://openai-compatible.invalid/v1"


class _OfficialResponse:
    status_code = 200
    text = '{"object":"response"}'

    @staticmethod
    def json():
        return {
            "id": "response_test",
            "object": "response",
            "output": [
                {
                    "type": "message",
                    "role": "assistant",
                    "content": [
                        {
                            "type": "output_text",
                            "text": '{"status":"ok"}',
                        }
                    ],
                }
            ],
            "usage": {"input_tokens": 7, "output_tokens": 3},
        }


def test_official_environment_resolution(monkeypatch):
    monkeypatch.setenv("SCMA_OFFLINE", "0")
    monkeypatch.setenv("OPENAI_API_KEY", FAKE_KEY)
    monkeypatch.delenv("OPENAI_BASE_URL", raising=False)
    key, endpoint, source = llm_client.resolve_api_with_source()
    assert key == FAKE_KEY
    assert endpoint == "https://api.openai.com/v1/responses"
    assert source == "environment"
    assert llm_client.provider_mode() == "official_openai"


def test_optional_compatible_base(monkeypatch):
    monkeypatch.setenv("SCMA_OFFLINE", "0")
    monkeypatch.setenv("OPENAI_API_KEY", FAKE_KEY)
    monkeypatch.setenv("OPENAI_BASE_URL", FAKE_BASE)
    key, endpoint, source = llm_client.resolve_api_with_source()
    assert key == FAKE_KEY
    assert endpoint == FAKE_BASE + "/responses"
    assert source == "environment"
    assert llm_client.provider_mode() == "custom_openai_base"


def test_offline_mode_forces_empty_credentials(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", FAKE_KEY)
    monkeypatch.setenv("OPENAI_BASE_URL", FAKE_BASE)
    monkeypatch.setenv("SCMA_OFFLINE", "1")
    assert llm_client.resolve_api_with_source() == (
        "",
        "",
        "offline",
    )


def test_official_responses_payload_and_audit(tmp_path, monkeypatch):
    calls = []
    log = tmp_path / "cold.jsonl"
    monkeypatch.setenv("SCMA_LLM_RAW_LOG", str(log))
    monkeypatch.setenv("OPENAI_API_KEY", FAKE_KEY)
    monkeypatch.delenv("OPENAI_BASE_URL", raising=False)

    def post(url, **kwargs):
        calls.append({"url": url, **kwargs})
        return _OfficialResponse()

    monkeypatch.setattr(llm_client.requests, "post", post)
    content, error = llm_client.call_llm(
        "audit prompt",
        "https://api.openai.com/v1/responses",
        FAKE_KEY,
        retries=1,
        reasoning_effort="xhigh",
    )
    assert content == '{"status":"ok"}'
    assert error is None
    payload = calls[0]["json"]
    assert payload == {
        "model": llm_client.MODEL,
        "input": "audit prompt",
        "reasoning": {"effort": "xhigh"},
    }
    forbidden = {"tem" + "perature", "se" + "ed"}
    assert forbidden.isdisjoint(payload)
    record = json.loads(log.read_text())
    assert record["request"] == payload
    assert record["parsed_json"] == {"status": "ok"}
    assert record["provider_response"]["object"] == "response"
    assert "Authorization" not in record
    assert FAKE_KEY not in log.read_text()


def test_transport_controls_are_shared_config_only():
    assert isinstance(DEFAULTS["llm"]["timeout_seconds"], int)
    assert isinstance(DEFAULTS["llm"]["retries"], int)
    assert isinstance(DEFAULTS["llm"]["threads"], int)


def test_a_prompt_keys_to_the_same_digest_on_both_arms(monkeypatch):
    """The same literals tests/test_r_llm_client.R asserts.

    Each arm computes the key from (model, effort, prompt) rather than reading the
    other's; pinning both to one hex string is what makes a run finished on one arm
    replay on the other instead of paying for every prompt twice. The last case is why
    a literal is needed at all: it fixes what "no effort named" hashes as, which no
    amount of agreement about the algorithm would settle.
    """
    monkeypatch.setattr(llm_client, "MODEL_OVERRIDE", "")
    model = "unit-test-model"
    assert (
        llm_client._vkey("ping", "high", model)
        == "d738f945517ec7cb7c8ffe1e671d929b9637e3d59b22e08e8514114bb26c896a"
    )
    assert (
        llm_client._vkey("β cell — 中文", "high", model)
        == "fb18c245d3bf92d71d00e774a3de264bf3fd26092adacee5166d46aa5435375f"
    )
    assert (
        llm_client._vkey("ping", None, model)
        == "16991eeca85954e48c4c92431a3ca4c082f0f6fe735fddc425a311ba479859a2"
    )
