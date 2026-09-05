from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = PACKAGE_ROOT / "config" / "defaults.json"
OUTPUT_SCHEMA_PATH = PACKAGE_ROOT / "schemas" / "output_schema.json"


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


CONFIG_PATH = Path(os.environ.get("SCMA_CONFIG", DEFAULT_CONFIG_PATH)).resolve()
DEFAULTS = _read_json(CONFIG_PATH)
if DEFAULTS.get("schema_version") != "scmarkeragent-defaults-v1":
    raise ValueError("unsupported scMarkerAgent defaults schema")
OUTPUT_SCHEMA = _read_json(OUTPUT_SCHEMA_PATH)
if OUTPUT_SCHEMA.get("schema_version") != "scmarkeragent-output-v1":
    raise ValueError("unsupported scMarkerAgent output schema")
for section in ("cluster_summary", "marker_evidence"):
    columns = OUTPUT_SCHEMA.get(section, {}).get("columns")
    if not isinstance(columns, list) or len(columns) != len(set(columns)):
        raise ValueError(f"invalid output columns for {section}")


@dataclass(frozen=True)
class LlmSettings:
    """Validated shared controls for Python and R LLM transports."""

    model: str
    threads: int
    timeout_seconds: int
    retries: int


_llm = DEFAULTS.get("llm")
_agents = DEFAULTS.get("cluster_annotation")
if not isinstance(_llm, dict) or not isinstance(_agents, dict):
    raise ValueError("defaults must define llm and cluster_annotation objects")
LLM_SETTINGS = LlmSettings(
    model=str(_llm["model"]),
    # Runtime-only transport controls. These do not change annotation semantics or
    # prompt contents; they let benchmark/release validation reduce concurrency or
    # extend a slow provider timeout without editing the frozen scientific config.
    threads=int(os.environ.get("SCMA_LLM_THREADS", _llm["threads"])),
    timeout_seconds=int(
        os.environ.get("SCMA_LLM_TIMEOUT_SECONDS", _llm["timeout_seconds"])
    ),
    retries=int(os.environ.get("SCMA_LLM_RETRIES", _llm["retries"])),
)
if (
    LLM_SETTINGS.threads < 1
    or LLM_SETTINGS.timeout_seconds < 1
    or LLM_SETTINGS.retries < 1
):
    raise ValueError("LLM controls must be positive integers")


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name, "").strip()
    return Path(value).expanduser().resolve() if value else default.resolve()


RESOURCE_DIR = _env_path("SCMA_RESOURCE_DIR", PACKAGE_ROOT / "resources" / "static")
STATIC_RESOURCE_DIR = (PACKAGE_ROOT / "resources" / "static").resolve()
WORK_DIR = _env_path("SCMA_WORK_DIR", Path.cwd() / ".scmarkeragent")
CACHE_DIR = _env_path("SCMA_CACHE", WORK_DIR / "cache")
RESULTS_DIR = _env_path("SCMA_RESULTS", WORK_DIR / "results")


def resource_path(key: str) -> str:
    relative = Path(DEFAULTS["resources"][key])
    candidate = RESOURCE_DIR / relative
    if candidate.exists():
        return str(candidate)
    fallback = STATIC_RESOURCE_DIR / relative
    return str(fallback)


DB_CSV = resource_path("markers")
DB_SOURCES = resource_path("sources")
# There is deliberately no Cell Ontology resource. Identity is the curated free-text
# cell-type string; a CL id is carried through to the output only as a display column
# read off the marker resource, and never resolved, expanded or compared.
OBO_UBERON = resource_path("uberon_ontology")
ORTHO_DIR = str(
    (RESOURCE_DIR / "ortholog")
    if (RESOURCE_DIR / "ortholog").exists()
    else (STATIC_RESOURCE_DIR / "ortholog")
)
PROMPT_DIR = str(PACKAGE_ROOT / "prompts")

QC_MIN_GENES = int(DEFAULTS["preprocessing"]["qc_min_genes"])
QC_MIN_CELLS = int(DEFAULTS["preprocessing"]["qc_min_cells"])
QC_MAX_MT = float(DEFAULTS["preprocessing"]["qc_max_mt_percent"])
LEIDEN_RESOLUTION = float(DEFAULTS["preprocessing"]["leiden_resolution"])
COMPUTE_UMAP = bool(DEFAULTS["features"]["compute_umap"])
EXCLUDE_IN_VITRO = bool(DEFAULTS["features"]["exclude_in_vitro"])
# Both arms read these from the shared configuration; neither hard-codes them, so a
# configuration sweep moves the Python and R arms together.
CROSS_SPECIES = bool(DEFAULTS["features"]["cross_species_markers"])
CORROBORATION_ONLY = bool(DEFAULTS["retrieval"]["corroboration_only"])
TISSUE_ROOT = str(DEFAULTS["eligibility"]["tissue_root"])
RU_MIN_ELIG_GENES = int(DEFAULTS["retrieval"]["min_eligible_genes"])
MIN_CORROBORATING_PUBLICATIONS = int(
    DEFAULTS["eligibility"]["min_corroborating_publications"]
)
CORROBORATING_TIERS = frozenset(
    str(value) for value in DEFAULTS["eligibility"]["corroborating_tiers"]
)
CLUSTER_ANNOTATION_ENABLED = bool(DEFAULTS["features"]["cluster_annotation"])
# Every table writes this token wherever a cell does not apply, so a blank never has to
# be read as either "absent" or "not applicable".
NOT_AVAILABLE = str(DEFAULTS["output"]["not_applicable_token"])
FIGURES_DIR = "figures"
FIGDATA_DIR = "figure_data"
TABLE_A_FILE = str(DEFAULTS["output"]["table_a_file"])
TABLE_B_FILE = str(DEFAULTS["output"]["table_b_file"])
METRICS_FILE = "metrics.csv"
STATS_FILE = "stats_reproducibility.txt"
CLUSTER_EVIDENCE_FILE = str(OUTPUT_SCHEMA["cluster_evidence"]["file"])


def na_display(value: Any) -> str:
    """Reader-facing text for one table cell, with the shared not-applicable token.

    Empty, missing and NaN all mean the same thing to a reader -- the field does not
    apply to this row -- so every table renders them as one explicit token instead of a
    blank that could equally be read as a lost value.
    """
    if value is None:
        return NOT_AVAILABLE
    if isinstance(value, float) and value != value:
        return NOT_AVAILABLE
    text = str(value).strip()
    return text if text and text.lower() not in {"nan", "none", "na"} else NOT_AVAILABLE


def ensure_workdirs() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def resolved_config(overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    from .llm_client import credential_source, provider_mode

    resolved_credential_source = credential_source()
    resolved_provider_mode = provider_mode()
    value = json.loads(json.dumps(DEFAULTS))
    value["resolved_paths"] = {
        "package_root": str(PACKAGE_ROOT),
        "resource_dir": str(RESOURCE_DIR),
        "cache_dir": str(CACHE_DIR),
        "results_dir": str(RESULTS_DIR),
        "markers": DB_CSV,
        "sources": DB_SOURCES,
        "uberon_ontology": OBO_UBERON,
        "ortholog_dir": ORTHO_DIR,
        "prompt_dir": PROMPT_DIR,
    }
    value["runtime"] = {
        "llm_model": os.environ.get("SCMA_LLM_MODEL", str(DEFAULTS["llm"]["model"])),
        # One agent decides a cluster, so one model name is recorded. `SCMA_LLM_MODEL`,
        # when set, pins it and is reflected here.
        "annotator_model": os.environ.get(
            "SCMA_LLM_MODEL", str(_agents["annotator_model"])
        ),
        "credential_source": resolved_credential_source,
        "provider_mode": resolved_provider_mode,
        "offline": os.environ.get("SCMA_OFFLINE", "0").lower()
        in {"1", "true", "yes", "on"},
    }
    if overrides:
        value["cli_overrides"] = overrides
    return value


def write_resolved_config(path: str | Path, overrides=None) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(resolved_config(overrides), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def summary() -> str:
    return json.dumps(resolved_config(), indent=2, ensure_ascii=False)


if __name__ == "__main__":
    print(summary())
