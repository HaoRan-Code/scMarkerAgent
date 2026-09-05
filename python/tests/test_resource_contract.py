from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scmarkeragent import config

ROOT = Path(__file__).resolve().parents[1]


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def test_resource_manifest_matches_resolved_resources():
    manifest = json.loads((ROOT / "resource_manifest.json").read_text())
    resolved = {
        "markers": Path(config.DB_CSV),
        "sources": Path(config.DB_SOURCES),
        "uberon_ontology": Path(config.OBO_UBERON),
        "ortholog_human_mouse": Path(config.ORTHO_DIR) / "ortho_Human_to_Mouse.csv",
        "ortholog_human_rat": Path(config.ORTHO_DIR) / "ortho_Human_to_Rat.csv",
    }
    if not resolved["markers"].exists():
        # The curated bundle is distributed separately and is absent from a plain
        # checkout or install; there is then nothing to hold the manifest against.
        # `download-resources --check-only` is the verification that applies there.
        pytest.skip(
            "curated resource bundle not present; set SCMA_RESOURCE_DIR to a "
            "downloaded bundle to verify it against the manifest"
        )
    for name, path in resolved.items():
        assert path.exists(), name
        record = manifest["files"][name]
        assert path.stat().st_size == record["size"]
        assert _sha(path) == record["sha256"]
    assert set(manifest["files"]) == set(resolved)


def test_no_cell_ontology_resource_is_declared_or_resolved():
    """The package must not depend on cl.obo in any form.

    A curated CL id reaches the output as a display column read off the marker resource;
    nothing resolves, expands or compares one, so declaring the ontology file would
    reintroduce a dependency the pipeline no longer has.
    """
    manifest = json.loads((ROOT / "resource_manifest.json").read_text())
    assert "cell_ontology" not in manifest["files"]
    assert "cell_ontology" not in config.DEFAULTS["resources"]
    assert not hasattr(config, "OBO")
    assert not (ROOT / "src/scmarkeragent/resources/static/ontology/cl.obo").exists()


def test_curated_database_is_the_only_execution_resource():
    # Pin the invariant, not the release: the shipped config and the shipped
    # manifest must name the same curated bundle, whichever rebuild is current.
    manifest = json.loads((ROOT / "resource_manifest.json").read_text())
    assert config.DEFAULTS["db_version"] == manifest["resource_version"]
    source = (ROOT / "src/scmarkeragent/marker_database.py").read_text(encoding="utf-8")
    assert "DB_" + "CURATION" not in source
    assert "load_curation_" + "blocklist" not in source
