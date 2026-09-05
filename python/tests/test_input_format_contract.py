"""The input format, and nothing else, selects the execution arm."""

from __future__ import annotations

import pytest

from scmarkeragent import input_format


def test_extension_selects_the_arm():
    assert input_format.arm_for_input("dataset.h5ad") == input_format.PYTHON_ARM
    assert input_format.arm_for_input("dataset.rds") == input_format.R_ARM


def test_selection_is_case_insensitive_and_path_independent():
    assert input_format.arm_for_input("/tmp/a/DATASET.H5AD") == input_format.PYTHON_ARM
    assert input_format.arm_for_input("/tmp/a/Dataset.Rds") == input_format.R_ARM


def test_an_unsupported_format_is_refused_rather_than_guessed():
    # Routing an unknown format to either arm would fail deep inside a reader, or
    # silently succeed against a partially parsed object.
    with pytest.raises(ValueError) as error:
        input_format.arm_for_input("dataset.loom")
    assert ".loom" in str(error.value)
    for expected in (".h5ad", ".rds"):
        assert expected in str(error.value)


def test_every_arm_has_exactly_one_native_format():
    # Neither format may be readable by both arms: the two arms read their own object
    # directly, so a shared format would mean one of them converts.
    arms = list(input_format.ARM_BY_SUFFIX.values())
    assert sorted(arms) == sorted({input_format.PYTHON_ARM, input_format.R_ARM})


def test_check_input_routes_a_seurat_object_to_the_other_arm(tmp_path):
    rds = tmp_path / "dataset.rds"
    rds.write_bytes(b"")
    with pytest.raises(ValueError) as error:
        input_format.check_input(rds, input_format.PYTHON_ARM)
    message = str(error.value)
    assert "R arm" in message
    # The message has to name the way out, not just the problem.
    assert input_format.ENTRY_POINT[input_format.R_ARM] in message


def test_check_input_refuses_an_unsupported_format_before_any_stage(tmp_path):
    loom = tmp_path / "dataset.loom"
    loom.write_bytes(b"")
    with pytest.raises(ValueError) as error:
        input_format.check_input(loom, input_format.PYTHON_ARM)
    message = str(error.value)
    assert ".loom" in message
    for expected in (".h5ad", ".rds"):
        assert expected in message


def test_check_input_reports_a_missing_file_as_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        input_format.check_input(tmp_path / "absent.h5ad", input_format.PYTHON_ARM)


def test_check_input_returns_the_resolved_path(tmp_path):
    h5ad = tmp_path / "dataset.h5ad"
    h5ad.write_bytes(b"")
    resolved = input_format.check_input(h5ad, input_format.PYTHON_ARM)
    assert resolved == h5ad.resolve()
