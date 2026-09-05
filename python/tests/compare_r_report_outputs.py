#!/usr/bin/env python3
"""Compare an R-native report package against the Python golden, field by field.

The acceptance bar for T8 is NOT byte-identity. Two things make byte-identity the wrong
test rather than merely a strict one:

* R's `as.numeric` is not correctly rounded -- it lands 1 ULP away from Python's `float`
  on about 0.014% of real values -- so "the shortest string that reads back" can legally
  differ between the arms while the double is the same. Comparing the STRINGS would fail
  on a difference no reader can see; comparing the DOUBLES catches every difference that
  matters and none that does not.
* Figures cannot match: matplotlib and R graphics do not produce the same bytes, ever.
  What can be asserted about them is the file set.

So: text cells must match exactly, numeric cells must be the same double, JSON key order
is free, figures are compared by name, and the ZIP is compared by member path. Anything
that records WHEN or WHERE the run happened is normalised away -- it says nothing about
whether the two arms computed the same thing.

  Usage:
    python tests/compare_r_report_outputs.py \\
      --expected validation/synthetic_golden/results \\
      --actual   /tmp/t8_native/results/synthetic_r
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import zipfile
from pathlib import Path
from typing import Any

csv.field_size_limit(10**9)

# Says when or where the run happened, not what it computed.
VOLATILE_KEYS = {"generated_at_utc", "generated_at", "resolved_paths", "runtime"}
# Written by the orchestrator, not by the report layer under test.
SKIP_FILES = {"run_manifest.json", "resolved_config.json"}
BINARY_SUFFIXES = {".pdf", ".png", ".jpg", ".jpeg", ".svg", ".zip"}


def as_number(text: str) -> float | None:
    """The double this cell denotes, or None if it does not denote one."""
    stripped = text.strip()
    if not stripped or stripped.lower() in {"n/a", "na", "nan", "none", "-", ""}:
        return None
    try:
        value = float(stripped)
    except ValueError:
        return None
    return value


def cells_agree(left: str, right: str) -> bool:
    if left == right:
        return True
    a, b = as_number(left), as_number(right)
    if a is None or b is None:
        return False
    if math.isnan(a) and math.isnan(b):
        return True
    return a == b


def read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def compare_csv(name: str, expected: Path, actual: Path, findings: list[str]) -> None:
    exp_header, exp_rows = read_table(expected)
    act_header, act_rows = read_table(actual)
    if exp_header != act_header:
        findings.append(
            f"{name}: header differs\n    expected {exp_header}\n    actual   {act_header}"
        )
        return
    if len(exp_rows) != len(act_rows):
        findings.append(f"{name}: {len(exp_rows)} data rows expected, {len(act_rows)} found")
        return
    mismatches = 0
    for index, (exp_row, act_row) in enumerate(zip(exp_rows, act_rows)):
        for column, (left, right) in enumerate(zip(exp_row, act_row)):
            if cells_agree(left, right):
                continue
            mismatches += 1
            if mismatches <= 5:
                findings.append(
                    f"{name}: row {index + 2}, column '{exp_header[column]}'\n"
                    f"    expected {left!r}\n    actual   {right!r}"
                )
    if mismatches > 5:
        findings.append(f"{name}: ...and {mismatches - 5} further differing cells")


def normalise(value: Any) -> Any:
    """Drop the keys that record when or where, and make numbers comparable as numbers."""
    if isinstance(value, dict):
        return {k: normalise(v) for k, v in value.items() if k not in VOLATILE_KEYS}
    if isinstance(value, list):
        return [normalise(item) for item in value]
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value)
    return value


def json_difference(expected: Any, actual: Any, path: str = "") -> list[str]:
    if isinstance(expected, dict) and isinstance(actual, dict):
        out = []
        for key in sorted(set(expected) | set(actual)):
            if key not in expected:
                out.append(f"{path}.{key}: only in actual")
            elif key not in actual:
                out.append(f"{path}.{key}: only in expected")
            else:
                out += json_difference(expected[key], actual[key], f"{path}.{key}")
        return out
    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            return [f"{path}: {len(expected)} items expected, {len(actual)} found"]
        out = []
        for index, (left, right) in enumerate(zip(expected, actual)):
            out += json_difference(left, right, f"{path}[{index}]")
        return out
    if isinstance(expected, float) and isinstance(actual, float):
        if expected == actual or (math.isnan(expected) and math.isnan(actual)):
            return []
        return [f"{path}: {expected!r} expected, {actual!r} found"]
    if expected != actual:
        return [f"{path}: {expected!r} expected, {actual!r} found"]
    return []


def compare_jsonl(name: str, expected: Path, actual: Path, findings: list[str]) -> None:
    def load(path: Path) -> list[Any]:
        with path.open(encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]

    exp, act = load(expected), load(actual)
    if len(exp) != len(act):
        findings.append(f"{name}: {len(exp)} records expected, {len(act)} found")
        return
    shown = 0
    for index, (left, right) in enumerate(zip(exp, act)):
        for difference in json_difference(normalise(left), normalise(right), f"record[{index}]"):
            shown += 1
            if shown <= 8:
                findings.append(f"{name}: {difference}")
    if shown > 8:
        findings.append(f"{name}: ...and {shown - 8} further differences")


def compare_zip(name: str, expected: Path, actual: Path, findings: list[str]) -> None:
    """Member paths only: the archive's timestamps and compression are not the contract.

    Absolute or escaping paths are checked here rather than in a separate test because a
    ZIP is the one artifact a user unpacks, and this is where its member list is already
    in hand.
    """
    with zipfile.ZipFile(expected) as handle:
        exp_names = sorted(handle.namelist())
    with zipfile.ZipFile(actual) as handle:
        act_names = sorted(handle.namelist())
    escaping = [
        member
        for member in act_names
        if member.startswith("/") or ".." in Path(member).parts
    ]
    if escaping:
        findings.append(f"{name}: ZIP members are not relative: {escaping[:5]}")
    if exp_names != act_names:
        only_expected = sorted(set(exp_names) - set(act_names))
        only_actual = sorted(set(act_names) - set(exp_names))
        findings.append(
            f"{name}: member set differs\n"
            f"    missing from actual: {only_expected[:8]}\n"
            f"    unexpected: {only_actual[:8]}"
        )


def relative_files(root: Path) -> set[str]:
    return {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.name not in SKIP_FILES
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", required=True, type=Path)
    parser.add_argument("--actual", required=True, type=Path)
    parser.add_argument("--report", type=Path, help="write a JSON verdict here")
    args = parser.parse_args()

    findings: list[str] = []
    exp_files = relative_files(args.expected)
    act_files = relative_files(args.actual)
    missing = sorted(exp_files - act_files)
    extra = sorted(act_files - exp_files)
    if missing:
        findings.append(f"files missing from actual: {missing}")
    if extra:
        findings.append(f"files not in the golden: {extra}")

    for name in sorted(exp_files & act_files):
        expected, actual = args.expected / name, args.actual / name
        suffix = expected.suffix.lower()
        if suffix == ".csv":
            compare_csv(name, expected, actual, findings)
        elif suffix == ".jsonl":
            compare_jsonl(name, expected, actual, findings)
        elif suffix == ".zip":
            compare_zip(name, expected, actual, findings)
        elif suffix == ".json":
            with expected.open(encoding="utf-8") as handle:
                left = normalise(json.load(handle))
            with actual.open(encoding="utf-8") as handle:
                right = normalise(json.load(handle))
            for difference in json_difference(left, right, name)[:8]:
                findings.append(difference)
        elif suffix in BINARY_SUFFIXES:
            # Present and non-empty is the whole claim: two graphics stacks never agree
            # on bytes, and a zero-byte figure is the failure worth catching.
            if actual.stat().st_size == 0:
                findings.append(f"{name}: written but empty")
        else:
            if expected.read_text(encoding="utf-8") != actual.read_text(encoding="utf-8"):
                findings.append(f"{name}: text differs")

    verdict = {
        "schema_version": "r-report-parity-v1",
        "expected": str(args.expected),
        "actual": str(args.actual),
        "files_compared": len(exp_files & act_files),
        "passed": not findings,
        "findings": findings,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(verdict, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

    print(f"compared {verdict['files_compared']} files")
    if findings:
        print(f"FAILED with {len(findings)} finding(s):")
        for finding in findings:
            print(f"  - {finding}")
        return 1
    print("PASSED: every field agrees")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
