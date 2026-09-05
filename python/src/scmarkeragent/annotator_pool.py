#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The measured evidence one annotating agent reasons over, and the tools it may query.

Everything here is a lookup over artifacts the pipeline already produced -- the curated
marker panels, the one-versus-rest DE table, and the curated source sentences. Nothing is
re-decided, nothing is scored, and no threshold removes a candidate: the agent compares
identities, and the code only supplies facts and, afterwards, checks the numbers it
quoted.

Two properties of the opening packet carry the design.

Each candidate's panel arrives WHOLE, ordered by curated publication support. Showing
only the up-regulated markers hides the two facts that decide a call. A claimed identity
whose best-corroborated markers are measured near zero is refuted by exactly those rows
-- an enteroendocrine claim whose CHGA (160 publications) sits at 2.3% in versus 4.1%
out is not visible from the handful of incidental genes that did rise. And a cell type
that dominates the dataset has its own defining markers driven to a flat fold change by
the one-versus-rest comparison, so an enterocyte cluster reads as having no enterocyte
evidence at all unless FABP2 at 92.4% in versus 68.5% out is on the page.

Publication support is the ordering, not a gate. It is the curated database's own record
of which gene defines the cell type, and sorting on it puts that gene first; whether the
cluster raises it is a separate question the measurements answer on the same row.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from . import evidence_gate
from .config import DEFAULTS, NOT_AVAILABLE

_A = DEFAULTS["cluster_annotation"]
SOURCES_PER_GENE = int(_A["sources_per_marker"])
SOURCE_BATCHES_PER_MARKER = int(_A["source_batches_per_marker"])
SOURCE_GENES_PER_QUERY = int(_A["source_genes_per_query"])
SOURCE_MAX_CHARS = int(_A["source_sentence_max_chars"])
UNMEASURED_SHOWN = int(_A["unmeasured_genes_shown"])
WARNING_GENES_SHOWN = int(_A["warning_genes_shown"])
SUBTYPE_EXCLUSIVE_REQUIRED = int(_A["subtype_exclusive_markers_required"])
# An exclusion the cluster actually detects is the one curated claim whose weight turns
# entirely on what the literature says, and it is not the kind of row an agent thinks to
# ask about. Those sentences ship with the packet; every other source is still fetched by
# name. The threshold selects which rows carry provenance and nothing else -- every
# negative marker appears in the panel either way, and the reading stays with the agent.
NEGATIVE_SOURCE_MIN_PCT_IN = float(_A["negative_source_min_pct_in"])
# How many rows of the cluster's own ranked evidence travel with the packet. It is a
# reading order, not a shortlist: every one of these rows is already on some candidate's
# panel, and the panels remain the complete record.
CLUSTER_MARKER_ROWS = int(_A["cluster_marker_rows"])
# How many of a judged label's genes arrive with their sentences already attached. It caps
# the SENTENCES, not the panel: the panel is whole, because a truncated panel let 10% of
# the genes an answer was actually built on -- C1QB for Kupffer cell, VSIG4 for alveolar
# macrophage -- reach the judge nowhere at all, and a check cannot weigh a row it was
# never shown. Sentences are the expensive half, so the best-published ones ship and the
# judge asks for the rest.
JUDGE_SOURCES_PUSHED = int(_A["judge_sources_pushed"])
JUDGE_SOURCES_PER_GENE = int(_A["judge_sources_per_gene"])
# How many cell types one `pool_search` answer names. Enough to show whether the retrieved
# menu is missing an identity, short enough that it cannot become a second candidate pool.
POOL_SEARCH_ROWS = int(_A["pool_search_rows"])

# The tools the annotator may call. `sources` is the one piece of evidence deliberately
# kept out of the opening packet -- the curated sentences for one dataset run to
# thousands of rows, and only a handful ever decide a cluster. The others answer
# questions the packet cannot: how a gene behaves in the OTHER clusters, which candidates
# claim a gene, and -- since retrieval hands over 15 of a median 67 admitted cell types
# and nothing downstream could see past that -- which cell types OUTSIDE the menu claim a
# gene this cluster raises.
TOOLS = ("sources", "gene_across_clusters", "candidates_with_gene", "pool_search")

# The tools a judge may call, and deliberately fewer. `candidates_with_gene` and
# `pool_search` both name identities other than the one under test, which is exactly what
# a judge must not be given: handed a menu it could prefer a different label, and that
# disagreement says nothing about whether THIS label is carried by the measurements.
# What it may ask for is more of its own label's sentences, and how a gene behaves in the
# other clusters -- neither of which names another identity.
JUDGE_TOOLS = ("sources", "gene_across_clusters")

# Detection fractions are carried on the DE table's own scale, percent, everywhere the
# agent sees them and everywhere the audit compares them. One scale, one place.
PCT_DECIMALS = 1
LOGFC_DECIMALS = 3
SPECIFICITY_DECIMALS = 3

# What the agent is shown of a marker row. `tier` is omitted because it carries nothing
# `n_pub` does not: the resource derives it from the publication count by a fixed rule
# (1 -> low, 2 -> medium, >=3 -> high), with no exception anywhere in a panel, and the
# packet's own field list never described it. The row keeps it for `marker_evidence.csv`,
# where a reader asked for a tier rather than a count.
PACKET_MARKER_FIELDS = (
    "gene",
    "polarity",
    "n_pub",
    "specificity",
    "pct_in",
    "pct_out",
    "avg_log2FC",
    "significant",
)


def _num(value, decimals):
    if value is None:
        return None
    number = float(value)
    return None if not np.isfinite(number) else round(number, decimals)


def _clip(text: str, limit: int) -> str:
    value = " ".join(str(text or "").split())
    return value if len(value) <= limit else value[: limit - 1] + "\u2026"


def _identifier(value) -> str:
    text = str(value or "").strip()
    return text if text and text.lower() not in {"nan", "none", "na"} else NOT_AVAILABLE


def is_raised(marker: dict) -> bool:
    """Whether the cluster raises this marker: up-regulated one-versus-rest.

    One predicate, one meaning, used by the packet's own wording and by the post-hoc
    audit. It is deliberately not a defensibility test -- a gene can be raised and
    uninformative, or unraised and decisive -- so nothing here removes a candidate.
    """
    value = marker.get("avg_log2FC")
    return value is not None and float(value) > 0.0


def _significant(measured: dict) -> bool:
    """The retrieval gate's verdict on this row, carried onto the row itself.

    `evidence_gate.sig_pass` is the predicate that decided whether this marker could put
    its candidate in the pool at all. Reusing it here rather than restating the
    thresholds is what keeps the packet from drifting away from retrieval a second time.
    """
    pct_in = measured.get("pct_in")
    padj = measured.get("adjusted_p_value")
    if measured.get("avg_log2FC") is None or pct_in is None or padj is None:
        return False
    return bool(
        evidence_gate.sig_pass(
            float(measured["avg_log2FC"]), float(pct_in) / 100.0, float(padj)
        )
    )


# ------------------------------------------------------------------- building ----
def build_pool(
    scoring: dict, prep: pd.DataFrame | dict, sources=None
) -> dict[str, Any]:
    """Every cluster's candidates with their complete measured panels.

    `prep` is the preprocessing payload; only its DE table is read. The result is a plain
    nested structure so the R arm can build the identical thing and the two can be
    compared field for field.

    `sources` supplies the curated sentences behind detected exclusions; without it those
    rows simply carry no provenance and everything else is unchanged.
    """
    de = (prep["de"] if isinstance(prep, dict) else prep).copy()
    de["group"] = de["group"].astype(str)
    de["gene_key"] = de["feature"].astype(str).str.upper()

    stats: dict[tuple[str, str], dict[str, float]] = {}
    by_gene: dict[str, list[tuple[str, float, float]]] = {}
    native: dict[str, str] = {}
    for row in de.itertuples(index=False):
        cluster, gene = str(row.group), str(row.gene_key)
        pct_in = _num(row.pct_in, PCT_DECIMALS)
        pct_out = _num(row.pct_out, PCT_DECIMALS)
        stats[(cluster, gene)] = {
            "pct_in": pct_in,
            "pct_out": pct_out,
            "avg_log2FC": _num(row.avg_log2FC, LOGFC_DECIMALS),
            "adjusted_p_value": _num(row.padj, 6),
        }
        by_gene.setdefault(gene, []).append((cluster, pct_in, pct_out))
        native.setdefault(gene, str(row.feature))

    panel_records = scoring["panel_records"]
    meta: dict[tuple[str, str], dict[str, Any]] = {}
    carriers: dict[str, list[tuple[str, int, str, str]]] = {}
    curated: dict[str, set[str]] = {}
    for row in panel_records.itertuples(index=False):
        candidate, gene = str(row.candidate), str(row.gene_key).upper()
        polarity = str(row.marker_polarity)
        # `pd.NA or x` raises rather than falling through, and a marker resource is not
        # required to carry either column: confidence_tier is empty in every database
        # except our own. The R arm already reads them as missing-tolerant, so guarding
        # here is also what keeps the two arms reading one resource the same way.
        n_pub = int(row.n_pub) if pd.notna(row.n_pub) else 0
        tier = str(row.tier) if pd.notna(row.tier) else NOT_AVAILABLE
        meta[(candidate, gene)] = {
            "n_pub": n_pub,
            "tier": tier,
            "polarity": polarity,
        }
        carriers.setdefault(gene, []).append((candidate, n_pub, tier, polarity))
        if polarity == "positive":
            curated.setdefault(candidate, set()).add(gene)

    context = scoring["context"]
    scored = scoring["scored"]
    positive_panels = scoring["measured_panels"]
    negative_panels = scoring["negative_panels"]
    specificity = scoring.get("marker_specificity") or {}

    clusters: dict[str, dict[str, Any]] = {}
    for cluster, state in scoring["clusters"].items():
        cluster = str(cluster)
        frame = scored[scored["cluster"] == cluster]
        candidates = []
        for name in state.get("candidates") or []:
            name = str(name)
            match = frame[frame["candidate"] == name]
            if match.empty:
                continue
            row = match.iloc[0]
            markers = _panel_rows(
                name,
                positive_panels.get(name, []),
                negative_panels.get(name, []),
                cluster,
                stats,
                meta,
                native,
                specificity,
            )
            measured = {marker["gene"].upper() for marker in markers}
            unmeasured = sorted(curated.get(name, set()) - measured)
            candidates.append(
                {
                    "cell_type": name,
                    "retrieval_rank": int(row["retrieval_rank"]),
                    "markers": markers,
                    "exclusion_sources": _exclusion_sources(name, markers, sources),
                    "unmeasured_curated_genes": {
                        "count": len(unmeasured),
                        "genes": [native.get(gene, gene) for gene in unmeasured][
                            :UNMEASURED_SHOWN
                        ],
                    },
                    "program": {
                        "median_in": _num(row["program_in_median"], 4),
                        "median_out": _num(row["program_out_median"], 4),
                    },
                }
            )
        clusters[cluster] = {
            "cluster_id": cluster,
            "status": state.get("status"),
            "n_cells": int(state.get("n_cells", 0)),
            "candidates": candidates,
        }

    return {
        "context": {
            "species": context["species"],
            "tissue": context["tissue"],
            "disease": context["disease"],
            "development_stage": context["development_stage"] or NOT_AVAILABLE,
            "clusters_in_dataset": int(context["n_clusters"]),
        },
        "clusters": clusters,
        "gene_clusters": by_gene,
        "gene_carriers": carriers,
        "native_gene": native,
    }


def _opening_records(sources, label: str, gene: str, k: int) -> list[dict]:
    """Sentences shipped with a packet, registered as delivered when a server is serving.

    Registering matters: the first thing an agent asks about is usually a gene it can
    already see sentences for, and a request that replied with those same sentences would
    read as a refusal to look further.
    """
    if sources is None:
        return []
    if hasattr(sources, "opening"):
        return sources.opening(label, gene, k=k)
    return sources.records_for_marker(label, gene, k=k)


def register_packet_sources(server, pool: dict, cluster: str) -> None:
    """Tell a cluster's server which sentences its opening packet already shows.

    `exclusion_sources` are resolved once, when the pool is built for every cluster at
    once, so the per-cluster server cannot have served them and would otherwise offer
    them again as though they were new. The packet is the first thing the agent reads;
    a request answered with what it just read looks like a refusal to go further.
    """
    if server is None or not hasattr(server, "note_delivered"):
        return
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        candidate = str(entry["cell_type"])
        for gene, records in (entry.get("exclusion_sources") or {}).items():
            server.note_delivered(candidate, gene, records)


def _exclusion_sources(candidate: str, markers: list[dict], sources) -> dict:
    """Curated sentences for the negative markers this cluster detects, keyed by gene."""
    out = {}
    for marker in markers:
        if marker["polarity"] != "negative":
            continue
        pct_in = marker["pct_in"]
        if pct_in is None or pct_in < NEGATIVE_SOURCE_MIN_PCT_IN:
            continue
        gene = str(marker["gene"]).upper()
        records = _opening_records(sources, candidate, gene, SOURCES_PER_GENE)
        out[gene] = [
            {
                "pmid": _identifier(record.get("pmid")),
                "pmcid": _identifier(record.get("pmcid")),
                "sentence": _clip(record.get("sentence"), SOURCE_MAX_CHARS),
            }
            for record in records
        ]
    return out


def _panel_rows(
    candidate: str,
    positive: list[str],
    negative: list[str],
    cluster: str,
    stats: dict,
    meta: dict,
    native: dict,
    specificity: dict,
) -> list[dict]:
    """One candidate's curated panel in this cluster, best-corroborated marker first.

    Positive and negative markers are one list with a polarity column rather than two
    blocks, so a reader of the packet meets the identity's strongest curated evidence
    first whichever direction it points in. Only genes this dataset measured appear; the
    rest are counted separately, because a gene with no measurement supports nothing and
    refutes nothing.

    `significant` and `specificity` carry the two quantities retrieval already decided
    this row on. Retrieval admits a marker on `sig_pass` and weighs it by `m_g`, while
    the row used to offer only `n_pub` -- which is measurably NOT a stand-in for either:
    18% of rows a reader would call raised fail `sig_pass`, and `n_pub` runs slightly
    against specificity (Spearman -0.25), so the best-published marker of an identity is
    on average the one more of its neighbours also claim. Neither field states what to
    conclude; both are the run's own numbers, arriving where the decision is made.
    """
    rows = []
    seen: set[str] = set()
    for gene, polarity in [(gene, "positive") for gene in positive] + [
        (gene, "negative") for gene in negative
    ]:
        gene = str(gene).upper()
        if gene in seen:
            continue
        measured = stats.get((cluster, gene))
        if measured is None:
            continue
        seen.add(gene)
        curated = meta.get((candidate, gene), {})
        rows.append(
            {
                "gene": native.get(gene, gene),
                "polarity": str(curated.get("polarity") or polarity),
                "n_pub": int(curated.get("n_pub", 0)),
                "tier": str(curated.get("tier") or NOT_AVAILABLE),
                "pct_in": measured["pct_in"],
                "pct_out": measured["pct_out"],
                "avg_log2FC": measured["avg_log2FC"],
                "significant": _significant(measured),
                "specificity": _num(specificity.get(gene), SPECIFICITY_DECIMALS),
            }
        )
    rows.sort(key=lambda row: (-row["n_pub"], row["gene"]))
    return rows


def candidate_names(pool: dict, cluster: str) -> list[str]:
    return [
        str(entry["cell_type"])
        for entry in pool["clusters"][str(cluster)]["candidates"]
    ]


def find_candidate(pool: dict, cluster: str, name: str) -> dict | None:
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        if str(entry["cell_type"]) == str(name):
            return entry
    return None


# --------------------------------------------------------------- opening packet ----
def _packet_marker(row: dict) -> dict:
    """A marker row as the agent sees it: the run's record minus what it cannot use."""
    return {field: row[field] for field in PACKET_MARKER_FIELDS if field in row}


def cluster_marker_table(pool: dict, cluster: str) -> list[dict]:
    """The cluster's measured evidence keyed by gene instead of by identity.

    The candidate blocks answer "is this identity's curated definition met here". They
    cannot answer "which identity does this cluster's evidence belong to", because that
    question is about a gene's claimants and each block shows one claimant. Answering it
    from the blocks means merging twenty-five panels by hand; the `candidates_with_gene`
    tool answers it one gene at a time, and over the two liver datasets the agent called
    it on 4% of clusters.

    Ordering by contrast rather than by `n_pub` is the other half. A panel is ordered by
    publication support, which is right for reading one definition and measurably wrong
    for comparing two: across the shipped library the two orders correlate at -0.08 within
    a panel, and half of all panels longer than the rows a reader meets first put their
    widest-contrast gene below that point. In one liver cluster the four Kupffer genes
    that separate it from hepatocyte sat at rows 78, 57, 35 and 76 of an 82-row panel,
    each carrying one or two publications, while the panel opened on a 107-publication
    gene measured lower inside the cluster than outside.

    Nothing is computed here that is not already on a panel row, nothing is scored and no
    row is a conclusion: the same measurements, ordered by the contrast [R4] compares on,
    with every claimant named so that how widely a gene is shared is visible where it is
    used. Significant rows only, on the pipeline's own gate, because a row that did not
    clear it is not evidence this cluster separates anything with.
    """
    rows: dict[str, dict] = {}
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        name = str(entry["cell_type"])
        for marker in entry["markers"]:
            if marker["polarity"] != "positive" or not marker["significant"]:
                continue
            gene = str(marker["gene"])
            row = rows.setdefault(
                gene,
                {
                    "gene": gene,
                    "pct_in": marker["pct_in"],
                    "pct_out": marker["pct_out"],
                    "specificity": marker["specificity"],
                    "claimed_by": [],
                },
            )
            row["claimed_by"].append({"cell_type": name, "n_pub": marker["n_pub"]})
    for row in rows.values():
        row["claimed_by"].sort(key=lambda claim: (-claim["n_pub"], claim["cell_type"]))
    ordered = sorted(
        rows.values(),
        key=lambda row: (-((row["pct_in"] or 0.0) - (row["pct_out"] or 0.0)), row["gene"]),
    )
    return ordered[:CLUSTER_MARKER_ROWS]


def cluster_packet(pool: dict, cluster: str) -> dict:
    """What the agent is given before it asks anything.

    Source sentences are the largest thing available and the smallest fraction of it that
    ever matters, so they are fetched by name -- except behind a detected exclusion, where
    the sentences are the evidence and arrive with the panel.
    """
    state = pool["clusters"][str(cluster)]
    return {
        "query": {
            **pool["context"],
            "cluster_id": str(cluster),
            "cells_in_cluster": int(state["n_cells"]),
        },
        "cluster_markers": cluster_marker_table(pool, cluster),
        "candidates": [
            {
                "cell_type": entry["cell_type"],
                "markers": [_packet_marker(row) for row in entry["markers"]],
                "exclusion_sources": entry["exclusion_sources"],
                "unmeasured_curated_genes": entry["unmeasured_curated_genes"],
                "program": entry["program"],
            }
            for entry in sorted(
                state["candidates"], key=lambda entry: str(entry["cell_type"])
            )
        ],
    }


# ------------------------------------------------------------------------ judge ----
def _judge_panel(entry: dict) -> list[dict]:
    """The identity's WHOLE curated definition as this cluster measures it.

    Ordered by publication support and NOT filtered to what the cluster raises, which is
    the whole difference between a check and a rubber stamp. A gene the literature ties
    most strongly to this identity, sitting near zero here, is the evidence AGAINST it;
    keep only the raised rows and every label arrives looking supported, because the only
    rows left are the ones that supported it. Each row carries its own measurement, so
    which of them the cluster actually raises is visible without pre-selecting them.

    Not cut to a length either, for the same reason one step further on. Cut at the
    fourteen best-published rows, 61 of the genes the annotator's own answers rested on --
    10.3% of them, over 29 of 104 clusters -- were absent from the judge's page entirely:
    C1QB and CSF1R under Kupffer cell, VSIG4 under alveolar macrophage, PKHD1 under
    cholangiocyte. A judge cannot weigh a row it was never shown, and it had no way to ask
    for one. Rows are cheap; it is the sentences that are not, and those are what it asks
    for.

    The detected exclusions come too: a curated negative marker the cluster expresses is
    the other half of the case against.
    """
    positives = [row for row in entry["markers"] if row["polarity"] == "positive"]
    detected_exclusions = [
        row
        for row in entry["markers"]
        if row["polarity"] == "negative"
        and float(row.get("pct_in") or 0.0) >= NEGATIVE_SOURCE_MIN_PCT_IN
    ]
    return [_packet_marker(row) for row in positives + detected_exclusions]


def _judge_sources(label: str, panel: list[dict], sources) -> dict:
    """The sentences that ship with the panel: the best-published, and every exclusion.

    A whole panel's sentences run to hundreds of rows and the judge reads a handful, so
    the packet opens with the ones a verdict usually turns on and the `sources` tool
    answers for any other gene on the panel. Every DETECTED EXCLUSION is pushed whatever
    its rank -- an exclusion's whole weight is what its sentence does, so a judge holding
    the row without the sentence holds the half that cannot be read.
    """
    if sources is None:
        return {}
    wanted: list[str] = []
    positives = 0
    for row in panel:
        gene = str(row["gene"]).upper()
        if row.get("polarity") == "negative":
            wanted.append(gene)
            continue
        if positives < JUDGE_SOURCES_PUSHED:
            positives += 1
            wanted.append(gene)
    found: dict[str, list[dict]] = {}
    for gene in wanted:
        if gene in found:
            continue
        records = _opening_records(sources, label, gene, JUDGE_SOURCES_PER_GENE)
        if records:
            found[gene] = [
                {
                    "pmid": _identifier(record.get("pmid")),
                    "pmcid": _identifier(record.get("pmcid")),
                    "sentence": _clip(record.get("sentence"), SOURCE_MAX_CHARS),
                }
                for record in records
            ]
    return found


def judge_packet(
    pool: dict, cluster: str, label: str, sources=None, parent: str = ""
) -> dict:
    """The evidence one judgment turns on, and deliberately nothing else.

    A judge handed the whole pool would be a second annotator: it could prefer a
    different candidate, and that disagreement says nothing about whether THIS label is
    carried by the measurements. So the packet narrows to the label under test -- the
    rows the cluster raises for it, the exclusions it detects against it, and the
    curated sentences behind those genes.

    Three things are absent by construction. The annotator's rationale, confidence and
    candidate ranking are the opinion under test, not evidence for it. The retrieval
    order would say which label the pipeline preferred. And no reference or author label
    exists anywhere in this package to leak in the first place.

    With `parent`, the question becomes whether the finer name is separable rather than
    whether it is plausible, so the parent's panel travels with it: a refinement whose
    raised markers all belong to the parent is the parent read at a narrower name.
    """
    entry = find_candidate(pool, cluster, label)
    if entry is None:
        return {}
    state = pool["clusters"][str(cluster)]
    panel = _judge_panel(entry)
    packet = {
        "query": {
            **pool["context"],
            "cluster_id": str(cluster),
            "cells_in_cluster": int(state["n_cells"]),
        },
        "label_under_test": str(label),
        "panel": panel,
        "unmeasured_curated_genes": entry["unmeasured_curated_genes"],
        "sources": _judge_sources(str(label), panel, sources),
    }
    parent_entry = find_candidate(pool, cluster, parent) if parent else None
    if parent_entry is not None:
        packet["parent_label"] = str(parent)
        packet["parent_panel"] = _judge_panel(parent_entry)
        packet["raised_and_absent_from_parent_panel"] = subtype_exclusive_markers(
            pool, cluster, str(parent), str(label)
        )
    return packet


# ------------------------------------------------------------------------ tools ----
def run_tool(pool: dict, cluster: str, tool: str, args: Any, sources=None) -> dict:
    """One tool call, answered from the artifacts. Read-only, and short by construction.

    An unknown tool or an unusable argument is answered with an error observation rather
    than an exception: the agent is mid-conversation, and telling it what went wrong lets
    it correct the call, where raising would lose the whole cluster.
    """
    args = args if isinstance(args, dict) else {}
    if tool == "sources":
        return _tool_sources(pool, cluster, args, sources)
    if tool == "gene_across_clusters":
        return _tool_gene_across_clusters(pool, args)
    if tool == "candidates_with_gene":
        return _tool_candidates_with_gene(pool, cluster, args)
    if tool == "pool_search":
        return _tool_pool_search(pool, cluster, args)
    return {"tool": str(tool), "error": f"unknown tool; available: {', '.join(TOOLS)}"}


def run_judge_tool(
    pool: dict, cluster: str, labels: tuple[str, ...], tool: str, args: Any, sources=None
) -> dict:
    """One tool call from a judge, answerable only about the labels it is judging.

    The same read-only answers the annotator gets, through a gate that keeps the judge
    from becoming a second annotator: `sources` only for a label under test, and no tool
    at all that names another identity. Without this a judge could ask for a rival's
    sentences and start choosing between them, and its verdict would no longer be about
    whether THIS label is carried.
    """
    args = args if isinstance(args, dict) else {}
    if tool not in JUDGE_TOOLS:
        return {
            "tool": str(tool),
            "error": f"unknown tool; available: {', '.join(JUDGE_TOOLS)}",
        }
    if tool == "gene_across_clusters":
        return _tool_gene_across_clusters(pool, args)
    asked = str(args.get("label") or args.get("candidate") or "")
    if asked and asked not in labels:
        return {
            "tool": "sources",
            "error": (
                f"'{asked}' is not under test here; "
                f"available: {', '.join(labels)}"
            ),
        }
    label = asked or (labels[0] if labels else "")
    return _tool_sources(pool, cluster, {**args, "candidate": label}, sources)


def _tool_sources(pool: dict, cluster: str, args: dict, sources) -> dict:
    candidate = str(args.get("candidate") or "")
    entry = find_candidate(pool, cluster, candidate)
    if entry is None:
        return {
            "tool": "sources",
            "error": f"'{candidate}' is not a candidate of this cluster",
        }
    requested = args.get("genes")
    requested = requested if isinstance(requested, list) else [requested]
    genes = [str(gene).strip() for gene in requested if str(gene or "").strip()]
    truncated = len(genes) > SOURCE_GENES_PER_QUERY
    result = {}
    left = {}
    for gene in genes[:SOURCE_GENES_PER_QUERY]:
        key = gene.upper()
        if sources is None:
            result[key] = []
            continue
        if hasattr(sources, "take"):
            answer = sources.take(candidate, key)
            records = answer["sources"]
            # Only the genes with something still behind them are reported, so the note
            # stays short and says something when it appears.
            if answer["remaining"] or answer["limit_reached"]:
                left[key] = (
                    "limit reached"
                    if answer["limit_reached"]
                    else f"{answer['remaining']} more"
                )
        else:
            records = sources.records_for_marker(candidate, key, k=SOURCES_PER_GENE)
        result[key] = [
            {
                "pmid": _identifier(record.get("pmid")),
                "pmcid": _identifier(record.get("pmcid")),
                "sentence": _clip(record.get("sentence"), SOURCE_MAX_CHARS),
            }
            for record in records
        ]
    answer = {
        "tool": "sources",
        "candidate": candidate,
        "sources": result,
        "truncated": truncated,
    }
    if left:
        answer["not_yet_shown"] = left
    return answer


def _tool_gene_across_clusters(pool: dict, args: dict) -> dict:
    gene = str(args.get("gene") or "").strip().upper()
    rows = pool["gene_clusters"].get(gene)
    if not rows:
        return {
            "tool": "gene_across_clusters",
            "gene": gene,
            "error": "not measured in this dataset",
        }
    return {
        "tool": "gene_across_clusters",
        "gene": pool["native_gene"].get(gene, gene),
        "clusters": [
            {"cluster_id": item[0], "pct_in": item[1], "pct_out": item[2]}
            for item in sorted(rows, key=lambda item: -(item[1] or 0.0))
        ],
    }


def _tool_candidates_with_gene(pool: dict, cluster: str, args: dict) -> dict:
    gene = str(args.get("gene") or "").strip().upper()
    here = set(candidate_names(pool, cluster))
    rows = [
        {"cell_type": name, "n_pub": n_pub, "tier": tier, "polarity": polarity}
        for name, n_pub, tier, polarity in pool["gene_carriers"].get(gene, [])
        if name in here
    ]
    rows.sort(key=lambda row: (-row["n_pub"], row["cell_type"]))
    return {
        "tool": "candidates_with_gene",
        "gene": pool["native_gene"].get(gene, gene),
        "candidates": rows,
    }


def _tool_pool_search(pool: dict, cluster: str, args: dict) -> dict:
    """The curated cell types claiming this gene that are NOT among the candidates.

    The one way anything downstream can see past retrieval. Retrieval admits a median of
    67 cell types per cluster and hands over 15, and the agent picked outside the top three
    in 45% of clusters and from the last four slots in 11% -- an answer rate that has not
    decayed at the boundary, so the rows beyond it are not empty, only invisible. Reading
    them is what lets a missing menu be reported instead of silently absorbed.

    What comes back is a name and a publication count, deliberately not a panel: it is
    evidence that the menu may be short, not a sixteenth candidate. Nothing here can be
    answered with -- `selected` is still checked against the supplied candidates -- and the
    agent reports it in `unlisted_identity`.
    """
    gene = str(args.get("gene") or "").strip().upper()
    here = set(candidate_names(pool, cluster))
    rows = [
        {"cell_type": name, "n_pub": n_pub, "tier": tier, "polarity": polarity}
        for name, n_pub, tier, polarity in pool["gene_carriers"].get(gene, [])
        if name not in here
    ]
    rows.sort(key=lambda row: (-row["n_pub"], row["cell_type"]))
    return {
        "tool": "pool_search",
        "gene": pool["native_gene"].get(gene, gene),
        "not_in_candidates": rows[:POOL_SEARCH_ROWS],
        "total": len(rows),
    }


# --------------------------------------------------------------------- warnings ----
def subtype_exclusive_markers(
    pool: dict, cluster: str, selected: str, subtype: str
) -> list[str]:
    """The raised markers that belong to the subtype's definition and not the parent's.

    A finer name is only a finer CALL when the cluster raises a marker the parent's own
    curated definition does not contain. A gene curated for both identities is the parent
    program read at a finer name, however wide its contrast, so the comparison is against
    the parent's WHOLE curated panel rather than the part of it this cluster raises.
    """
    parent = find_candidate(pool, cluster, selected)
    finer = find_candidate(pool, cluster, subtype)
    if parent is None or finer is None:
        return []
    parent_genes = {
        str(marker["gene"]).upper()
        for marker in parent["markers"]
        if marker["polarity"] == "positive"
    }
    return sorted(
        {
            str(marker["gene"]).upper()
            for marker in finer["markers"]
            if marker["polarity"] == "positive"
            and is_raised(marker)
            and str(marker["gene"]).upper() not in parent_genes
        }
    )


def _group_of(identity_groups, name: str) -> set[str]:
    """The candidates the answer put in one group with `name`, `name` included."""
    group = {str(name)}
    for entry in identity_groups or []:
        members = {str(member) for member in (entry.get("candidates") or [])}
        if str(name) in members:
            group |= members
    return group


def _detected_exclusions(pool: dict, cluster: str, group: set[str]) -> list[dict]:
    found: dict[str, dict] = {}
    for name in sorted(group):
        entry = find_candidate(pool, cluster, name)
        if entry is None:
            continue
        for marker in entry["markers"]:
            if marker["polarity"] != "negative":
                continue
            pct_in = marker["pct_in"]
            if pct_in is None or pct_in < NEGATIVE_SOURCE_MIN_PCT_IN:
                continue
            gene = str(marker["gene"])
            previous = found.get(gene.upper())
            if previous is None or pct_in > previous["pct_in"]:
                found[gene.upper()] = {
                    "gene": gene,
                    "curated_negative_for": name,
                    "pct_in": pct_in,
                    "pct_out": marker["pct_out"],
                    "n_pub": marker["n_pub"],
                }
    return [found[key] for key in sorted(found)]


def binding_exclusions(
    pool: dict, cluster: str, selected: str, identity_groups, co_occurring=None
) -> list[dict]:
    """Detected exclusions bound to the identity the answer selected AND to each it
    reported as merely co-occurring, one block per identity.

    `claim_warnings` performs the mirror of this check on the positive half of the panel,
    and both read the same packet the same way. The difference is when: that one runs
    after the answer is fixed and is filed for a human reader, while an exclusion on the
    identity the agent is about to assert is still answerable, and over the shipped 809
    clusters it usually went unanswered.

    Both sides are returned because one side alone cannot be read. A block whose
    `detected_exclusions` is empty says so explicitly rather than by absence: the fact
    that decides a cluster like a phagocytic macrophage called hepatocyte is not that the
    pick carries exclusions, it is that the identity it demoted carries none. Ambient
    RNA, engulfed cargo and doublets all raise a positive marker, and none of them can
    make a cell express the gene whose absence defines it -- so this asymmetry is the one
    quantity in the packet that contamination cannot manufacture.

    Whether the reason already discusses the gene is deliberately NOT consulted. Naming a
    gene is not the same as weighing it against the chosen identity -- the failure this
    check exists for quotes an exclusion while rejecting a rival and never applies it to
    its own pick, and a text match cannot tell those apart.

    Returns facts only: the gene, the candidate the resource curates it under, and the
    measurement. Nothing here says a detected exclusion refutes anything, or ranks one
    against another; that reading stays with the agent, answered against the
    `exclusion_sources` already in its packet.
    """
    selected = str(selected)
    blocks = [
        {
            "identity": selected,
            "role": "selected",
            "detected_exclusions": _detected_exclusions(
                pool, cluster, _group_of(identity_groups, selected)
            ),
        }
    ]
    seen = {selected}
    for name in co_occurring or []:
        name = str(name)
        if not name or name in seen:
            continue
        seen.add(name)
        blocks.append(
            {
                "identity": name,
                "role": "co_occurring",
                "detected_exclusions": _detected_exclusions(
                    pool, cluster, _group_of(identity_groups, name)
                ),
            }
        )
    return blocks if any(block["detected_exclusions"] for block in blocks) else []


def claim_warnings(
    pool: dict, cluster: str, claims: list[tuple[str, str]]
) -> list[tuple[str, str]]:
    """One short line per claimed identity whose curated markers are not all raised.

    The warning is computed over the identity's WHOLE measured positive panel, not over
    the gene the agent chose to quote. That is the point of it: a claim can be written on
    a single incidental marker with a wide contrast while the genes the database says
    define that cell sit flat, and a check that only looked at the quoted gene would see
    nothing wrong.

    It removes nothing and renames nothing. Publication support only orders which
    unraised genes are worth naming in the line; being unraised is the whole trigger.

    Each line travels WITH its identity rather than the identity being parsed back out
    of the line: 212 curated cell-type names carry a colon of their own (`BCR::ABL1
    positive primitive cell`), so cutting the line at its first colon filed the warning
    against an identity that does not exist.
    """
    lines = []
    for role, name in claims:
        entry = find_candidate(pool, cluster, name)
        if entry is None:
            continue
        positive = [
            marker for marker in entry["markers"] if marker["polarity"] == "positive"
        ]
        unraised = [marker for marker in positive if not is_raised(marker)]
        if not positive or not unraised:
            continue
        shown = unraised[:WARNING_GENES_SHOWN]
        listed = ",".join(
            f"{marker['gene']}({marker['pct_in']}/{marker['pct_out']})"
            for marker in shown
        )
        lines.append((
            name,
            f"{role} {name}: not_raised {len(unraised)}/{len(positive)}"
            f" | top_n_pub: {listed} | +{len(unraised) - len(shown)} more",
        ))
    return lines
