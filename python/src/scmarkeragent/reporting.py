#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Arm-agnostic production reporting for scMarkerAgent.

The module emits a concise reader-facing cluster summary, a long marker-evidence
table, and a structured per-cluster audit sidecar. It reads the retrieval artifacts and
the per-cluster annotations and never re-decides one.

Every table writes an explicit not-applicable token where a field does not apply, so a
blank cell in a delivered file always means a lost value and never an intended one.
"""

import os
import sys
import json
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import (
    CACHE_DIR as CACHE,
    RESULTS_DIR,
    FIGDATA_DIR,
    STATS_FILE,
    DEFAULTS,
    NOT_AVAILABLE,
    OUTPUT_SCHEMA,
    CLUSTER_EVIDENCE_FILE,
    na_display,
)
from . import evidence_gate
from .marker_sources import SourceDB

# ---- fixed formatting policy -------------------------------------------------
_D = DEFAULTS["output"]["decimals"]
ND_CONF = int(_D["confidence"])
ND_FRAC = int(_D["fraction"])
ND_S = int(_D["specificity"])

TABLE_A_COLS = list(OUTPUT_SCHEMA["cluster_summary"]["columns"])
TABLE_B_COLS = list(OUTPUT_SCHEMA["marker_evidence"]["columns"])
ND_LOGFC = int(_D["logfc"])
# Table A shows a shortlist; the full panels stay in marker_evidence.csv. Twelve matches
# the long-standing key_markers budget. On the 2026-08-06 freeze, 57% of assigned panels
# and 66% of co-occurring panels already fit in twelve; the rest are readable in Table B.
SUMMARY_MARKER_CAP = 12

# ================================================================= helpers ====
def _cluster_sort_key(c):
    s = str(c)
    return (0, int(s)) if s.lstrip("-").isdigit() else (1, 0, s)


def _fmt_frac(x):
    return f"{float(x):.{ND_FRAC}f}"


def _fmt_pvalue(value):
    """An adjusted p-value for a table cell, or the not-applicable token.

    Kept in the float's own notation rather than rounded to a fixed number of places,
    which would read every strongly raised marker as exactly zero. This is the scale
    `markers_all_by_cluster.csv` already reports padj on, so the two tables can be
    compared gene by gene.
    """
    if value is None:
        return NOT_AVAILABLE
    try:
        number = float(value)
    except (TypeError, ValueError):
        return NOT_AVAILABLE
    return NOT_AVAILABLE if not np.isfinite(number) else repr(number)


def _fmt(value, decimals):
    """A number formatted for a table cell, or the not-applicable token."""
    if value is None:
        return NOT_AVAILABLE
    try:
        number = float(value)
    except (TypeError, ValueError):
        return NOT_AVAILABLE
    return NOT_AVAILABLE if not np.isfinite(number) else f"{number:.{decimals}f}"


def _fmt_alternative_candidates(names):
    """The other identities the cluster also carries, ';'-joined.

    What a reader needs from a summary row is which OTHER cells are in this cluster. The
    full candidate set, each with its claim role and retrieval rank, is in
    `cluster_evidence.jsonl`, where it can be queried instead of widening this column past
    the point a reader can use it.
    """
    parts = [str(name) for name in names or [] if str(name).strip()]
    return "; ".join(parts) if parts else NOT_AVAILABLE


def _fmt_claim_warnings(lines):
    """The audit's short note on claims whose curated markers this cluster does not raise.

    It is a flag for a reader, not a correction: the annotation is exactly what the agent
    decided, and this column says which of its claims rest on an identity whose
    best-corroborated markers sit unraised here. An empty column means every claimed
    identity had its whole measured positive panel up-regulated.
    """
    parts = [str(line) for line in lines or [] if str(line).strip()]
    return " || ".join(parts) if parts else NOT_AVAILABLE


def _json_safe(value):
    """Convert pandas/numpy values to strict, portable JSON values."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, np.ndarray)):
        return [_json_safe(item) for item in value]
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


def _load_run_artifacts(tag):
    """Load the frozen run artifacts every builder reads. Nothing is re-decided here."""
    with open(os.path.join(CACHE, f"{tag}_de_meta.pkl"), "rb") as handle:
        dm = pickle.load(handle)
    with open(os.path.join(CACHE, f"{tag}_candidate_scoring.pkl"), "rb") as handle:
        scoring = pickle.load(handle)
    with open(os.path.join(CACHE, f"{tag}_annotations.pkl"), "rb") as handle:
        annotations = pickle.load(handle)
    norm = sp.load_npz(os.path.join(CACHE, f"{tag}_norm.npz")).tocsr()
    return dict(
        tag=tag,
        dm=dm,
        scoring=scoring,
        annotations=annotations,
        res=annotations["results"],
        norm=norm,
    )


DOTPLOT_COLS = [
    "gene",
    "marker_slot",
    "gene_order",
    "gene_group",
    "gene_group_order",
    "cluster",
    "cluster_order",
    "cluster_celltype",
    "pct_exp",
    "avg_exp_scaled",
]


# ================================================================ builders ====
def build_clustermap(fr):
    """Per-cluster annotation intermediate shared with the figure bundle.

    Both arms consume these labels rather than re-deriving them, so no plot can disagree
    with the delivered table.
    """
    res = fr["res"]
    cl_id = fr["scoring"].get("candidate_cl_id", {})
    rows = []
    for cluster in sorted(res.keys(), key=_cluster_sort_key):
        record = res[cluster]
        label = str(record.get("annotation") or NOT_AVAILABLE)
        subtype = str(record.get("subtype") or "")
        # Finest evidence-backed identity: subtype when present, otherwise the type.
        # Plots (UMAP / identity-marker dotplot) and benchmark headline scoring all
        # read this through cluster_celltype / primary_annotation; gene panels stay
        # unchanged (marker_evidence still carries both selected levels).
        primary = subtype or label
        rows.append(
            dict(
                cluster=str(cluster),
                n_cells=int(fr["scoring"]["clusters"][str(cluster)]["n_cells"]),
                cell_type_annotation=label,
                cell_subtype_annotation=na_display(record.get("subtype")),
                # The one label to read, plot and score against: the finest level the
                # evidence established. The two columns above keep the levels apart for a
                # reader who needs to see which was which; this one answers "what is this
                # cluster" without the reader having to combine them.
                primary_annotation=primary,
                cell_state=na_display(record.get("state")),
                cell_ontology=na_display(cl_id.get(label)),
                annotation_confidence=na_display(record.get("confidence")),
                annotation_rationale=na_display(record.get("rationale")),
                resolution_status=str(record.get("resolution_status") or ""),
                resolution_detail=str(record.get("resolution_detail") or ""),
                annotation_source=str(record.get("annotation_source") or ""),
                llm_status=str(record.get("llm_status") or ""),
                annotation_qc=str(record.get("annotation_qc") or ""),
                cluster_celltype=f"{cluster}: {primary}",
                alternative_candidates=_fmt_alternative_candidates(
                    record.get("co_occurring_identities")
                ),
                claim_warnings=_fmt_claim_warnings(record.get("claim_warnings")),
                unlisted_identity=na_display(record.get("unlisted_identity")),
                identity_groups_json=json.dumps(
                    record.get("identity_groups") or [], ensure_ascii=False
                ),
            )
        )
    clustermap = pd.DataFrame(rows)
    clustermap.insert(1, "cluster_order", range(len(clustermap)))
    return clustermap


def _by_evidence(frame):
    """Rank a cluster's marker rows so a cap keeps the strongest, best-corroborated ones.

    Table B is written gene-ascending because it is a full record and a reader looks genes
    up in it; nothing is lost there. Table A shows twelve, and until 2026-08-06 it took the
    first twelve of that alphabetical order. On `human_bladder` cluster 1, called
    `endothelial cell`, that dropped PECAM1 (2,640 publications, detected in 91.9%, log2FC
    3.20) and VWF (1,169, 91.3%, 4.16) -- the two textbook markers of the identity -- and
    kept F11R, which has one publication and a log2FC of 0.18, because P and V sort after
    K. The summary is the column a reader meets first, so it is ordered the way every other
    panel in this package is: publication support first, contrast to break ties.
    """
    ranked = frame
    if "marker_polarity" in ranked.columns:
        # `key_markers` is the support for the call. Detected exclusions are in the
        # evidence table and in the rationale that weighs them; a reader meeting one here,
        # in a list of reasons to believe the label, would read it as the opposite of what
        # it is.
        ranked = ranked[ranked["marker_polarity"].astype(str) != "negative"]
    ranked = ranked.copy()
    ranked["_support"] = pd.to_numeric(
        ranked["publication_support_count"], errors="coerce"
    ).fillna(0.0)
    ranked["_contrast"] = pd.to_numeric(
        ranked["cluster_detection_fraction"], errors="coerce"
    ).fillna(0.0) - pd.to_numeric(
        ranked["out_of_cluster_detection_fraction"], errors="coerce"
    ).fillna(0.0)
    return (
        ranked.sort_values(
            ["_support", "_contrast", "gene"],
            ascending=[False, False, True],
            kind="stable",
        )
        .drop(columns=["_support", "_contrast"])
    )


def _summary_marker_rows(group, tag_identity: bool):
    """Shortlist rows for one summary column, fair-shared across identities in `group`.

    A single global corroboration sort lets a high-n_pub parent (or one co-occurring
    identity) fill the whole cap and hide the markers that actually define a subtype or a
    second population. Each identity in the group therefore gets a share of
    SUMMARY_MARKER_CAP first; leftovers fill in the same corroboration order.
    """
    if group is None or len(group) == 0:
        return "", ""
    ranked_by_identity = {}
    for identity, own in group.groupby("candidate_annotation", sort=False):
        ranked_by_identity[str(identity)] = list(
            _by_evidence(own).itertuples(index=False)
        )
    identities = list(ranked_by_identity)
    share = max(1, SUMMARY_MARKER_CAP // len(identities))
    chosen = []
    seen = set()

    def push(row) -> bool:
        # One slot per gene, not per (identity, gene): a gene curated for both a parent
        # and its subtype -- or for two co-occurring identities -- used to be printed
        # twice, which reads as two pieces of evidence and spends two of the twelve slots
        # on one measurement. The first occurrence wins, so the identity holding the
        # earlier fair-share slot keeps it.
        gene = str(row.gene)
        if gene in seen:
            return False
        seen.add(gene)
        chosen.append(row)
        return len(chosen) >= SUMMARY_MARKER_CAP

    for identity in identities:
        for row in ranked_by_identity[identity][:share]:
            if push(row):
                return _format_summary_markers(chosen, tag_identity)
    for identity in identities:
        for row in ranked_by_identity[identity][share:]:
            if push(row):
                return _format_summary_markers(chosen, tag_identity)
    return _format_summary_markers(chosen, tag_identity)


def _format_summary_markers(rows, tag_identity: bool):
    if not rows:
        return "", ""
    if tag_identity:
        # Pipe fields, not gene(identity;fraction): curated free-text names routinely
        # contain parentheses (e.g. "Mac-2 (Atf3)"), which made the older spelling
        # unparseable. Gene symbols and the 581 names in the freeze never contain `|`.
        markers = "; ".join(
            f"{row.gene}|{row.candidate_annotation}|{row.cluster_detection_fraction}"
            for row in rows
        )
    else:
        markers = "; ".join(
            f"{row.gene}({row.cluster_detection_fraction})" for row in rows
        )
    # Publications behind exactly the markers this column shows; order follows the
    # markers and duplicates are dropped, because one paper often backs several of them.
    pmcids = "; ".join(
        dict.fromkeys(
            str(row.pmcid)
            for row in rows
            if str(getattr(row, "pmcid", "")) not in {"", NOT_AVAILABLE}
        )
    )
    return markers, pmcids


def _cooccurring_names(alternative_candidates) -> set[str]:
    """Names listed in alternative_candidates — the only identities cooc columns may show.

    marker_evidence also keeps rejected_subtype rows (and similar non-selected claims).
    Those are audit trail, not co-occurring populations; putting them in
    cooccurring_markers would invent a mix the delivery never declared.
    """
    text = str(alternative_candidates or "").strip()
    if not text or text == NOT_AVAILABLE:
        return set()
    return {part.strip() for part in text.split(";") if part.strip() and part.strip() != NOT_AVAILABLE}


def build_table_a(tag, clustermap, marker_evidence):
    """Table A -- concise reader-facing summary, one row per dataset and cluster."""
    by_cluster = {
        str(cluster): group
        for cluster, group in marker_evidence.groupby("cluster_id", sort=False)
    }
    out = []
    for _, row in clustermap.iterrows():
        cluster = str(row["cluster"])
        group = by_cluster.get(cluster)
        if group is None or group.empty:
            key_markers = pmcids = co_markers = co_pmcids = ""
        else:
            selected = group[group["is_selected_annotation"].astype(str) == "True"]
            source = selected if len(selected) else group
            key_markers, pmcids = _summary_marker_rows(source, tag_identity=False)
            cooc_names = _cooccurring_names(row.get("alternative_candidates"))
            if cooc_names:
                cooccurring = group[
                    group["candidate_annotation"].astype(str).isin(cooc_names)
                ]
                co_markers, co_pmcids = _summary_marker_rows(
                    cooccurring, tag_identity=True
                )
            else:
                co_markers = co_pmcids = ""
        out.append(
            {
                "dataset": tag,
                "cluster_id": cluster,
                "n_cells": int(row["n_cells"]),
                "primary_annotation": row["primary_annotation"],
                "cell_type_annotation": row["cell_type_annotation"],
                "cell_subtype_annotation": row["cell_subtype_annotation"],
                "cell_state": row["cell_state"],
                "cell_ontology": row["cell_ontology"],
                "annotation_confidence": row["annotation_confidence"],
                "annotation_rationale": row["annotation_rationale"],
                "resolution_status": row["resolution_status"],
                "annotation_source": row["annotation_source"],
                "llm_status": row["llm_status"],
                "annotation_qc": row["annotation_qc"],
                "key_markers": na_display(key_markers),
                "pmcid": na_display(pmcids),
                "cooccurring_markers": na_display(co_markers),
                "cooccurring_pmcid": na_display(co_pmcids),
                "alternative_candidates": row["alternative_candidates"],
                "claim_warnings": row["claim_warnings"],
                "unlisted_identity": row["unlisted_identity"],
            }
        )
    return pd.DataFrame(out, columns=TABLE_A_COLS)


def build_cluster_evidence(tag, fr, marker_evidence):
    """Build the structured per-cluster audit sidecar kept out of the summary CSV.

    It carries what the summary cannot: every retrieved candidate with its three relative
    scores and its complete measured panel, and the grouping the agent used to tell
    "several names for one cell" from "several different cells".
    """
    marker_records = {
        str(cluster): group.to_dict("records")
        for cluster, group in marker_evidence.groupby("cluster_id", sort=False)
    }
    scored = fr["scoring"]["scored"]
    cl_id = fr["scoring"].get("candidate_cl_id", {})
    records = []
    for cluster in sorted(fr["res"], key=_cluster_sort_key):
        cluster_id = str(cluster)
        result = fr["res"][cluster]
        label = str(result.get("annotation") or NOT_AVAILABLE)
        pool = scored[scored["cluster"] == cluster_id]
        state = fr["scoring"]["clusters"].get(cluster_id) or {}
        record = {
            "schema_version": "scmarkeragent-cluster-evidence-v2",
            "dataset": tag,
            "cluster_id": cluster_id,
            "annotation": {
                "primary_annotation": str(result.get("subtype") or "") or label,
                "cell_type_annotation": label,
                "cell_subtype_annotation": na_display(result.get("subtype")),
                "cell_state": na_display(result.get("state")),
                "co_occurring_identities": list(
                    result.get("co_occurring_identities") or []
                ),
                "cell_ontology": na_display(cl_id.get(label)),
                "annotation_confidence": na_display(result.get("confidence")),
                "resolution_status": str(result.get("resolution_status") or ""),
                "annotation_source": str(result.get("annotation_source") or ""),
                "llm_status": str(result.get("llm_status") or ""),
                "support_markers": list(result.get("support_markers") or []),
                "claim_warnings": list(result.get("claim_warnings") or []),
            },
            "evidence": {
                "annotation_rationale": na_display(result.get("rationale")),
                "key_markers": marker_records.get(cluster_id, []),
                "candidates": result.get("candidate_entries") or [],
                "identity_groups": result.get("identity_groups") or [],
            },
            "audit": {
                "annotator_turns": int(result.get("turns", 0) or 0),
                "tool_calls": list(result.get("tool_calls") or []),
                # Exclusions returned to the agent before its answer was accepted; empty
                # when the first answer already addressed every one that binds its pick.
                "exclusion_audit": list(result.get("exclusion_audit") or []),
                # Both quality-control verdicts in full, with the genes each turned on.
                # The summary column says only how the label left the check; a reader
                # asking why it was demoted or what the judge disputed reads it here.
                "judge": dict(result.get("quality_control") or {}),
                # The reader-facing column has three values; the condition that produced
                # one of them is kept here rather than splintering that column.
                "resolution_detail": str(result.get("resolution_detail") or ""),
                "candidate_pool_size": int(len(pool)),
                # How many candidates the significant-hit gate saw and which threshold it
                # ended up applying. Without both, a small pool cannot be told apart from
                # a thin curated menu, and a relaxed gate would be invisible.
                "candidate_pool_size_before_hits_gate": int(
                    state.get("pool_size_before_hits_gate", len(pool))
                ),
                "hits_threshold_applied": int(state.get("hits_threshold_applied", 1)),
                "candidates_shown": len(result.get("candidates") or []),
                "retrieval_order": [
                    {
                        "candidate": str(row.candidate),
                        "retrieval_rank": int(row.retrieval_rank),
                        "retrieval_score": round(float(row.retrieval_score), ND_CONF),
                        "marker_level": round(float(row.marker_level), ND_CONF),
                        "cluster_level": round(float(row.cluster_level), ND_CONF),
                        "single_cell_level": round(
                            float(row.single_cell_level), ND_CONF
                        ),
                        "significant_markers": int(row.hits),
                        "measured_panel_size": int(row.panel_size),
                    }
                    for row in pool.head(
                        int(fr["scoring"]["top_candidates"])
                    ).itertuples(index=False)
                ],
            },
        }
        records.append(_json_safe(record))
    return records


def write_cluster_evidence(records, outdir):
    path = os.path.join(outdir, CLUSTER_EVIDENCE_FILE)
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            handle.write(
                json.dumps(
                    record,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
                + "\n"
            )
    return os.path.abspath(path)


def _marker_role(entry: dict, result: dict, gene: str) -> str:
    """What part this gene played in the claim the annotating agent made.

    `decisive` is the gene the agent quoted for that identity, `support` a gene it named
    behind the selected label, and `panel` a curated marker of the identity that the
    cluster raises without either having been cited.
    """
    upper = gene.upper()
    evidence = entry.get("claim_evidence") or {}
    if str(evidence.get("decisive_gene", "")).upper() == upper:
        return "decisive"
    if upper in {
        str(value).upper() for value in result.get("support_markers", []) or []
    }:
        return "support"
    return "panel"


def build_marker_evidence(tag, fr):
    """Table B -- the measured, source-backed markers behind every identity claimed.

    A row exists for each raised curated marker of every identity the agent asserted the
    cluster carries: the assigned one, its subtype, and each co-occurring identity. The
    candidates that were considered and not claimed keep their complete panels in
    `cluster_evidence.jsonl`, which is where the reason a cluster was NOT called
    something now lives -- the whole pool is on that record rather than only the
    survivors of a filter.

    One kind of unclaimed candidate does put rows here: a sibling in the assigned
    identity's own group, and only its DETECTED EXCLUSIONS. The audit turn binds those to
    the assigned identity -- a cluster called `basal cell of prostate epithelium` is asked
    to answer the KLK3 curated as absent from `basal cell`, because the agent placed the
    two names on one identity -- and the delivered reason then argues from them. Of the
    306 detected exclusions a delivered reason named, 216 were curated under a group
    sibling rather than under the assigned name, so leaving them out kept two thirds of
    the exclusion argument unverifiable. They carry the sibling in
    `candidate_annotation`, so a reader can see which curated entry the exclusion came
    from, and `marker_polarity` tells them apart from anything the cluster is claimed to
    be.
    """
    context = fr["scoring"]["context"]
    sources = SourceDB().context(
        context["species"], context["tissue"], context["disease"]
    )
    specificity = fr["scoring"]["marker_specificity"]
    # padj is decided during differential expression and never travels on the marker rows
    # the agent reads, so it is read back here from the same DE table the significance
    # gate used. Both arms hand this builder the same frame, so both fill the column.
    de_frame = fr["dm"]["de"]
    padj_of = {
        (str(group), str(feature).upper()): padj
        for group, feature, padj in zip(
            de_frame["group"], de_frame["feature"], de_frame["padj"]
        )
    }
    rows = []
    for cluster in sorted(fr["res"], key=_cluster_sort_key):
        result = fr["res"][cluster]
        # The assigned identity is now reported at two levels, and the markers behind the
        # finer one are the markers behind the call. Both count as selected rows so the
        # table and the summary's key markers cannot disagree about what was assigned.
        assigned = {
            str(result.get("annotation") or ""),
            str(result.get("subtype") or ""),
        } - {""}
        # The names the agent grouped with the assigned identity. Their exclusions are the
        # ones the audit turn made it answer, so they belong with its evidence.
        grouped = set()
        for group in result.get("identity_groups") or []:
            names = {str(name) for name in (group.get("candidates") or [])}
            if names & assigned:
                grouped |= names
        for entry in result.get("candidate_entries") or []:
            candidate = str(entry["cell_type"])
            is_selected = candidate in assigned
            claimed = is_selected or bool(str(entry.get("claim_role") or ""))
            if not claimed and candidate not in grouped:
                continue
            for marker in entry.get("decisive_marker_measurements") or []:
                # An unclaimed sibling is here only for what it excludes.
                if not claimed and str(marker.get("polarity") or "") != "negative":
                    continue
                gene = str(marker["gene"])
                records = sources.records_for_marker(candidate, gene, k=1)
                source = records[0] if records else {}
                rows.append(
                    {
                        "dataset": tag,
                        "cluster_id": str(cluster),
                        "candidate_annotation": candidate,
                        "is_selected_annotation": str(bool(is_selected)),
                        "gene": gene,
                        "marker_polarity": str(marker.get("polarity") or "positive"),
                        "marker_role": _marker_role(entry, result, gene),
                        "marker_provenance": "db_cited",
                        "cluster_detection_fraction": _fmt(
                            marker.get("detection_fraction_in"), ND_FRAC
                        ),
                        "out_of_cluster_detection_fraction": _fmt(
                            marker.get("detection_fraction_out"), ND_FRAC
                        ),
                        "average_log2_fold_change": _fmt(
                            marker.get("avg_log2FC"), ND_LOGFC
                        ),
                        "adjusted_p_value": _fmt_pvalue(
                            padj_of.get((str(cluster), gene.upper()))
                        ),
                        "cross_cluster_percentile": _fmt(
                            marker.get("cross_cluster_percentile"), ND_FRAC
                        ),
                        "marker_specificity_weight": _fmt(
                            specificity.get(gene.upper()), ND_S
                        ),
                        "publication_support_count": marker.get(
                            "publication_support", NOT_AVAILABLE
                        ),
                        "evidence_tier": na_display(marker.get("evidence_tier")),
                        "pmid": na_display(source.get("pmid")),
                        "pmcid": na_display(source.get("pmcid")),
                        "source_sentence": na_display(source.get("sentence")),
                    }
                )
    table = pd.DataFrame(rows, columns=TABLE_B_COLS)
    if len(table):
        table["_k"] = table["cluster_id"].map(_cluster_sort_key)
        table = (
            table.sort_values(
                ["_k", "is_selected_annotation", "candidate_annotation", "gene"],
                ascending=[True, False, True, True],
                kind="stable",
            )
            .drop(columns="_k")
            .drop_duplicates(["cluster_id", "candidate_annotation", "gene"])
        )
    return table.reset_index(drop=True)


def _supporting_panel(marker_evidence, owner_of):
    """The markers that support each owner's assigned identity, one slot per gene.

    The panel is what `marker_evidence` already records for the assigned identity: every
    curated marker of that identity this cluster raises. It is taken whole rather than
    narrowed to the most specific few, because specificity alone does not support a call
    -- the identity's best-corroborated markers are what do, and a reader has to see the
    same rows the table lists.

    Genes are deduplicated across owners: the first owner in display order keeps the slot
    and later owners read that same column, so a marker shared by several identities is
    one column rather than one per identity. `owner_of` maps a row of the table to the
    block it belongs to; publication support orders the genes inside a block, matching the
    order the agent read its panel in.
    """
    if marker_evidence.empty:
        return [], {}
    panel = marker_evidence[
        marker_evidence["is_selected_annotation"].astype(str).str.lower() == "true"
    ].copy()
    if "marker_polarity" in panel.columns:
        # A dot on this plot reads as evidence for the identity. A curated exclusion the
        # cluster detects is evidence against it, and plotting the two in one grid would
        # say the opposite of what the row means.
        panel = panel[panel["marker_polarity"].astype(str) != "negative"]
    if panel.empty:
        return [], {}
    panel["owner"] = panel["cluster_id"].map(owner_of)
    panel = panel[panel["owner"].notna()]
    panel["gene_upper"] = panel["gene"].astype(str).str.upper()
    panel["support"] = pd.to_numeric(
        panel["publication_support_count"], errors="coerce"
    ).fillna(0)
    slots, publications, seen = [], {}, set()
    for owner in dict.fromkeys(owner_of.values()):
        group = panel[panel["owner"] == owner]
        if group.empty:
            continue
        ranked = (
            group.groupby("gene_upper", as_index=False)
            .agg(support=("support", "max"))
            .sort_values(["support", "gene_upper"], ascending=[False, True], kind="stable")
        )
        for gene, support in zip(ranked["gene_upper"], ranked["support"]):
            publications[(owner, str(gene))] = int(support)
            if gene in seen:
                continue
            seen.add(gene)
            slots.append((owner, str(gene)))
    return slots, publications


def build_identity_marker_dotplot(fr, clustermap, marker_evidence):
    """Build a cluster-aligned dotplot from source-backed identity markers.

    Every summary row remains a dotplot row, and every marker block belongs to
    one cluster.
    """
    dm = fr["dm"]
    menu = list(dm["menu_genes"])
    gene_index = {g.upper(): i for i, g in enumerate(menu)}
    meta = dm["meta"].copy()
    cell_cluster = meta["cluster"].astype(str).values
    order = clustermap.sort_values("cluster_order")
    cluster_order = order["cluster"].astype(str).tolist()
    cluster_labels = dict(zip(cluster_order, order["cluster_celltype"].astype(str)))
    row_indices = {
        cluster: np.where(cell_cluster == cluster)[0] for cluster in cluster_order
    }

    candidates = marker_evidence.copy()
    if candidates.empty:
        return pd.DataFrame(columns=DOTPLOT_COLS)
    candidates["cluster_id"] = candidates["cluster_id"].astype(str)
    candidates = candidates[
        candidates["gene"].astype(str).str.upper().isin(gene_index)
    ].copy()
    marker_slots, _ = _supporting_panel(
        candidates, {cluster: cluster for cluster in cluster_order}
    )

    rows = []
    cluster_rank = {cluster: index for index, cluster in enumerate(cluster_order)}
    for gene_order, (owner_cluster, gene) in enumerate(marker_slots):
        expression = fr["norm"][:, gene_index[gene]].toarray().ravel()
        averages = np.array(
            [
                (
                    float(np.mean(np.expm1(expression[row_indices[cluster]])))
                    if len(row_indices[cluster])
                    else 0.0
                )
                for cluster in cluster_order
            ]
        )
        percentages = np.array(
            [
                (
                    100.0 * float(np.mean(expression[row_indices[cluster]] > 0))
                    if len(row_indices[cluster])
                    else 0.0
                )
                for cluster in cluster_order
            ]
        )
        standard_deviation = averages.std()
        scaled = (
            (averages - averages.mean()) / standard_deviation
            if standard_deviation > 0
            else np.zeros_like(averages)
        )
        scaled = np.clip(scaled, -2.5, 2.5)
        marker_slot = f"{owner_cluster}\x1f{gene}"
        owner_label = cluster_labels[owner_cluster]
        for row_order, cluster in enumerate(cluster_order):
            rows.append(
                {
                    "gene": gene,
                    "marker_slot": marker_slot,
                    "gene_order": gene_order,
                    "gene_group": owner_label,
                    "gene_group_order": cluster_rank[owner_cluster],
                    "cluster": cluster,
                    "cluster_order": row_order,
                    "cluster_celltype": cluster_labels[cluster],
                    "pct_exp": round(float(percentages[row_order]), 3),
                    "avg_exp_scaled": round(float(scaled[row_order]), 4),
                }
            )
    return pd.DataFrame(rows, columns=DOTPLOT_COLS)


CELLTYPE_DOTPLOT_COLS = [
    "gene",
    "marker_slot",
    "gene_order",
    "gene_group",
    "gene_group_order",
    "cell_type",
    "cell_type_order",
    "n_clusters",
    "n_cells",
    "pct_exp",
    "avg_exp_scaled",
    "n_pub",
]


def build_celltype_marker_dotplot(fr, clustermap, marker_evidence):
    """The identity dotplot with cell types, not clusters, as rows.

    Clusters that resolved to the SAME cell type are one row here, with detection and
    mean expression recomputed over their pooled cells rather than averaged across
    clusters, so a cell type split across a large and a small cluster is not
    misrepresented by the small one. Genes are the union of the identity markers those
    clusters contributed, deduplicated to a single slot.

    `n_pub` is the publication support `marker_evidence.csv` already reports for that
    (cell type, gene) pair, so the support-weighted twin and the table cannot disagree, and
    a free-word or finer label -- which no curated candidate name matches -- still carries
    the count of the panel it was actually called on. A gene curated for two cell types
    carries its own count under each; a pair with no row there is 0, not blank.
    """
    dm = fr["dm"]
    menu = list(dm["menu_genes"])
    gene_index = {g.upper(): i for i, g in enumerate(menu)}
    cell_cluster = dm["meta"]["cluster"].astype(str).values

    order = clustermap.sort_values("cluster_order")
    # The figure carries the label a reader is meant to take away, which is the finest level
    # the evidence established rather than the coarser of the two identity columns.
    label_of = {
        str(row.cluster): str(row.primary_annotation) or str(row.cluster)
        for row in order.itertuples(index=False)
    }
    # First appearance in cluster order fixes the row order, so the cell-type plot reads
    # in the same sequence as the cluster plot.
    cell_types = list(dict.fromkeys(label_of[str(c)] for c in order["cluster"]))
    clusters_of = {
        cell_type: [c for c in order["cluster"].astype(str) if label_of[c] == cell_type]
        for cell_type in cell_types
    }
    row_indices = {
        cell_type: np.where(np.isin(cell_cluster, clusters_of[cell_type]))[0]
        for cell_type in cell_types
    }

    candidates = marker_evidence.copy()
    if candidates.empty:
        return pd.DataFrame(columns=CELLTYPE_DOTPLOT_COLS)
    candidates["cluster_id"] = candidates["cluster_id"].astype(str)
    candidates = candidates[
        candidates["gene"].astype(str).str.upper().isin(gene_index)
    ].copy()
    marker_slots, publications = _supporting_panel(candidates, label_of)

    rows = []
    type_rank = {cell_type: index for index, cell_type in enumerate(cell_types)}
    for gene_order, (owner, gene) in enumerate(marker_slots):
        expression = fr["norm"][:, gene_index[gene]].toarray().ravel()
        averages = np.array(
            [
                (
                    float(np.mean(np.expm1(expression[row_indices[cell_type]])))
                    if len(row_indices[cell_type])
                    else 0.0
                )
                for cell_type in cell_types
            ]
        )
        percentages = np.array(
            [
                (
                    100.0 * float(np.mean(expression[row_indices[cell_type]] > 0))
                    if len(row_indices[cell_type])
                    else 0.0
                )
                for cell_type in cell_types
            ]
        )
        standard_deviation = averages.std()
        scaled = (
            (averages - averages.mean()) / standard_deviation
            if standard_deviation > 0
            else np.zeros_like(averages)
        )
        scaled = np.clip(scaled, -2.5, 2.5)
        for row_order, cell_type in enumerate(cell_types):
            rows.append(
                {
                    "gene": gene,
                    "marker_slot": f"{owner}\x1f{gene}",
                    "gene_order": gene_order,
                    "gene_group": owner,
                    "gene_group_order": type_rank[owner],
                    "cell_type": cell_type,
                    "cell_type_order": row_order,
                    "n_clusters": len(clusters_of[cell_type]),
                    "n_cells": int(len(row_indices[cell_type])),
                    "pct_exp": round(float(percentages[row_order]), 3),
                    "avg_exp_scaled": round(float(scaled[row_order]), 4),
                    "n_pub": int(publications.get((cell_type, gene), 0)),
                }
            )
    return pd.DataFrame(rows, columns=CELLTYPE_DOTPLOT_COLS)


MARKER_LIST_COLS = [
    "cluster",
    "cluster_celltype",
    "gene",
    "avg_log2FC",
    "pct_in",
    "pct_out",
    "auc",
    "mean_log_expression",
    "log_expression_difference",
    "pval",
    "padj",
    "in_marker_menu",
    "significant",
]


def build_marker_lists(tag, clustermap):
    """Per-cluster marker tables: every tested gene, and the significant subset.

    The source is the genome-wide FindAllMarkers-equivalent table written by
    preprocessing (both arms), NOT the menu-restricted DE that drives annotation. The two
    differ on purpose: the annotation gate is applied inside the marker menu so the
    per-cluster BH family matches gene for gene across arms, while this table asks the
    broader question of what separates a cluster in the whole transcriptome and so carries
    its own genome-wide BH family. `in_marker_menu` marks the rows that also belong to the
    menu universe, so a reader can tell the two apart at a glance.

    `significant` applies the shared marker gate to THIS table's genome-wide padj.
    """
    path = os.path.join(CACHE, f"{tag}_markers_all.pkl")
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path}: the genome-wide marker table is missing; rebuild prep"
        )
    with open(path, "rb") as handle:
        markers = pickle.load(handle)
    if not len(markers):
        return pd.DataFrame(columns=MARKER_LIST_COLS), pd.DataFrame(
            columns=MARKER_LIST_COLS
        )
    with open(os.path.join(CACHE, f"{tag}_de_meta.pkl"), "rb") as handle:
        menu_genes = {str(g).upper() for g in pickle.load(handle)["menu_genes"]}
    labels = dict(
        zip(
            clustermap["cluster"].astype(str),
            clustermap["cluster_celltype"].astype(str),
        )
    )
    out = pd.DataFrame(
        {
            "cluster": markers["group"].astype(str),
            "gene": markers["feature"].astype(str),
            "avg_log2FC": markers["avg_log2FC"].astype(float).round(ND_LOGFC),
            "pct_in": (markers["pct_in"].astype(float) / 100.0).round(ND_FRAC),
            "pct_out": (markers["pct_out"].astype(float) / 100.0).round(ND_FRAC),
            "auc": markers["auc"].astype(float).round(ND_FRAC),
            "mean_log_expression": markers["avgExpr"].astype(float).round(ND_LOGFC),
            "log_expression_difference": markers["logFC"].astype(float).round(ND_LOGFC),
            "pval": markers["pval"].astype(float),
            "padj": markers["padj"].astype(float),
        }
    )
    out["cluster_celltype"] = out["cluster"].map(labels).fillna(out["cluster"])
    out["in_marker_menu"] = out["gene"].str.upper().isin(menu_genes)
    out["significant"] = [
        evidence_gate.sig_pass(lfc, pin, padj)
        for lfc, pin, padj in zip(
            markers["avg_log2FC"].astype(float),
            markers["pct_in"].astype(float) / 100.0,
            markers["padj"].astype(float),
        )
    ]
    out["_k"] = out["cluster"].map(_cluster_sort_key)
    out = (
        out.sort_values(
            ["_k", "significant", "avg_log2FC"],
            ascending=[True, False, False],
            kind="stable",
        )
        .drop(columns="_k")
        .reset_index(drop=True)[MARKER_LIST_COLS]
    )
    return out, out[out["significant"]].reset_index(drop=True)


def write_stats_note(tag, outdir, clustermap, cells, table_a, table_b):
    """Write a concise, user-facing run and provenance summary."""
    statuses = clustermap["resolution_status"].value_counts().to_dict()
    lines = [
        f"scMarkerAgent result summary: {tag}",
        "",
        "Pipeline stages",
        "1. Preprocessing: QC, normalization, Leiden partitioning, genome-wide DE",
        "2. Candidate retrieval: marker-level, cluster-level and single-cell-level",
        "   relative evidence, combined into one retrieval order",
        "3. One annotating agent per cluster, opening with every candidate's complete",
        "   measured marker panel and free to query sources and cross-cluster expression",
        "4. Post-hoc audit: quoted measurements checked against the DE table, and claims",
        "   whose curated markers are unraised flagged without being altered",
        "",
        f"Clusters: {len(clustermap)}",
        f"Cells scored: {len(cells)}",
        f"Cluster summary rows: {len(table_a)}",
        f"Marker evidence rows: {len(table_b)}",
        "Resolution: "
        + ", ".join(f"{name} {count}" for name, count in sorted(statuses.items())),
        "",
        f"Structured cluster audit: {CLUSTER_EVIDENCE_FILE}",
        "Every retrieved candidate, its three relative scores, its complete measured",
        "marker panel and the agent's identity grouping are kept in that JSON Lines",
        "sidecar rather than embedded in cluster_summary.csv.",
        "",
        "Marker evidence is measured in this dataset and carries its curated source. The",
        "retrieval order is a search order over evidence, not a probability and not a",
        "confidence. No author reference labels are used anywhere in this package.",
        f"A field that does not apply to a row is written as {NOT_AVAILABLE}.",
    ]
    path = os.path.join(outdir, STATS_FILE)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return os.path.abspath(path)


def build_cells(fr):
    """Reference-free per-cell figure data: stable ID, cluster, and UMAP."""
    dm = fr["dm"]
    meta = dm["meta"].copy()
    umap = np.asarray(dm.get("umap")) if dm.get("umap") is not None else None
    cells = pd.DataFrame(
        {
            "cell": meta["cell"].astype(str).values,
            "cluster": meta["cluster"].astype(str).values,
        }
    )
    if umap is not None:
        cells["umap_x"] = np.round(umap[:, 0].astype(float), 6)
        cells["umap_y"] = np.round(umap[:, 1].astype(float), 6)
    else:
        cells["umap_x"] = np.nan
        cells["umap_y"] = np.nan
    return cells


# ================================================================= generate ====
def generate_report(
    tag: str, outdir: str | Path, bundle_only: bool = False
) -> dict[str, str]:
    """Build public annotation tables and an auditable evidence-data bundle."""
    os.makedirs(outdir, exist_ok=True)
    bundle = os.path.join(outdir, FIGDATA_DIR)
    os.makedirs(bundle, exist_ok=True)

    fr = _load_run_artifacts(tag)
    clustermap = build_clustermap(fr)
    table_b = build_marker_evidence(tag, fr)
    identity_markers = build_identity_marker_dotplot(fr, clustermap, table_b)
    celltype_markers = build_celltype_marker_dotplot(fr, clustermap, table_b)
    cells = build_cells(fr)

    paths = {}

    def _w(df, path):
        df.to_csv(path, index=False, encoding="utf-8")
        paths[os.path.basename(path)] = os.path.abspath(path)

    # ---- figure-data (language-agnostic; read by BOTH plots.py & plots.R) ----------
    # De-prefixed, professional layout: results/<name>/figure_data/<file>.csv. File names
    # are dataset-independent (no <tag>_ prefix) so every dataset's directory looks the same.
    _w(cells, os.path.join(bundle, "cells.csv"))
    figure_clustermap_columns = [
        "cluster",
        "cluster_order",
        "n_cells",
        "cell_type_annotation",
        "cell_ontology",
        "annotation_confidence",
        "resolution_status",
        "annotation_source",
        "cluster_celltype",
    ]
    _w(
        clustermap[figure_clustermap_columns],
        os.path.join(bundle, "clustermap.csv"),
    )
    # Two views of the same identity markers: one row per cluster, and one row per
    # resolved cell type with the clusters sharing a label pooled. The cell-type table
    # additionally carries n_pub, which the support-weighted plot reads.
    _w(
        identity_markers,
        os.path.join(bundle, "dotplot_celltype_markers.csv"),
    )
    _w(
        celltype_markers,
        os.path.join(bundle, "dotplot_celltype_markers_by_celltype.csv"),
    )

    if not bundle_only:
        table_a = build_table_a(tag, clustermap, table_b)
        cluster_evidence = build_cluster_evidence(tag, fr, table_b)
        markers_all, markers_significant = build_marker_lists(tag, clustermap)
        _w(table_a, os.path.join(outdir, "cluster_summary.csv"))
        _w(table_b, os.path.join(outdir, "marker_evidence.csv"))
        _w(markers_all, os.path.join(outdir, "markers_all_by_cluster.csv"))
        _w(
            markers_significant,
            os.path.join(outdir, "markers_significant_by_cluster.csv"),
        )
        evidence_path = write_cluster_evidence(cluster_evidence, outdir)
        paths[os.path.basename(evidence_path)] = evidence_path
        note = write_stats_note(tag, outdir, clustermap, cells, table_a, table_b)
        paths[os.path.basename(note)] = note
        from . import viewer, visualization

        visualization.plot_dataset(tag, str(outdir))
        paths.update(viewer.build_dataset_viewer(tag, outdir))

    print(f"[report] {tag}: wrote reference-free result package -> {outdir}")
    return paths


if __name__ == "__main__":
    tag = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RESULTS_DIR, tag)
    generate_report(tag, outdir)
