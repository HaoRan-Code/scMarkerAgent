#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from .marker_sources import CompositeSources, SourceDB
from .ortho_map import OrthoMap

_D = DEFAULTS["output"]["decimals"]
ND_CONF = int(_D["confidence"])
ND_FRAC = int(_D["fraction"])
ND_S = int(_D["specificity"])

TABLE_A_COLS = list(OUTPUT_SCHEMA["cluster_summary"]["columns"])
TABLE_B_COLS = list(OUTPUT_SCHEMA["marker_evidence"]["columns"])
ND_LOGFC = int(_D["logfc"])
SUMMARY_MARKER_CAP = 12

def _cluster_sort_key(c):
    s = str(c)
    return (0, int(s)) if s.lstrip("-").isdigit() else (1, 0, s)


def _fmt_frac(x):
    return f"{float(x):.{ND_FRAC}f}"


def _fmt_pvalue(value):
    if value is None:
        return NOT_AVAILABLE
    try:
        number = float(value)
    except (TypeError, ValueError):
        return NOT_AVAILABLE
    return NOT_AVAILABLE if not np.isfinite(number) else repr(number)


def _fmt(value, decimals):
    if value is None:
        return NOT_AVAILABLE
    try:
        number = float(value)
    except (TypeError, ValueError):
        return NOT_AVAILABLE
    return NOT_AVAILABLE if not np.isfinite(number) else f"{number:.{decimals}f}"


def _fmt_alternative_candidates(names):
    parts = [str(name) for name in names or [] if str(name).strip()]
    return "; ".join(parts) if parts else NOT_AVAILABLE


def _fmt_claim_warnings(lines):
    parts = [str(line) for line in lines or [] if str(line).strip()]
    return " || ".join(parts) if parts else NOT_AVAILABLE


def _json_safe(value):
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


def build_clustermap(fr):
    res = fr["res"]
    cl_id = fr["scoring"].get("candidate_cl_id", {})
    rows = []
    for cluster in sorted(res.keys(), key=_cluster_sort_key):
        record = res[cluster]
        label = str(record.get("annotation") or NOT_AVAILABLE)
        subtype = str(record.get("subtype") or "")
        primary = subtype or label
        rows.append(
            dict(
                cluster=str(cluster),
                n_cells=int(fr["scoring"]["clusters"][str(cluster)]["n_cells"]),
                cell_type_annotation=label,
                cell_subtype_annotation=na_display(record.get("subtype")),
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
                evidence_review_rounds=len(
                    (record.get("review") or {}).get("rounds") or []
                ),
                delivered_tier=int(record.get("delivered_tier") or 0),
                arbitrated=bool(record.get("arbitration")),
            )
        )
    clustermap = pd.DataFrame(rows)
    clustermap.insert(1, "cluster_order", range(len(clustermap)))
    return clustermap


def _by_evidence(frame):
    ranked = frame
    if "marker_polarity" in ranked.columns:
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
        markers = "; ".join(
            f"{row.gene}|{row.candidate_annotation}|{row.cluster_detection_fraction}"
            for row in rows
        )
    else:
        markers = "; ".join(
            f"{row.gene}({row.cluster_detection_fraction})" for row in rows
        )
    pmcids = "; ".join(
        dict.fromkeys(
            str(row.pmcid)
            for row in rows
            if str(getattr(row, "pmcid", "")) not in {"", NOT_AVAILABLE}
        )
    )
    return markers, pmcids


def _cooccurring_names(alternative_candidates) -> set[str]:
    text = str(alternative_candidates or "").strip()
    if not text or text == NOT_AVAILABLE:
        return set()
    return {part.strip() for part in text.split(";") if part.strip() and part.strip() != NOT_AVAILABLE}


def build_table_a(tag, clustermap, marker_evidence):
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
                "evidence_review_rounds": row["evidence_review_rounds"],
                "delivered_tier": row["delivered_tier"],
                "arbitrated": row["arbitrated"],
            }
        )
    return pd.DataFrame(out, columns=TABLE_A_COLS)


def build_cluster_evidence(tag, fr, marker_evidence):
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
                "candidate_facts": dict(result.get("candidate_facts") or {}),
            },
            "audit": {
                "annotator_turns": int(result.get("turns", 0) or 0),
                "tool_calls": list(result.get("tool_calls") or []),
                "evidence_review": dict(result.get("review") or {}),
                "review_flags": list(result.get("review_flags") or []),
                "program_reading": dict(result.get("program_reading") or {}),
                "definer_audit": dict(result.get("definer_audit") or {}),
                "contested_rivals": list(result.get("contested_rivals") or []),
                "contested": list(result.get("contested") or []),
                "confidence_cap": result.get("confidence_cap") or {},
                "structural_notes": list(result.get("structural_notes") or []),
                "consistency_note": result.get("consistency_note") or {},
                "reopened_by_consistency": result.get("reopened_by_consistency") or {},
                "arbitration": result.get("arbitration") or {},
                "delivered_tier": int(result.get("delivered_tier") or 0),
                "tiers_available": int(result.get("tiers_available") or 0),
                "resolution_detail": str(result.get("resolution_detail") or ""),
                "candidate_pool_size": int(len(pool)),
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
    upper = gene.upper()
    evidence = entry.get("claim_evidence") or {}
    if str(evidence.get("decisive_gene", "")).upper() == upper:
        return "decisive"
    if upper in {
        str(value).upper() for value in result.get("support_markers", []) or []
    }:
        return "support"
    return "panel"


def _borrowed_context_names(fr) -> tuple[set[str], dict[str, set[str]]]:
    names: set[str] = set()
    donor_species_by_name: dict[str, set[str]] = {}
    for result in fr["res"].values():
        for entry in result.get("candidate_entries") or []:
            borrowed = entry.get("borrowed_context")
            if not borrowed:
                continue
            name = str(entry["cell_type"])
            names.add(name)
            donor = borrowed.get("donor_species")
            if donor:
                donor_species_by_name.setdefault(name, set()).add(str(donor))
    return names, donor_species_by_name


def build_marker_evidence(tag, fr):
    context = fr["scoring"]["context"]
    source_db = SourceDB()
    native_sources = source_db.context(
        context["species"], context["tissue"], context["disease"]
    )
    borrowed_names, donor_species_by_name = _borrowed_context_names(fr)
    if borrowed_names:
        ortho = OrthoMap() if donor_species_by_name else None
        extra = source_db.cell_types_across_tissues(
            context["species"], borrowed_names, context["disease"],
            donor_species_by_name=donor_species_by_name, ortho=ortho,
        )
        sources = CompositeSources(native_sources, extra, borrowed_names)
    else:
        sources = native_sources
    specificity = fr["scoring"]["marker_specificity"]
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
        assigned = {
            str(result.get("annotation") or ""),
            str(result.get("subtype") or ""),
        } - {""}
        reported = {str(name) for name in (result.get("possible_refinements") or [])}
        for entry in result.get("candidate_entries") or []:
            candidate = str(entry["cell_type"])
            is_selected = candidate in assigned
            claimed = is_selected or bool(str(entry.get("claim_role") or ""))
            if not claimed and candidate not in reported:
                continue
            for marker in entry.get("decisive_marker_measurements") or []:
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
                        "auc": _fmt(marker.get("auc"), ND_FRAC),
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
    if marker_evidence.empty:
        return [], {}
    panel = marker_evidence[
        marker_evidence["is_selected_annotation"].astype(str).str.lower() == "true"
    ].copy()
    if "marker_polarity" in panel.columns:
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
    dm = fr["dm"]
    menu = list(dm["menu_genes"])
    gene_index = {g.upper(): i for i, g in enumerate(menu)}
    cell_cluster = dm["meta"]["cluster"].astype(str).values

    order = clustermap.sort_values("cluster_order")
    label_of = {
        str(row.cluster): str(row.primary_annotation) or str(row.cluster)
        for row in order.itertuples(index=False)
    }
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
        evidence_gate.sig_pass(lfc, pin, padj, auc)
        for lfc, pin, padj, auc in zip(
            markers["avg_log2FC"].astype(float),
            markers["pct_in"].astype(float) / 100.0,
            markers["padj"].astype(float),
            markers["auc"].astype(float),
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


def generate_report(
    tag: str, outdir: str | Path, bundle_only: bool = False
) -> dict[str, str]:
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
