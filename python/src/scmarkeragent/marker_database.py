#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Curated marker database access, mirrored by the R arm's candidate_scoring.R.

Owns context gating (species, tissue, normalized disease, gene QC, valid gene and
cell type, marker polarity), the per-candidate marker panel, and the candidate
eligibility gates: enough measured positive genes, the in-vitro origin gate, the
cross-species gate, and publication corroboration.

A candidate IS its curated free-text cell-type string for a species. The curated CL id
is carried only as a display passthrough (`candidate_cl_id`) and never resolved,
expanded or compared.

Reads the curated marker export of the authoritative database, whose columns are:
  species, tissue_type, disease_normalized, gene_symbol, gene_qc_pass,
  cell_type, cell_type_cl_id, cell_type_canonical_cl_id,
  marker_polarity, n_pub_support, confidence_tier, is_recommended_marker
Scoring consumes only n_pub_support and confidence_tier as reliability signals;
is_recommended_marker is curation metadata that neither arm reads.
"""

import pandas as pd

from .config import (
    DB_CSV,
    OBO_UBERON,
    TISSUE_ROOT,
    EXCLUDE_IN_VITRO,
    RU_MIN_ELIG_GENES,
    MIN_CORROBORATING_PUBLICATIONS,
    CORROBORATING_TIERS,
)

ALL_SPECIES = ("Human", "Mouse", "Rat")


def normalize_disease_name(value):
    """Canonical comparison key for disease labels.

    Policy is deliberately narrow: trim surrounding whitespace and compare
    case-insensitively.  Do not perform substring, fuzzy, synonym, or ontology
    expansion; those operations can silently broaden the marker context.
    """
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().lower()


def normalize_disease_query(disease):
    """Return normalized exact-match terms, or None for the explicit ALL mode."""
    if disease is None:
        return None
    values = [disease] if isinstance(disease, str) else list(disease)
    terms = {normalize_disease_name(value) for value in values}
    terms.discard("")
    return None if not terms or "all" in terms else terms


# ---- tissue subsumption via UBERON (replaces the old hand-written tissue grouping) ---
# A registry `tissue` is resolved to a UBERON id; the per-dataset menu pools every
# DB row whose tissue_type is that tissue OR (per tissue_root) a UBERON descendant
# of it. This generalizes the old intestine-only hand list to ALL organs WITHOUT
# cross-ORGAN borrowing: is_a/part_of stay within the organ (digestive_tract is NOT
# an ancestor of intestine, so stomach/esophagus are never pulled in). Matching is
# by tissue_type_uberon_id in the closure OR normalized tissue_type name in the
# closure names, so the ~10% of rows lacking a uberon_id are still caught by name
# (no coverage regression vs the old exact/hand-list gate).
_UB = None


def _uberon():
    global _UB
    if _UB is None:
        from .uberon_ontology import UberonOntology

        _UB = UberonOntology(OBO_UBERON)
    return _UB


def tissue_closure(tissue, tissue_root=None):
    """(uberon_id set, normalized-name set) of DB tissues to pool for a registry
    `tissue`, per the tissue_root switch (default from config: 'self' = subsume the
    dataset's own tissue; 'exact'/'none' = no subsumption; else = explicit root)."""
    from .uberon_ontology import norm_name

    ub = _uberon()
    raw = TISSUE_ROOT if tissue_root is None else tissue_root
    tr = str(raw or "self").strip().lower()
    ids, names = set(), {norm_name(tissue)}
    rid = ub.resolve(tissue)
    if rid:
        ids.add(rid)
    if tr in ("exact", "none", "off"):
        return ids, names  # exact tissue match: no subsumption
    root = tissue if tr in ("self", "") else raw
    r = ub.resolve(root)
    if r:
        cids = ub.closure_self(r)
        ids |= cids
        names |= {norm_name(ub.name_of(x)) for x in cids}
    return ids, names


def tissue_members(tissue):
    """Back-compat: the set of DB tissue_type NAMES to pool (now UBERON-derived)."""
    return tissue_closure(tissue)[1]


class MarkerDatabase:
    def __init__(self, csv=DB_CSV):
        self.csv = csv
        self._df = None

    # ---- load -----------------------------------------------------------------
    def load(self):
        if self._df is not None:
            return self._df
        dt = {
            "species": "string",
            "tissue_type": "string",
            "disease_normalized": "string",
            "gene_symbol": "string",
            "gene_qc_pass": "string",
            "cell_type": "string",
            "cell_type_cl_id": "string",
            "cell_type_canonical_cl_id": "string",
            "marker_polarity": "string",
            "confidence_tier": "string",
        }
        df = pd.read_csv(self.csv, dtype=dt, low_memory=False)
        # n_pub_support: db_export already coerced NA->1 (integer); be defensive
        df["n_pub_support"] = (
            pd.to_numeric(df["n_pub_support"], errors="coerce").fillna(1).astype(int)
        )
        # normalized tissue_type for UBERON name matching (cache once)
        from .uberon_ontology import norm_name

        df["_tt_norm"] = df["tissue_type"].astype(str).map(norm_name)
        # Disease matching is normalized exact equality on BOTH sides.  Keeping
        # the key in the loaded frame makes every query use the same policy.
        df["_disease_norm"] = df["disease_normalized"].map(normalize_disease_name)
        if "is_in_vitro" not in df.columns:
            raise ValueError(
                "curated marker resource lacks required is_in_vitro column"
            )
        df["is_in_vitro"] = (
            df["is_in_vitro"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(["true", "1", "yes", "t"])
        )
        self._df = df
        return df

    # ---- context subset (candidate_scoring.R sub_db) ------------------------------
    @staticmethod
    def _derive(sub):
        """Add derived columns g(NATIVE gene symbol)/cell/clid/src_species to a context
        slice. Gene casing is preserved (species nomenclature) and matching against the
        measured matrix genes is by EXACT symbol at every comparison site."""
        sub["g"] = sub["gene_symbol"].astype(str)
        sub["cell"] = sub["cell_type"].astype(str)
        canon = sub["cell_type_canonical_cl_id"]
        clid = sub["cell_type_cl_id"]
        sub["clid"] = canon.where(canon.notna() & (canon != ""), clid)
        sub["src_species"] = sub["species"]
        return sub

    def context_subset(self, species, tissue, disease, cross_species=False):
        """Filter to (species, tissue, disease names) + gene_qc + valid gene/cell
        + polarity in {positive,negative}; add derived g (upper), cell, clid.

        disease=None or "ALL" (or a list containing "ALL") = tissue-level Layer1
        (species+tissue, NO disease narrowing), per the confirmed Q6 decision
        that cell identity markers are disease-independent. Otherwise filter to the
        given disease name(s) (e.g. ["Normal"], ["Normal","glioblastoma"]) -- this
        matches candidate_scoring.R's disease_normalized gate for the R datasets.

        cross_species=True: MUTUALLY pool the OTHER species' markers for the SAME
        tissue + SAME disease, mapping their gene symbols into `species`'s symbol
        space via ortholog (ortho_map.OrthoMap). Bidirectional (any query species
        borrows from any other), strictly within identical tissue+disease (no
        cross-condition). OFF by default => identical to the within-species path."""
        df = self.load()
        disease_terms = normalize_disease_query(disease)
        no_disease = disease_terms is None
        # Normalized EXACT equality is the intended policy.  In particular,
        # "glioblastoma" does not include "IDH-wildtype glioblastoma".
        qc_ok = (df["gene_qc_pass"] == "TRUE") | (df["gene_qc_pass"].isna())
        g_ok = df["gene_symbol"].notna() & (
            ~df["gene_symbol"].str.lower().isin(["", "unknown"])
        )
        c_ok = df["cell_type"].notna() & (df["cell_type"] != "")
        dmask = True if no_disease else df["_disease_norm"].isin(disease_terms)
        pol = df["marker_polarity"].isin(["positive", "negative"])
        # context mask WITHOUT species (so the SAME tissue+disease gate applies to
        # every species when pooling cross-species). tissue is expanded via UBERON
        # subsumption (tissue_root switch): a row matches if its tissue_type_uberon_id
        # is in the closure OR its normalized tissue_type name is in the closure names
        # (name fallback keeps the ~10% rows without a uberon_id -> no regression).
        ids, names = tissue_closure(tissue)
        tmask = df["_tt_norm"].isin(names)
        if "tissue_type_uberon_id" in df.columns:
            tmask = tmask | df["tissue_type_uberon_id"].isin(ids)
        if EXCLUDE_IN_VITRO:
            tmask = tmask & (~df["is_in_vitro"])
        base = tmask & dmask & qc_ok & g_ok & c_ok & pol
        sub = self._derive(df.loc[(df["species"] == species) & base].copy())
        if cross_species:
            sub = self._pool_cross_species(sub, df, base, species)
        return sub

    def _pool_cross_species(self, sub, df, base, species):
        """Append ortholog-mapped markers from the OTHER species (same tissue+disease
        via `base`). 1:many orthologs all kept; cell_name/CL id are species-agnostic
        so a candidate (e.g. hepatocyte) merges across species and its panel grows."""
        from .ortho_map import OrthoMap

        om = OrthoMap()
        parts = [sub]
        for sp in ALL_SPECIES:
            if sp == species:
                continue
            ssub = df.loc[(df["species"] == sp) & base].copy()
            if ssub.empty:
                continue
            pr = om.pairs(sp, species)  # src_sym(sp) -> tgt_sym(query)
            if pr is None or pr.empty:
                continue
            ssub["src_sym"] = ssub["gene_symbol"].str.upper()
            mg = ssub.merge(pr, on="src_sym", how="inner")
            if mg.empty:
                continue
            mg["gene_symbol"] = mg["tgt_sym"]  # remap into query symbol space
            mg = self._derive(mg)
            mg["src_species"] = sp
            parts.append(mg[sub.columns])
        if len(parts) == 1:
            return sub
        out = pd.concat(parts, ignore_index=True)
        # drop exact duplicate marker rows (same candidate+gene+polarity+pub+tier)
        out = out.drop_duplicates(
            subset=["cell", "g", "marker_polarity", "n_pub_support", "confidence_tier"]
        )
        return out.reset_index(drop=True)

    # ---- candidate -> curated CL id, for DISPLAY ONLY -------------------------
    @staticmethod
    def candidate_cl_id(sub):
        """{free-text cell type: curated CL id} for the ids the resource happens to carry.

        This is a passthrough for one output column and nothing else. No id is resolved
        against an ontology, expanded to ancestors or descendants, compared between two
        candidates, or used to merge, gate or rank anything: most curated free-text cell
        types have no defensible CL id, so any operation built on these ids would be
        systematically wrong exactly where the free-text vocabulary is richest. A
        candidate missing from this map is reported as not applicable.
        """
        v = sub[sub["clid"].notna() & (sub["clid"] != "")]
        if v.empty:
            return {}
        cnt = (
            v.groupby(["cell", "clid"])
            .size()
            .reset_index(name="N")
            .sort_values(["cell", "N"], ascending=[True, False])
        )
        first = cnt.drop_duplicates("cell", keep="first")
        return dict(zip(first["cell"], first["clid"]))

    # ---- per-candidate panel: exactly one row per (cell, g, marker_polarity) ---------
    @staticmethod
    def panel(sub):
        """One row per (cell, g, marker_polarity).

        Deduplicating on n_pub or tier as well would let the SAME (cell, gene, polarity)
        curated with different citation counts survive as several rows and be
        DOUBLE-COUNTED downstream: in the shortlist evidence mass, in the call gate
        k_markers/k_co, and in the displayed key markers. Collapse to exactly one row per
        (cell, g, marker_polarity), keeping the MAX n_pub and the best tier
        (high > medium > low):
        add a tier-rank column, sort by [n_pub desc, tier_rank asc], keep the first row per key.
        """
        panel = (
            sub[
                [
                    "cell",
                    "g",
                    "marker_polarity",
                    "n_pub_support",
                    "confidence_tier",
                    "is_in_vitro",
                ]
            ]
            .rename(columns={"n_pub_support": "n_pub", "confidence_tier": "tier"})
            .copy()
        )
        panel["n_pub"] = panel["n_pub"].fillna(1).astype(int)
        tier_rank = {"high": 0, "medium": 1, "low": 2}
        panel["_tr"] = (
            panel["tier"].astype(str).str.lower().map(tier_rank).fillna(3).astype(int)
        )
        panel = (
            panel.sort_values(["n_pub", "_tr"], ascending=[False, True], kind="stable")
            .drop_duplicates(["cell", "g", "marker_polarity"], keep="first")
            .drop(columns="_tr")
            .reset_index(drop=True)
        )
        return panel

    # ---- eligible candidate gating (candidate_scoring.R lines 143-155) ------------
    @staticmethod
    def eligible_candidates(panel, measured_genes, corrob_only=True):
        """Return (pos_df, eligible_cells). pos_df is positive markers whose gene is
        measured, restricted to eligible candidates. Gates (verbatim):
          - >=RU_MIN_ELIG_GENES distinct measured positive genes (default 2; was 3),
          - curated in-vitro origin exclusion,
          - corroboration: max n_pub>=2 OR >=1 high/medium-tier marker.
        """
        # EXACT symbol membership: this is the second place a curated symbol meets an
        # uploaded one (the first is the DE menu in preprocessing), and the two arms
        # must gate it the same way. Case-folding here let a mouse panel qualify on a
        # human matrix and vice versa.
        meas = {str(x) for x in measured_genes}
        pos = panel[
            (panel["marker_polarity"] == "positive") & (panel["g"].isin(meas))
        ].copy()
        if EXCLUDE_IN_VITRO:
            pos = pos[~pos["is_in_vitro"]].copy()
        if pos.empty:
            return pos, []
        npg = pos.groupby("cell")["g"].nunique()
        elig = set(npg[npg >= RU_MIN_ELIG_GENES].index)
        if corrob_only:
            grp = pos.groupby("cell")
            mx = grp["n_pub"].max()
            himed = grp["tier"].apply(
                lambda values: int(values.isin(CORROBORATING_TIERS).sum())
            )
            reli = set(mx[(mx >= MIN_CORROBORATING_PUBLICATIONS)].index) | set(
                himed[himed >= 1].index
            )
            elig = elig & reli
        pos = pos[pos["cell"].isin(elig)].copy()
        return pos, sorted(elig)

    # ---- convenience: menu genes for a context (for the prep DE / expression union) ---
    def menu_genes(self, species, tissue, disease, cross_species=False, polarity="positive"):
        """Curated gene symbols for this context at the requested marker polarity.

        Negative markers are measured alongside the positive ones so that a cluster's
        actual detection fraction for a curated negative marker is observable. They
        never enter candidate retrieval or any score; only the positive polarity does.
        """
        sub = self.context_subset(species, tissue, disease, cross_species=cross_species)
        if polarity != "any":
            sub = sub[sub["marker_polarity"] == polarity]
        return sorted(set(sub["g"].dropna()))
