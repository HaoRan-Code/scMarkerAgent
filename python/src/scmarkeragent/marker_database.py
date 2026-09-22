#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from .config import (
    DB_CSV,
    OBO_UBERON,
    TISSUE_ROOT,
    EXCLUDE_IN_VITRO,
    RU_MIN_ELIG_GENES,
    MIN_CORROBORATING_PUBLICATIONS,
    CORROBORATING_TIERS,
    NOT_AVAILABLE,
)

ALL_SPECIES = ("Human", "Mouse", "Rat")


def normalize_disease_name(value):
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().lower()


def normalize_disease_query(disease):
    if disease is None:
        return None
    values = [disease] if isinstance(disease, str) else list(disease)
    terms = {normalize_disease_name(value) for value in values}
    terms.discard("")
    return None if not terms or "all" in terms else terms


_UB = None

TISSUE_CONTEXT_SEP = ";"


def tissue_list(tissue):
    if tissue is None:
        return []
    values = [tissue] if isinstance(tissue, str) else list(tissue)
    out = []
    for value in values:
        name = str(value).strip()
        if name and name not in out:
            out.append(name)
    return out


def parse_tissue_arg(value):
    text = str(value or "")
    if "|" not in text:
        return text
    names = tissue_list(text.split("|"))
    if not names:
        raise ValueError("--tissue lists no tissue name")
    return names[0] if len(names) == 1 else names


def _uberon():
    global _UB
    if _UB is None:
        from .uberon_ontology import UberonOntology

        _UB = UberonOntology(OBO_UBERON)
    return _UB


def _one_tissue_closure(tissue, tissue_root=None):
    from .uberon_ontology import norm_name

    ub = _uberon()
    raw = TISSUE_ROOT if tissue_root is None else tissue_root
    tr = str(raw or "self").strip().lower()
    ids, names = set(), {norm_name(tissue)}
    rid = ub.resolve(tissue)
    if rid:
        ids.add(rid)
    if tr in ("exact", "none", "off"):
        return ids, names
    root = tissue if tr in ("self", "") else raw
    r = ub.resolve(root)
    if r:
        cids = ub.closure_self(r)
        ids |= cids
        names |= {norm_name(ub.name_of(x)) for x in cids}
    return ids, names


def tissue_closures(tissue, tissue_root=None):
    return {name: _one_tissue_closure(name, tissue_root) for name in tissue_list(tissue)}


def tissue_closure(tissue, tissue_root=None):
    ids, names = set(), set()
    for cids, cnames in tissue_closures(tissue, tissue_root).values():
        ids |= cids
        names |= cnames
    return ids, names


def tissue_members(tissue):
    return tissue_closure(tissue)[1]


def merge_tissue_contexts(labels, order=None):
    seen = []
    for label in labels:
        if label is None or (isinstance(label, float) and pd.isna(label)):
            continue
        for name in str(label).split(TISSUE_CONTEXT_SEP):
            if name and name not in seen:
                seen.append(name)
    if order:
        rank = {name: index for index, name in enumerate(order)}
        seen.sort(key=lambda name: (rank.get(name, len(rank)), name))
    return TISSUE_CONTEXT_SEP.join(seen)


class MarkerDatabase:
    def __init__(self, csv=DB_CSV):
        self.csv = csv
        self._df = None

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
        df["n_pub_support"] = (
            pd.to_numeric(df["n_pub_support"], errors="coerce").fillna(1).astype(int)
        )
        from .uberon_ontology import norm_name

        df["_tt_norm"] = df["tissue_type"].astype(str).map(norm_name)
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

    @staticmethod
    def _derive(sub):
        sub["g"] = sub["gene_symbol"].astype(str)
        sub["cell"] = sub["cell_type"].astype(str)
        canon = sub["cell_type_canonical_cl_id"]
        clid = sub["cell_type_cl_id"]
        sub["clid"] = canon.where(canon.notna() & (canon != ""), clid)
        sub["src_species"] = sub["species"]
        return sub

    def context_subset(self, species, tissue, disease, cross_species=False):
        df = self.load()
        disease_terms = normalize_disease_query(disease)
        no_disease = disease_terms is None
        qc_ok = (df["gene_qc_pass"] == "TRUE") | (df["gene_qc_pass"].isna())
        g_ok = df["gene_symbol"].notna() & (
            ~df["gene_symbol"].str.lower().isin(["", "unknown"])
        )
        c_ok = df["cell_type"].notna() & (df["cell_type"] != "")
        dmask = True if no_disease else df["_disease_norm"].isin(disease_terms)
        pol = df["marker_polarity"].isin(["positive", "negative"])
        tmask, context_label = self._tissue_masks(df, tissue)
        if EXCLUDE_IN_VITRO:
            tmask = tmask & (~df["is_in_vitro"])
        base = tmask & dmask & qc_ok & g_ok & c_ok & pol
        sub = df.loc[(df["species"] == species) & base].copy()
        sub["tissue_context"] = context_label.loc[sub.index]
        sub = self._derive(sub)
        if cross_species:
            sub = self._pool_cross_species(sub, df, base, species, context_label)
        return sub

    @staticmethod
    def _tissue_masks(df, tissue):
        closures = tissue_closures(tissue)
        names_in_order = list(closures)
        bits = np.zeros(len(df), dtype=np.int64)
        union = np.zeros(len(df), dtype=bool)
        for position, (ids, names) in enumerate(closures.values()):
            mask = df["_tt_norm"].isin(names).to_numpy()
            if "tissue_type_uberon_id" in df.columns:
                mask = mask | df["tissue_type_uberon_id"].isin(ids).to_numpy()
            bits |= mask.astype(np.int64) << position
            union |= mask
        label_of = {
            int(value): TISSUE_CONTEXT_SEP.join(
                name for k, name in enumerate(names_in_order) if (int(value) >> k) & 1
            )
            for value in np.unique(bits[union])
        }
        label = pd.Series(bits, index=df.index).map(label_of).fillna("").astype(str)
        return pd.Series(union, index=df.index), label

    def _pool_cross_species(self, sub, df, base, species, context_label=None):
        from .ortho_map import OrthoMap

        om = OrthoMap()
        parts = [sub]
        for sp in ALL_SPECIES:
            if sp == species:
                continue
            ssub = df.loc[(df["species"] == sp) & base].copy()
            if ssub.empty:
                continue
            if context_label is not None:
                ssub["tissue_context"] = context_label.loc[ssub.index]
            pr = om.pairs(sp, species)
            if pr is None or pr.empty:
                continue
            ssub["src_sym"] = ssub["gene_symbol"].str.upper()
            mg = ssub.merge(pr, on="src_sym", how="inner")
            if mg.empty:
                continue
            mg["gene_symbol"] = mg["tgt_sym"]
            mg = self._derive(mg)
            mg["src_species"] = sp
            parts.append(mg[sub.columns])
        if len(parts) == 1:
            return sub
        out = pd.concat(parts, ignore_index=True)
        out = out.drop_duplicates(
            subset=["cell", "g", "marker_polarity", "n_pub_support", "confidence_tier"]
        )
        return out.reset_index(drop=True)

    @staticmethod
    def candidate_cl_id(sub):
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

    @staticmethod
    def panel(sub, tissue_contexts=None):
        columns = [
            "cell",
            "g",
            "marker_polarity",
            "n_pub_support",
            "confidence_tier",
            "is_in_vitro",
        ]
        has_context = "tissue_context" in sub.columns
        if has_context:
            columns.append("tissue_context")
        panel = (
            sub[columns]
            .rename(columns={"n_pub_support": "n_pub", "confidence_tier": "tier"})
            .copy()
        )
        panel["n_pub"] = panel["n_pub"].fillna(1).astype(int)
        tier_rank = {"high": 0, "medium": 1, "low": 2}
        panel["_tr"] = (
            panel["tier"].astype(str).str.lower().map(tier_rank).fillna(3).astype(int)
        )
        panel = panel.sort_values(["n_pub", "_tr"], ascending=[False, True], kind="stable")
        merged_context = None
        if has_context:
            key = ["cell", "g", "marker_polarity"]
            multi = panel.groupby(key, sort=False)["tissue_context"].nunique()
            if (multi > 1).any():
                merged_context = (
                    panel.groupby(key, sort=False)["tissue_context"]
                    .agg(lambda values: merge_tissue_contexts(values, tissue_contexts))
                )
        panel = (
            panel.drop_duplicates(["cell", "g", "marker_polarity"], keep="first")
            .drop(columns="_tr")
            .reset_index(drop=True)
        )
        if merged_context is not None:
            index = pd.MultiIndex.from_frame(panel[["cell", "g", "marker_polarity"]])
            panel["tissue_context"] = merged_context.reindex(index).to_numpy()
        return panel

    @staticmethod
    def eligible_candidates(panel, measured_genes, corrob_only=True):
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

    @staticmethod
    def species_positive_rows(df, species):
        qc_ok = (df["gene_qc_pass"] == "TRUE") | (df["gene_qc_pass"].isna())
        g_ok = df["gene_symbol"].notna() & (
            ~df["gene_symbol"].str.lower().isin(["", "unknown"])
        )
        mask = (
            (df["species"] == species)
            & (df["marker_polarity"] == "positive")
            & qc_ok
            & g_ok
            & df["cell_type"].notna()
        )
        if EXCLUDE_IN_VITRO:
            mask = mask & (~df["is_in_vitro"])
        return df.loc[mask, ["cell_type", "gene_symbol", "n_pub_support",
                             "confidence_tier", "is_recommended_marker", "tissue_type"]]

    @staticmethod
    def species_definers(rows, names, k=8):
        sub = rows[rows["cell_type"].isin(set(names))]
        if sub.empty:
            return {}
        rec = sub["is_recommended_marker"].astype(str).str.upper() == "TRUE"
        tier_rank = {"high": 0, "medium": 1, "low": 2}
        sub = sub.assign(
            rec_flag=rec.astype(int),
            tier_rank=sub["confidence_tier"].astype(str).str.lower().map(tier_rank).fillna(3).astype(int),
        )
        grp = (
            sub.groupby(["cell_type", "gene_symbol"], sort=False)
            .agg(
                n_pub=("n_pub_support", "max"),
                rec_flag=("rec_flag", "max"),
                tier_rank=("tier_rank", "min"),
                tissue_contexts=("tissue_type", "nunique"),
            )
            .reset_index()
        )
        inv_tier = {0: "high", 1: "medium", 2: "low", 3: NOT_AVAILABLE}
        out = {}
        for name, block in grp.groupby("cell_type", sort=False):
            block = block.sort_values(
                ["n_pub", "rec_flag", "gene_symbol"], ascending=[False, False, True], kind="stable"
            )
            head = block.head(k)
            extra = block[(block["rec_flag"] == 1) & (~block.index.isin(head.index))].head(4)
            chosen = pd.concat([head, extra])
            out[str(name)] = [
                {
                    "gene": str(r.gene_symbol),
                    "n_pub": int(r.n_pub),
                    "tier": inv_tier.get(int(r.tier_rank), NOT_AVAILABLE),
                    "recommended": bool(r.rec_flag),
                    "tissue_contexts": int(r.tissue_contexts),
                }
                for r in chosen.itertuples(index=False)
            ]
        return out

    @staticmethod
    def species_gene_claimants(rows):
        grp = (
            rows.groupby([rows["gene_symbol"].str.upper(), "cell_type"], sort=False)["n_pub_support"]
            .max()
            .reset_index()
        )
        out: dict[str, list[tuple[str, int]]] = {}
        for r in grp.itertuples(index=False):
            out.setdefault(str(r[0]), []).append((str(r.cell_type), int(r.n_pub_support)))
        for gene in out:
            out[gene].sort(key=lambda item: (-item[1], item[0]))
        return out

    @staticmethod
    def known_genes(df):
        g = df["gene_symbol"].dropna().astype(str)
        return set(g.str.upper())

    def menu_genes(self, species, tissue, disease, cross_species=False, polarity="positive"):
        sub = self.context_subset(species, tissue, disease, cross_species=cross_species)
        if polarity != "any":
            sub = sub[sub["marker_polarity"] == polarity]
        return sorted(set(sub["g"].dropna()))
