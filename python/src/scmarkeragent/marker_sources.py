#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""marker_sources.py -- analysis-consistent literature SOURCE lookup.

Serves DB source sentences for (species, tissue, disease, cell_type, gene) using the
SAME context gate as the annotation analysis (marker_database.context_subset / candidate_scoring.R):
  - species exact,
  - tissue expanded to its same-organ sub-region group (tissue_members pooling),
  - disease name(s) matched case-insensitively (None/"ALL" => no disease narrowing),
  - positive markers + gene_qc already enforced in the export (db_export_sources.R).

Per (cell_type, gene) at the resolved context, up to `k` (default 4) DEDUP source
sentences are returned, ordered by citation support (n_pub desc). Every marker carries
its OWN sources so a per-marker judgment is possible.

Both the judge input builder and the main-flow LLM decision draw from this one
consistent lookup (they may still surface different sentences because they query
different gene sets -- sources need not be shared, only the FILTER must be consistent).
"""

import hashlib
import os
import pandas as pd

from .config import CROSS_SPECIES, DB_CSV, DB_SOURCES, DEFAULTS
from .marker_database import (
    ALL_SPECIES,
    normalize_disease_name,
    normalize_disease_query,
    tissue_members,
)
from .uberon_ontology import norm_name

SOURCES_CSV = DB_SOURCES
MAX_PER_MARKER = int(DEFAULTS["cluster_annotation"]["sources_per_marker"])
SOURCE_ORDER_SEED = str(DEFAULTS["cluster_annotation"]["source_order_seed"])


def order_key(seed: str, pmcid: str, sentence: str) -> str:
    """The tie-break among sentences a marker's publication count cannot separate.

    `n_pub_support` is an aggregate of the (cell type, gene) pair, so every source row of
    one pair carries the SAME count -- sorting by it leaves the rows of a pair completely
    tied, and which of them a reader or an agent sees was decided by the sort's internal
    order. That is not reproducible: it varies with the sort implementation, and it
    differs between this arm and the R one, which is why the R arm used to take the order
    across a language bridge rather than compute it.

    A hash of the citation and the sentence gives an order that is arbitrary in the way a
    random draw is arbitrary, and identical wherever it is computed. `seed` selects which
    arbitrary order, so the same run can be repeated under a different draw to show
    whether an annotation depended on which equally-corroborated sentences it happened to
    be shown.
    """
    payload = f"{seed}\x1f{pmcid}\x1f{sentence}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _norm_disease(disease):
    """Mirror marker_database.context_subset: str|list|None; None or 'ALL' => no narrowing."""
    return normalize_disease_query(disease)


class SourceDB:
    """Lazy, species-scoped loader over the curated source export."""

    def __init__(self, csv=SOURCES_CSV):
        self.csv = csv
        self._by_species = {}

    def _load_species(self, species):
        if species in self._by_species:
            return self._by_species[species]
        if not os.path.exists(self.csv):
            raise FileNotFoundError(
                f"curated source export not found: {self.csv}. Point SCMA_RESOURCE_DIR at a "
                "resource bundle containing it."
            )
        cols = [
            "species",
            "tissue_type",
            "disease_normalized",
            "cell_type",
            "gene_symbol",
            "n_pub_support",
            "pmcid",
            "pmid",
            "source",
        ]
        it = pd.read_csv(
            self.csv, usecols=cols, dtype=str, keep_default_na=False, chunksize=200_000
        )
        parts = [c[c["species"] == species] for c in it]
        df = (
            pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
        )
        df["n_pub_support"] = (
            pd.to_numeric(df["n_pub_support"], errors="coerce").fillna(1).astype(int)
        )
        # `norm_name` on BOTH sides, as the marker filter does. Lower-casing only the
        # database side left the comparison asymmetric: the query names arrive normalized,
        # so a curated tissue written "epiblast (generic)" could never match a query that
        # resolved to "epiblast". Three of the 371 curated tissue names are written that
        # way, none of them a tissue this project measures, so this closes the gap without
        # moving any result -- but the two filters now agree by construction.
        df["_tl"] = df["tissue_type"].map(norm_name)
        df["_dl"] = df["disease_normalized"].map(normalize_disease_name)
        df["_gu"] = df["gene_symbol"].str.upper()
        self._by_species[species] = df
        return df

    def context(self, species, tissue, disease, seed=SOURCE_ORDER_SEED, cross_species=None):
        """Return an analysis-consistent context slice (a ContextSources helper).

        `cross_species=None` follows the run configuration
        (`features.cross_species_markers`), so the sentence lookup tracks the SAME switch
        that pools the marker menu (marker_database.context_subset). When the switch is
        ON, the OTHER species' source rows at the SAME tissue + disease are pooled in,
        their gene symbols ortholog-mapped into `species`' symbol space by the same
        OrthoMap the marker pooling uses -- so a translated marker row carries its donor
        species' sentences instead of an empty lookup. OFF => identical to the
        within-species path.
        """
        if cross_species is None:
            cross_species = CROSS_SPECIES
        ts = {norm_name(t) for t in tissue_members(tissue)}
        dl = _norm_disease(disease)

        def scope(frame):
            mask = frame["_tl"].isin(ts)
            if dl is not None:
                mask = mask & frame["_dl"].isin(dl)
            return frame.loc[mask]

        sub = scope(self._load_species(species))
        if cross_species:
            from .ortho_map import OrthoMap

            om = OrthoMap()
            # Native rows first: where a donor row duplicates a native sentence for the
            # same (cell type, gene), the dedup in ContextSources keeps the native one.
            parts = [sub]
            for sp in ALL_SPECIES:
                if sp == species:
                    continue
                donor = scope(self._load_species(sp))
                if donor.empty:
                    continue
                pairs = om.pairs(sp, species)  # src_sym(sp, upper) -> tgt_sym(query)
                if pairs is None or pairs.empty:
                    continue
                donor = donor.rename(columns={"_gu": "src_sym"}).merge(
                    pairs, on="src_sym", how="inner"
                )
                if donor.empty:
                    continue
                donor["gene_symbol"] = donor["tgt_sym"]
                donor["_gu"] = donor["tgt_sym"].str.upper()
                parts.append(donor[sub.columns])
            if len(parts) > 1:
                sub = pd.concat(parts, ignore_index=True)
        return ContextSources(sub, seed=seed)


class ContextSources:
    """A (species,tissue,disease)-scoped slice; serves <=k dedup sources per (cell,gene)."""

    def __init__(self, sub, seed=SOURCE_ORDER_SEED):
        # Ordering happens per (cell type, gene) in `_record_index`, where the tie among
        # a pair's equally-corroborated rows can actually be broken. Sorting the whole
        # slice here by publication count alone did nothing to that tie: it is one number
        # per pair, so the rows of a pair stayed in whatever order the sort left them.
        self.sub = sub
        self.seed = str(seed)
        self._by_cell_gene = None
        self._records = None

    def _index(self):
        """Sentences alone, in the same order as the records they come from."""
        if self._by_cell_gene is None:
            self._by_cell_gene = {
                key: [item["sentence"] for item in bucket]
                for key, bucket in self._record_index().items()
            }
        return self._by_cell_gene

    def _record_index(self):
        """{(cell_type, GENE): [{pmid, pmcid, sentence}]} deduplicated on the sentence.

        A judgement about whether a source really establishes an identity needs the
        citation alongside the sentence, so the identifiers travel with the text rather
        than being dropped at the lookup boundary.

        The list is the pair's WHOLE ordered evidence, not a shortlist: publication count
        descending, then `order_key` within a count. Callers take a prefix of it, and a
        caller that has already taken one can come back for the next -- which is only
        meaningful because the order is fixed rather than incidental.
        """
        if self._records is None:
            index: dict[tuple[str, str], list[dict[str, str]]] = {}
            for cell_type, gene, sentence, pmid, pmcid, n_pub in zip(
                self.sub["cell_type"],
                self.sub["_gu"],
                self.sub["source"],
                self.sub["pmid"],
                self.sub["pmcid"],
                self.sub["n_pub_support"],
            ):
                bucket = index.setdefault((str(cell_type), str(gene)), [])
                text = str(sentence)
                if any(item["sentence"] == text for item in bucket):
                    continue
                bucket.append(
                    {
                        "sentence": text,
                        "pmid": str(pmid),
                        "pmcid": str(pmcid),
                        "_n_pub": int(n_pub),
                    }
                )
            for bucket in index.values():
                bucket.sort(
                    key=lambda item: (
                        -item["_n_pub"],
                        order_key(self.seed, item["pmcid"], item["sentence"]),
                    )
                )
            self._records = index
        return self._records

    def for_marker(self, cell_type, gene, k=MAX_PER_MARKER):
        """<=k dedup source sentences backing (cell_type, gene) at this context."""
        srcs = self._index().get((cell_type, str(gene).upper()), [])
        return srcs[:k]

    def all_records(self, cell_type, gene):
        """Every dedup record backing (cell_type, gene), in the served order."""
        bucket = self._record_index().get((cell_type, str(gene).upper()), [])
        return [
            {"sentence": item["sentence"], "pmid": item["pmid"], "pmcid": item["pmcid"]}
            for item in bucket
        ]

    def records_for_marker(self, cell_type, gene, k=MAX_PER_MARKER):
        """<=k dedup {sentence, pmid, pmcid} records backing (cell_type, gene)."""
        return self.all_records(cell_type, gene)[:k]


class SourceServer:
    """Serves a cluster's source sentences without handing back what it already gave.

    A marker's sentences are not ranked evidence to be truncated at the best few: a pair's
    rows share one publication count, so the first three are simply three of however many
    exist -- sometimes three of several hundred. Fixing the agent's view at those three
    means a claim that turns on the literature is settled by an arbitrary draw, and the
    agent has no way to say that the draw did not settle it.

    So the sentences are a stream. Each request returns the next batch in the fixed order,
    the ones already delivered -- including those that arrived unrequested with the
    opening packet -- are skipped, and the reply says whether anything is left. The cap is
    on how many times one pair can be drawn from, not on what exists: a runaway loop is a
    cost problem, but a hidden second half of the literature is an evidence problem.
    """

    def __init__(self, context, batch=MAX_PER_MARKER, max_batches=0):
        self.context = context
        self.batch = int(batch)
        self.max_batches = int(max_batches)
        self._served: dict[tuple[str, str], set[str]] = {}
        self._draws: dict[tuple[str, str], int] = {}

    def _key(self, cell_type, gene):
        return (str(cell_type), str(gene).upper())

    def note_delivered(self, cell_type, gene, records):
        """Record sentences the agent received without asking, so a request skips them."""
        served = self._served.setdefault(self._key(cell_type, gene), set())
        for record in records:
            served.add(str(record.get("sentence", "")))

    def opening(self, cell_type, gene, k=None):
        """The batch that ships with the packet, marked as delivered."""
        records = self.context.records_for_marker(
            cell_type, gene, k=self.batch if k is None else int(k)
        )
        self.note_delivered(cell_type, gene, records)
        return records

    def take(self, cell_type, gene):
        """The next batch for one pair, with what is left and why it stopped."""
        key = self._key(cell_type, gene)
        every = self.context.all_records(*key)
        served = self._served.setdefault(key, set())
        drawn = self._draws.get(key, 0)
        remaining = [item for item in every if item["sentence"] not in served]

        if self.max_batches and drawn >= self.max_batches:
            return {
                "sources": [],
                "remaining": len(remaining),
                "exhausted": not remaining,
                "limit_reached": True,
            }

        batch = remaining[: self.batch]
        self.note_delivered(cell_type, gene, batch)
        self._draws[key] = drawn + 1
        return {
            "sources": batch,
            "remaining": len(remaining) - len(batch),
            "exhausted": len(remaining) <= len(batch),
            "limit_reached": False,
        }


def ordered_context_records(species, tissue, disease, seed=SOURCE_ORDER_SEED):
    """Language-neutral source-order bridge for the R arm.

    Emits each pair's deduplicated records in the served order, so the arm taking a
    prefix of them takes the same sentences this one would. The order is now `order_key`
    rather than a sort's internal arrangement, which is what makes it reproducible in the
    other language rather than only transferable to it -- this bridge exists to be
    deleted once that arm computes the order itself.
    """
    context = SourceDB().context(species, tissue, disease, seed=seed)
    records = []
    for (cell_type, gene), bucket in context._record_index().items():
        for item in bucket:
            records.append(
                {
                    "cell_type": str(cell_type),
                    "gene": str(gene),
                    "source": item["sentence"],
                    "pmid": item["pmid"],
                    "pmcid": item["pmcid"],
                }
            )
    return records


def ordered_negative_records(species, tissue, min_pub=2):
    """Exact Python Evidence-Constrained LLM Refinement negative-marker ordering for the R bridge."""
    if not os.path.exists(DB_CSV):
        return []
    data = pd.read_csv(DB_CSV, dtype=str, keep_default_na=False)
    data = data[
        (data["marker_polarity"] == "negative")
        & (data["gene_qc_pass"].str.upper() == "TRUE")
    ].copy()
    data["npub"] = (
        pd.to_numeric(data["n_pub_support"], errors="coerce").fillna(0).astype(int)
    )
    tissues = {value.lower() for value in tissue_members(tissue)}
    data = data[
        (data["species"] == species)
        & (data["tissue_type"].str.lower().isin(tissues))
        & (data["npub"] >= int(min_pub))
    ]
    output = []
    for cell_type, group in data.groupby("cell_type"):
        ordered = group.sort_values("npub", ascending=False).drop_duplicates(
            "gene_symbol"
        )
        for row in ordered.itertuples(index=False):
            output.append(
                {
                    "cell_type": str(cell_type),
                    "gene": str(row.gene_symbol),
                    "npub": int(row.npub),
                }
            )
    return output


if __name__ == "__main__":
    db = SourceDB()
    for sp, ti, dis, ct in [
        ("Human", "lung", "Normal", "club cell"),
        (
            "Mouse",
            "pancreas",
            ["Normal", "type 2 diabetes mellitus"],
            "type B pancreatic cell",
        ),
        ("Rat", "lung", "Normal", "T cell"),
    ]:
        cx = db.context(sp, ti, dis)
        print(f"\n=== {sp}/{ti}/{dis} :: {ct} ===  context rows={len(cx.sub):,}")
        idx = cx._index()
        genes = [g for (c, g) in idx if c == ct][:5]
        for g in genes:
            ss = cx.for_marker(ct, g)
            print(f"  {g}: {len(ss)} src | e.g. {ss[0][:90] if ss else '(none)'}")
