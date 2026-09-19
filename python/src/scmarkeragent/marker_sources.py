#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    payload = f"{seed}\x1f{pmcid}\x1f{sentence}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _norm_disease(disease):
    return normalize_disease_query(disease)


class SourceDB:

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
        df["_tl"] = df["tissue_type"].map(norm_name)
        df["_dl"] = df["disease_normalized"].map(normalize_disease_name)
        df["_gu"] = df["gene_symbol"].str.upper()
        self._by_species[species] = df
        return df

    def context(self, species, tissue, disease, seed=SOURCE_ORDER_SEED, cross_species=None):
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
            parts = [sub]
            for sp in ALL_SPECIES:
                if sp == species:
                    continue
                donor = scope(self._load_species(sp))
                if donor.empty:
                    continue
                pairs = om.pairs(sp, species)
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

    def cell_types_across_tissues(self, species, cell_types, disease, seed=SOURCE_ORDER_SEED,
                                  donor_species_by_name=None, ortho=None):
        names = {str(n) for n in cell_types}
        frame = self._load_species(species)
        mask = frame["cell_type"].isin(names)
        dl = _norm_disease(disease)
        if dl is not None:
            mask = mask & frame["_dl"].isin(dl)
        parts = [frame.loc[mask]]
        if donor_species_by_name and ortho is not None:
            by_donor: dict[str, set[str]] = {}
            for name, donors in donor_species_by_name.items():
                for donor in donors:
                    by_donor.setdefault(str(donor), set()).add(str(name))
            for donor, donor_names in by_donor.items():
                gene_map = ortho.one_to_one(donor, species)
                if not gene_map:
                    continue
                donor_frame = self._load_species(donor)
                donor_mask = donor_frame["cell_type"].isin(donor_names)
                if dl is not None:
                    donor_mask = donor_mask & donor_frame["_dl"].isin(dl)
                donor_sub = donor_frame.loc[donor_mask]
                if donor_sub.empty:
                    continue
                mapped = donor_sub["_gu"].map(gene_map)
                keep = mapped.notna()
                if not keep.any():
                    continue
                donor_sub = donor_sub[keep].copy()
                donor_sub["gene_symbol"] = mapped[keep]
                donor_sub["_gu"] = donor_sub["gene_symbol"].str.upper()
                parts.append(donor_sub[frame.columns])
        combined = pd.concat(parts, ignore_index=True) if len(parts) > 1 else parts[0]
        return ContextSources(combined, seed=seed)


class CompositeSources:

    def __init__(self, native, borrowed, borrowed_names):
        self.native = native
        self.borrowed = borrowed
        self.borrowed_names = {str(n) for n in borrowed_names}
        self.sub = native.sub

    def all_records(self, cell_type, gene):
        records = list(self.native.all_records(cell_type, gene))
        if str(cell_type) in self.borrowed_names:
            seen = {r["sentence"] for r in records}
            records += [r for r in self.borrowed.all_records(cell_type, gene) if r["sentence"] not in seen]
        return records

    def records_for_marker(self, cell_type, gene, k=MAX_PER_MARKER):
        return self.all_records(cell_type, gene)[:k]

    def for_marker(self, cell_type, gene, k=MAX_PER_MARKER):
        return [r["sentence"] for r in self.records_for_marker(cell_type, gene, k=k)]


class ContextSources:

    def __init__(self, sub, seed=SOURCE_ORDER_SEED):
        self.sub = sub
        self.seed = str(seed)
        self._by_cell_gene = None
        self._records = None

    def _index(self):
        if self._by_cell_gene is None:
            self._by_cell_gene = {
                key: [item["sentence"] for item in bucket]
                for key, bucket in self._record_index().items()
            }
        return self._by_cell_gene

    def _record_index(self):
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
        srcs = self._index().get((cell_type, str(gene).upper()), [])
        return srcs[:k]

    def all_records(self, cell_type, gene):
        bucket = self._record_index().get((cell_type, str(gene).upper()), [])
        return [
            {
                "sentence": item["sentence"],
                "pmid": item["pmid"],
                "pmcid": item["pmcid"],
                "n_pub": int(item["_n_pub"]),
            }
            for item in bucket
        ]

    def records_for_marker(self, cell_type, gene, k=MAX_PER_MARKER):
        return self.all_records(cell_type, gene)[:k]


class SourceServer:

    def __init__(self, context, batch=MAX_PER_MARKER, max_batches=0):
        self.context = context
        self.batch = int(batch)
        self.max_batches = int(max_batches)
        self._served: dict[tuple[str, str], set[str]] = {}
        self._draws: dict[tuple[str, str], int] = {}

    def _key(self, cell_type, gene):
        return (str(cell_type), str(gene).upper())

    def note_delivered(self, cell_type, gene, records):
        served = self._served.setdefault(self._key(cell_type, gene), set())
        for record in records:
            served.add(str(record.get("sentence", "")))

    def opening(self, cell_type, gene, k=None):
        records = self.context.records_for_marker(
            cell_type, gene, k=self.batch if k is None else int(k)
        )
        self.note_delivered(cell_type, gene, records)
        return records

    def take(self, cell_type, gene):
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
