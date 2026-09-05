#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ortho_map.py -- bidirectional ortholog SYMBOL mapping across Human/Mouse/Rat
for the optional cross-species marker-sharing switch.

Built OFFLINE from the two human-hub tables exported by ortholog.R:
  cache/ortho_Human_to_Mouse.csv, cache/ortho_Human_to_Rat.csv  (cols hum,tgt, each
  in ITS OWN SPECIES' NATIVE CASING: A1BG,A1bg)

Direction support (all 6 ordered pairs):
  Human<->Mouse, Human<->Rat : direct table (reversed when query is the hub).
  Mouse<->Rat                : COMPOSED via Human (Mouse->Human->Rat). Composition
                               is the standard hub approach; it is slightly more
                               permissive on 1:many but only used when neither side
                               is Human.
1:many orthologs are ALL kept (a source symbol can expand to several target
symbols), matching ortholog.R. ASCII-only source.

WHY THE TWO SIDES ARE CASED DIFFERENTLY. `tgt_sym` is written into the query species'
symbol space and is then matched to that object's measured genes EXACTLY, so it has to
carry the spelling that species is written in: an all-upper PECAM1 handed to a mouse
matrix matches nothing. `src_sym` is only ever the left side of a merge against curated
symbols that the caller uppercases first, so it stays an uppercase key -- which is also
what lets one table serve a direction and its reverse.
"""

import os
import pandas as pd

from .config import ORTHO_DIR

SPECIES = ("Human", "Mouse", "Rat")


class OrthoMap:
    def __init__(self, cache=ORTHO_DIR):
        self.cache = cache
        self._hub = {}  # other -> DataFrame[hum, tgt]
        self._pairs = {}  # (src,tgt) -> DataFrame[src_sym, tgt_sym]

    def _load_hub(self, other):
        if other in self._hub:
            return self._hub[other]
        f = os.path.join(self.cache, f"ortho_Human_to_{other}.csv")
        if not os.path.exists(f):
            raise FileNotFoundError(
                f"ortholog table missing: {f} (run the ortho export; see ortholog.R)"
            )
        d = pd.read_csv(f, dtype=str)
        d = d.iloc[:, :2]
        d.columns = ["hum", "tgt"]
        d = d.dropna()
        # Both spellings are kept as the table stores them, and an uppercase key is
        # derived beside each. Either column can be the TARGET of a direction, so
        # neither may be folded in place.
        d["hum_key"] = d["hum"].str.upper()
        d["tgt_key"] = d["tgt"].str.upper()
        d = d.drop_duplicates().reset_index(drop=True)
        self._hub[other] = d
        return d

    def pairs(self, src, tgt):
        """DataFrame[src_sym, tgt_sym] (deduped) mapping src species symbols to tgt
        species symbols. `src_sym` is an UPPERCASE key, `tgt_sym` the target species'
        NATIVE spelling. Returns None when src == tgt (identity)."""
        if src == tgt:
            return None
        if (src, tgt) in self._pairs:
            return self._pairs[(src, tgt)]
        if src not in SPECIES or tgt not in SPECIES:
            raise ValueError(f"unsupported species pair: {src}->{tgt}")
        if src == "Human" and tgt in ("Mouse", "Rat"):
            d = self._load_hub(tgt).rename(
                columns={"hum_key": "src_sym", "tgt": "tgt_sym"}
            )
        elif tgt == "Human" and src in ("Mouse", "Rat"):
            d = self._load_hub(src).rename(
                columns={"tgt_key": "src_sym", "hum": "tgt_sym"}
            )
        else:
            # Mouse <-> Rat: compose via Human. The hub is crossed on the uppercase key,
            # because the first leg now delivers Human in its native spelling and the
            # second leg keys on the folded one -- the 23 HGNC symbols carrying lower
            # case (C1orf43) would otherwise drop out of the composition.
            a = self.pairs(src, "Human")  # src -> Human
            b = self.pairs("Human", tgt)  # Human -> tgt
            a = a.assign(hub_key=a["tgt_sym"].str.upper())
            d = a.merge(
                b, left_on="hub_key", right_on="src_sym", suffixes=("_a", "_b")
            ).rename(columns={"src_sym_a": "src_sym", "tgt_sym_b": "tgt_sym"})
        d = d[["src_sym", "tgt_sym"]].dropna().drop_duplicates().reset_index(drop=True)
        self._pairs[(src, tgt)] = d
        return d


if __name__ == "__main__":
    om = OrthoMap()
    for a, b in [
        ("Human", "Mouse"),
        ("Mouse", "Human"),
        ("Human", "Rat"),
        ("Rat", "Human"),
        ("Mouse", "Rat"),
        ("Rat", "Mouse"),
    ]:
        pr = om.pairs(a, b)
        print(
            f"{a:5s} -> {b:5s} : {len(pr):6d} ortholog symbol pairs; "
            f"e.g. {list(zip(pr['src_sym'][:3], pr['tgt_sym'][:3]))}"
        )
