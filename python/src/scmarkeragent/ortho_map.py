#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

from .config import ORTHO_DIR

SPECIES = ("Human", "Mouse", "Rat")


class OrthoMap:
    def __init__(self, cache=ORTHO_DIR):
        self.cache = cache
        self._hub = {}
        self._pairs = {}
        self._one_to_one = {}

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
        d["hum_key"] = d["hum"].str.upper()
        d["tgt_key"] = d["tgt"].str.upper()
        d = d.drop_duplicates().reset_index(drop=True)
        self._hub[other] = d
        return d

    def pairs(self, src, tgt):
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
            a = self.pairs(src, "Human")
            b = self.pairs("Human", tgt)
            a = a.assign(hub_key=a["tgt_sym"].str.upper())
            d = a.merge(
                b, left_on="hub_key", right_on="src_sym", suffixes=("_a", "_b")
            ).rename(columns={"src_sym_a": "src_sym", "tgt_sym_b": "tgt_sym"})
        d = d[["src_sym", "tgt_sym"]].dropna().drop_duplicates().reset_index(drop=True)
        self._pairs[(src, tgt)] = d
        return d

    def one_to_one(self, src, tgt) -> dict[str, str]:
        key = (src, tgt)
        if key in self._one_to_one:
            return self._one_to_one[key]
        if src == tgt:
            self._one_to_one[key] = {}
            return self._one_to_one[key]
        d = self.pairs(src, tgt).copy()
        d["tgt_key"] = d["tgt_sym"].str.upper()
        fwd_ok = d.groupby("src_sym")["tgt_key"].nunique() == 1
        rev_ok = d.groupby("tgt_key")["src_sym"].nunique() == 1
        d = d[d["src_sym"].isin(fwd_ok[fwd_ok].index) & d["tgt_key"].isin(rev_ok[rev_ok].index)]
        out = dict(zip(d["src_sym"], d["tgt_sym"]))
        self._one_to_one[key] = out
        return out


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
