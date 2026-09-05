#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""UBERON anatomy-ontology helpers: tissue subsumption for the tissue_root switch.

Mirrors cl_ontology.py in style and normalization, but:
  - parses uberon-basic.obo (id: UBERON:..., is_a: UBERON:..., and
    relationship: part_of UBERON:...) so a *parent* edge = is_a OR part_of.
    part_of is what makes "colon" reachable from "intestine"; is_a covers
    anatomical subtype edges. Non-UBERON targets (CARO/BFO) are ignored.
  - adds descendants(id) / members(root) so a query tissue can pool every
    same-organ anatomical sub-region under a root (e.g. intestine -> colon,
    ileum, Peyer's patch, vermiform appendix ...).

Deterministic, file-only: no network, no LLM. Both arms parse the obo at
runtime, so the ontology file is the single source of truth (no committed
derived artifact to drift).
"""

import re
from .config import OBO_UBERON as OBO_DEFAULT

_DASH = re.compile(r"[\u2010-\u2015\u2212]")
_PAREN_TAIL = re.compile(r"\s*\([^)]*\)\s*$")
_WS = re.compile(r"\s+")
_SP_HYPHEN = re.compile(r" ?- ?")
_SYN = re.compile(r'^synonym: "([^"]*)" *([A-Z]*).*$')
_ISA = re.compile(r"^is_a:\s*(UBERON:\d+)")
_PARTOF = re.compile(r"^relationship:\s*part_of\s+(UBERON:\d+)")


def norm_name(x):
    if x is None:
        return ""
    x = str(x).strip().lower()
    x = _PAREN_TAIL.sub("", x)
    x = _DASH.sub("-", x)
    x = _WS.sub(" ", x)
    x = _SP_HYPHEN.sub("-", x)
    return x.strip()


class UberonOntology:
    def __init__(self, obo=OBO_DEFAULT):
        self.id2name = {}
        self.parents = {}  # id -> [is_a + part_of parents]
        self._children = {}  # id -> [children] (inverse of parents)
        self._anc = {}  # ancestor cache
        self._desc = {}  # descendant-set cache
        self.namemap = {}  # normalized name/EXACT-synonym -> id
        self._parse(obo)
        self._invert_children()
        self._build_namemap_from_parsed()

    # ---- parse ----------------------------------------------------------------
    def _parse(self, obo):
        prim: dict[str, set[str]] = {}
        syn: dict[str, set[str]] = {}
        cid = None
        name = None
        obs = False
        syns: list[str] = []

        def flush():
            if cid is None or obs:
                return
            if name is not None:
                prim.setdefault(norm_name(name), set()).add(cid)
            for s in syns:
                syn.setdefault(norm_name(s), set()).add(cid)

        with open(obo, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                x = raw.rstrip("\n")
                if x == "[Term]":
                    flush()
                    cid = None
                    name = None
                    obs = False
                    syns = []
                    continue
                if x.startswith("id: UBERON:"):
                    cid = x[4:].strip()
                elif x.startswith("name: "):
                    name = x[6:]
                    if cid is not None:
                        self.id2name[cid] = name
                elif x.startswith("is_obsolete: true"):
                    obs = True
                elif x.startswith("is_a:"):
                    m = _ISA.match(x)
                    if m and cid is not None:
                        self.parents.setdefault(cid, []).append(m.group(1))
                elif x.startswith("relationship: part_of"):
                    m = _PARTOF.match(x)
                    if m and cid is not None:
                        self.parents.setdefault(cid, []).append(m.group(1))
                elif x.startswith("synonym: "):
                    m = _SYN.match(x)
                    if m and m.group(1) and m.group(2) == "EXACT":
                        syns.append(m.group(1))
        flush()
        self._prim = prim
        self._syn = syn

    def _build_namemap_from_parsed(self):
        prim_v = {k: next(iter(v)) for k, v in self._prim.items() if len(v) == 1}
        syn_v = {
            k: next(iter(v))
            for k, v in self._syn.items()
            if len(v) == 1 and k not in prim_v
        }
        nm = dict(prim_v)
        nm.update(syn_v)
        self.namemap = nm

    def _invert_children(self):
        for k, ps in self.parents.items():
            for p in set(ps):
                self._children.setdefault(p, []).append(k)

    # ---- graph ops ------------------------------------------------------------
    def ancestors(self, uid):
        if not uid:
            return set()
        if uid in self._anc:
            return self._anc[uid]
        out = set()
        frontier = list(self.parents.get(uid, []))
        while frontier:
            nxt = []
            for f in frontier:
                if f not in out:
                    out.add(f)
                    nxt.extend(self.parents.get(f, []))
            frontier = nxt
        self._anc[uid] = out
        return out

    def descendants(self, uid):
        """All ids reachable downward via is_a/part_of children (excludes uid)."""
        if not uid:
            return set()
        if uid in self._desc:
            return self._desc[uid]
        seen = set()
        frontier = list(self._children.get(uid, []))
        while frontier:
            nxt = []
            for f in frontier:
                if f not in seen:
                    seen.add(f)
                    nxt.extend(self._children.get(f, []))
            frontier = nxt
        self._desc[uid] = seen
        return seen

    def is_desc(self, a, b):
        """a is a descendant of b (b is ancestor of a)."""
        return bool(a) and bool(b) and (b in self.ancestors(a))

    def closure_self(self, uid):
        """{uid} + all descendants -- the set to pool for a root tissue."""
        if not uid:
            return set()
        return {uid} | self.descendants(uid)

    # ---- name <-> id ----------------------------------------------------------
    def name_of(self, uid):
        return self.id2name.get(uid, uid)

    def resolve(self, name):
        """Normalized name / EXACT synonym -> UBERON id (None if unresolved).
        Accepts a raw 'UBERON:xxxxxxx' id directly."""
        if not name:
            return None
        s = str(name).strip()
        if s.upper().startswith("UBERON:"):
            return s.upper()
        return self.namemap.get(norm_name(s))

    def members(self, root):
        """Resolve a root tissue (name or id) and return {root_id} + descendant
        ids. Empty set if the root name cannot be resolved."""
        rid = self.resolve(root)
        if not rid:
            return set()
        return self.closure_self(rid)


if __name__ == "__main__":
    import time

    t0 = time.time()
    ub = UberonOntology()
    print(f"parsed uberon-basic.obo in {time.time() - t0:.1f}s")
    print(f"  terms (id2name): {len(ub.id2name)}")
    print(f"  parents entries: {len(ub.parents)}")
    print(f"  namemap entries: {len(ub.namemap)}")
    for nm in ["intestine", "brain", "kidney", "heart", "lung", "pancreas", "stomach"]:
        rid = ub.resolve(nm)
        mem = ub.members(nm)
        kids = sorted(ub.name_of(d) for d in ub.descendants(rid)) if rid else []
        print(
            f"  {nm!r:12s} -> {rid} ({ub.name_of(rid) if rid else 'NA'}); "
            f"members={len(mem)}; e.g. {kids[:6]}"
        )
    # subsumption sanity: intestine must NOT pull in stomach/esophagus
    it = ub.resolve("intestine")
    for probe in [
        "stomach",
        "esophagus",
        "colon",
        "ileum",
        "Peyer's patch",
        "vermiform appendix",
    ]:
        pid = ub.resolve(probe)
        print(f"  is_desc({probe!r:20s}, intestine) = {ub.is_desc(pid, it)}")
