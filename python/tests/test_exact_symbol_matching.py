"""Curated marker symbols meet uploaded gene symbols EXACTLY, in both arms.

WHY THIS EXISTS. The two sides of the match are written in different nomenclatures: HGNC
puts human symbols in upper case (PECAM1), MGI and RGD put mouse and rat symbols in title
case (Pecam1). Both crossings used to fold case before comparing, which meant a matrix
could match a menu curated for another species: on one mouse matrix (14,877 symbols,
93.3% title case) the Human/blood/Normal menu matched 1,029 genes under case folding and
4 under exact matching -- and those 4 are the real human/mouse homographs (C2, C3, F3,
F8). The run did not fail; it annotated confidently against 1,025 genes it had never
measured under those names.

Exact matching removes that path. It also means the engine analyses strictly the species
the caller chose and never second-guesses it: nothing here reads the object's symbols and
decides what species they "look like". A menu that matches nothing fails on the existing
zero-hit guard, whose message names both namespaces so the reader can see why.

A regression would not crash anything -- it would quietly restore confident wrong answers
-- so the two crossings and the guard's wording are pinned here rather than left to review.

WHAT IS NOT CASE-SENSITIVE, deliberately: `gene_key` inside scoring and reporting is still
an uppercased key. By the time it is built, both sides are already the matrix's own
symbols, and measured symbols are unique case-insensitively (preprocessing sums case
duplicates), so uppercasing is an injective relabelling within one namespace rather than a
second crossing. Those sites are left alone.

Offline; the R arm is read as source, since it has no importable entry point.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "scmarkeragent"
PREP_PY = PACKAGE / "preprocessing.py"
SCORING_PY = PACKAGE / "candidate_scoring.py"
PREP_R = PACKAGE / "rflow" / "preprocessing.R"
SCORING_R = PACKAGE / "rflow" / "candidate_scoring.R"
RSCRIPT = shutil.which("Rscript")


def _uncommented(path: Path) -> str:
    """Source with `#` comments stripped, so a claim in prose cannot satisfy a test."""
    return "\n".join(
        re.sub(r"#.*$", "", line)
        for line in path.read_text(encoding="utf-8").splitlines()
    )


class TestTheDeMenuCrossing(unittest.TestCase):
    """Crossing 1: the curated menu intersected with what the object measured."""

    @classmethod
    def setUpClass(cls):
        from scmarkeragent import preprocessing

        cls.menu_measured = staticmethod(preprocessing.menu_measured_genes)
        cls.sample = staticmethod(preprocessing._sample_symbols)

    def test_a_case_difference_is_not_a_match(self):
        """The regression: a mouse matrix against a human menu."""
        positive, negative = self.menu_measured(
            measured={"Pecam1", "Cd3e", "Actb"},
            menu_positive={"PECAM1", "CD3E"},
            menu_negative={"ACTB"},
        )
        self.assertEqual(positive, [])
        self.assertEqual(negative, [])

    def test_the_reverse_direction_is_also_refused(self):
        """A human matrix against a mouse menu."""
        positive, _ = self.menu_measured(
            measured={"PECAM1", "CD3E"},
            menu_positive={"Pecam1", "Cd3e"},
            menu_negative=set(),
        )
        self.assertEqual(positive, [])

    def test_an_exact_match_is_kept_in_the_matrix_namespace(self):
        positive, negative = self.menu_measured(
            measured={"Pecam1", "Cd3e", "Actb", "Gapdh"},
            menu_positive={"Pecam1", "Cd3e", "Absent1"},
            menu_negative={"Actb", "Absent2"},
        )
        self.assertEqual(positive, ["Cd3e", "Pecam1"])
        self.assertEqual(negative, ["Actb"])

    def test_a_gene_curated_both_ways_counts_as_positive_only(self):
        self.assertEqual(
            self.menu_measured({"Cd3e"}, {"Cd3e"}, {"Cd3e"}), (["Cd3e"], [])
        )

    def test_the_homographs_still_match(self):
        """C2, C3, F3 and F8 are spelled the same in HGNC and MGI. Exact matching keeps
        exactly these, which is the whole of the honest overlap."""
        positive, _ = self.menu_measured(
            measured={"C2", "C3", "F3", "F8", "Pecam1"},
            menu_positive={"C2", "C3", "F3", "F8", "PECAM1"},
            menu_negative=set(),
        )
        self.assertEqual(positive, ["C2", "C3", "F3", "F8"])

    def test_the_hgnc_mixed_case_symbols_still_match_themselves(self):
        """HGNC keeps 23 symbols with lower case in them (C1orf43, RNASEK-C17orf49).
        Exact matching must not cost a human run its own curated symbols."""
        positive, _ = self.menu_measured(
            measured={"C1orf43", "RNASEK-C17orf49"},
            menu_positive={"C1orf43", "RNASEK-C17orf49"},
            menu_negative=set(),
        )
        self.assertEqual(positive, ["C1orf43", "RNASEK-C17orf49"])

    def test_a_measured_gene_outside_the_menu_is_not_invented(self):
        positive, negative = self.menu_measured({"Pecam1", "Xist"}, {"Pecam1"}, set())
        self.assertEqual((positive, negative), (["Pecam1"], []))

    def test_the_sample_shown_in_a_failure_is_bounded_and_readable(self):
        self.assertEqual(self.sample({"B", "A", "D", "C", "E"}), "A, B, C, D")
        self.assertEqual(self.sample(set()), "nothing")


class TestTheEligibilityCrossing(unittest.TestCase):
    """Crossing 2: a curated candidate's panel against the measured genes."""

    @classmethod
    def setUpClass(cls):
        import pandas

        from scmarkeragent import marker_database

        cls.pd = pandas
        cls.eligible = staticmethod(marker_database.MarkerDatabase.eligible_candidates)

    def _panel(self, genes):
        return self.pd.DataFrame(
            [
                {
                    "cell": "endothelial cell",
                    "g": gene,
                    "marker_polarity": "positive",
                    "n_pub": 40,
                    "tier": "high",
                    "is_in_vitro": False,
                }
                for gene in genes
            ]
        )

    def test_a_candidate_does_not_qualify_on_case_folded_genes(self):
        _, eligible = self.eligible(
            self._panel(["Pecam1", "Cdh5", "Vwf"]), {"PECAM1", "CDH5", "VWF"}
        )
        self.assertEqual(eligible, [])

    def test_a_candidate_qualifies_on_exactly_matching_genes(self):
        positive, eligible = self.eligible(
            self._panel(["Pecam1", "Cdh5", "Vwf"]), {"Pecam1", "Cdh5", "Actb"}
        )
        self.assertEqual(eligible, ["endothelial cell"])
        self.assertEqual(
            sorted(positive["g"]),
            ["Cdh5", "Pecam1"],
            "only the genes this object measured may enter a panel",
        )


class TestBothArmsMatchTheSameWay(unittest.TestCase):
    """A caller picks an arm by input format alone, so the two must gate identically.
    Static over the R source, which has no importable entry point."""

    def test_the_r_de_menu_intersects_exactly(self):
        code = _uncommented(PREP_R)
        self.assertNotIn(
            "meas_by_key",
            code,
            "the case-folded menu-vs-measured index is back in the R arm",
        )
        self.assertNotIn("toupper(symbols)", code)
        self.assertIn("intersect(unique(as.character(symbols)), meas)", code)

    def test_the_python_de_menu_intersects_exactly(self):
        code = _uncommented(PREP_PY)
        self.assertNotIn(
            "meas_by_key", code, "the case-folded menu-vs-measured index is back"
        )
        self.assertIn("measured.intersection", code)

    def test_the_r_eligibility_gate_matches_exactly(self):
        code = _uncommented(SCORING_R)
        self.assertNotIn("toupper(g) %in% measured_upper", code)
        self.assertIn("g %in% measured_genes", code)

    def test_the_python_eligibility_gate_matches_exactly(self):
        code = _uncommented(PACKAGE / "marker_database.py")
        self.assertNotIn('panel["g"].str.upper().isin(meas)', code)
        self.assertIn('panel["g"].isin(meas)', code)

    def test_the_r_negative_panel_is_gated_like_the_positive_one(self):
        code = _uncommented(SCORING_R)
        negatives = re.search(r"negatives <- panel\[[^\]]*\]", code)
        self.assertIsNotNone(negatives, "the negative panel moved; re-point this test")
        self.assertIn(
            "g %in% measured_genes",
            negatives.group(0),
            "negatives cross the same boundary and must be gated the same way",
        )

    def test_the_python_negative_panel_is_gated_like_the_positive_one(self):
        code = _uncommented(SCORING_PY)
        self.assertIn('panel["g"].astype(str).isin(measured)', code)

    def test_neither_arm_reads_the_symbols_to_guess_a_species(self):
        """The engine analyses the species the caller chose. Nomenclature is not
        decidable well enough to overrule them (mouse and rat are written alike, and a
        legitimately re-symbolised object would be misread), so it is not attempted."""
        for path in (
            PREP_PY,
            PREP_R,
            SCORING_R,
            SCORING_PY,
            PACKAGE / "marker_database.py",
        ):
            with self.subTest(path=path.name):
                code = _uncommented(path)
                for token in (
                    "looks_like",
                    "species_guard",
                    "upperFraction",
                    "upper_fraction",
                    "DECISIVE",
                ):
                    self.assertNotIn(token, code)


class TestTheZeroHitGuardSaysWhatToDo(unittest.TestCase):
    """When nothing matches, the run stops on a message the reader can act on -- both
    namespaces named, and the two things that can be wrong. It reports; it does not
    diagnose the species."""

    # Both arms build the message from adjacent string literals, so it is read as the
    # caller will see it rather than as it is typed.
    PHRASES = (
        "this object measures",
        "the menu curates",
        "matched exactly, case included",
        "select the species these symbols belong to",
        "rewrite the symbols in the selected species' nomenclature",
        "tissue or disease does not describe this object",
    )

    def _message(self, path):
        source = re.sub(r'",?\s*\n\s*f?"', "", path.read_text(encoding="utf-8"))
        start = source.index("not one of the")
        return source[start : start + 900]

    def test_the_python_message_names_both_namespaces_and_the_fix(self):
        text = self._message(PREP_PY)
        for phrase in self.PHRASES:
            self.assertIn(phrase, text)

    def test_the_r_message_says_the_same_thing(self):
        text = self._message(PREP_R)
        for phrase in self.PHRASES:
            self.assertIn(phrase, text)

    def test_the_message_names_no_server_path(self):
        """A caller reads this message; it must not carry a filesystem path."""
        for path in (PREP_PY, PREP_R):
            with self.subTest(path=path.name):
                self.assertIsNone(
                    re.search(r"(?:/[A-Za-z0-9_.+-]+){2,}", self._message(path))
                )

    def test_the_python_guard_actually_fires_on_an_empty_intersection(self):
        from scmarkeragent import preprocessing

        positive, _ = preprocessing.menu_measured_genes({"Pecam1"}, {"PECAM1"}, set())
        self.assertEqual(positive, [], "the precondition the guard reacts to")
        self.assertIn("if not positive_genes:", PREP_PY.read_text(encoding="utf-8"))


class TestTheOrthologTablesSpeakEachSpeciesNomenclature(unittest.TestCase):
    """An ortholog-mapped symbol crosses the same boundary, one step earlier.

    Cross-species pooling rewrites another species' curated symbol into the query
    species' symbol space, and that rewritten symbol is then matched to measured genes
    exactly. While both ortholog columns were stored upper-cased the rewrite produced
    `PECAM1` for a mouse target, which no MGI matrix carries, so the switch was inert for
    mouse and rat. The tables were recast on 2026-08-12; upper-casing either column again
    would silently switch the feature back off while the run still reports the pooling.
    """

    VENDORED = PACKAGE / "resources" / "static" / "ortholog"

    @classmethod
    def setUpClass(cls):
        import pandas  # noqa: F401

        from scmarkeragent import config, ortho_map

        cls.ortho_map = ortho_map
        cls.ortho = Path(config.ORTHO_DIR)
        if not (cls.ortho / "ortho_Human_to_Mouse.csv").is_file():
            raise unittest.SkipTest("ortholog tables not present on this host")

    def _map(self):
        return self.ortho_map.OrthoMap(str(self.ortho))

    def test_the_stored_tables_carry_the_target_species_own_spelling(self):
        for target in ("Mouse", "Rat"):
            with self.subTest(target=target):
                rows = (
                    (self.ortho / f"ortho_Human_to_{target}.csv")
                    .read_text(encoding="utf-8")
                    .splitlines()[1:]
                )
                lowered = sum(
                    1 for row in rows if row.split(",")[1] != row.split(",")[1].upper()
                )
                self.assertGreater(
                    lowered,
                    len(rows) * 0.9,
                    f"{target} targets look upper-cased again",
                )

    def test_a_rodent_target_is_handed_out_in_title_case(self):
        pairs = self._map().pairs("Human", "Mouse")
        row = pairs[pairs["src_sym"] == "PECAM1"]
        self.assertEqual(list(row["tgt_sym"]), ["Pecam1"])

    def test_a_symbol_a_capitalisation_rule_would_get_wrong(self):
        """`HLA-DQA1` is `H2-Aa`, not `H2-aa`. Recovery was a lookup, not a rule."""
        pairs = self._map().pairs("Human", "Mouse")
        self.assertEqual(
            list(pairs[pairs["src_sym"] == "HLA-DQA1"]["tgt_sym"]), ["H2-Aa"]
        )
        self.assertIn("S100a8", set(pairs["tgt_sym"]))
        self.assertIn("C1qa", set(pairs["tgt_sym"]))

    def test_the_source_side_stays_an_uppercase_key(self):
        """`_pool_cross_species` upper-cases the curated symbol before merging on it, so
        a native source column would match nothing."""
        for src, tgt in (("Human", "Mouse"), ("Mouse", "Human"), ("Rat", "Mouse")):
            with self.subTest(pair=f"{src}->{tgt}"):
                pairs = self._map().pairs(src, tgt)
                self.assertTrue(
                    (pairs["src_sym"] == pairs["src_sym"].str.upper()).all()
                )

    def test_a_human_target_keeps_the_hgnc_spelling_including_its_lower_case(self):
        pairs = self._map().pairs("Mouse", "Human")
        self.assertIn("C1orf43", set(pairs["tgt_sym"]))

    def test_composing_mouse_to_rat_does_not_lose_the_mixed_case_hub_symbols(self):
        """The composition crosses the hub on the folded key. Crossing on the spelling
        would drop every pair routed through an HGNC symbol carrying lower case."""
        composed = self._map().pairs("Mouse", "Rat")
        self.assertGreater(len(composed), 16000)
        self.assertTrue((composed["src_sym"] == composed["src_sym"].str.upper()).all())
        self.assertIn("S100a8", set(composed["tgt_sym"]))

    def test_the_vendored_fallback_copy_is_the_same_bytes(self):
        """Both roots carry a copy and `resource_path` falls back between them, so a
        bundle that failed to mount must not answer with the old casing."""
        if self.ortho.resolve() == self.VENDORED.resolve():
            self.skipTest("the vendored copy IS the resolved one here")
        for name in ("ortho_Human_to_Mouse.csv", "ortho_Human_to_Rat.csv"):
            with self.subTest(name=name):
                self.assertEqual(
                    (self.ortho / name).read_bytes(),
                    (self.VENDORED / name).read_bytes(),
                )

    def test_neither_arm_folds_a_target_column_on_load(self):
        for path, banned in (
            (
                PACKAGE / "ortho_map.py",
                ('d["tgt"] = d["tgt"].str.upper()', 'd["hum"] = d["hum"].str.upper()'),
            ),
            (
                PACKAGE / "rflow" / "ortho_map.R",
                ("hum := toupper(hum)", "tgt := toupper(tgt)"),
            ),
        ):
            with self.subTest(path=path.name):
                code = _uncommented(path)
                for line in banned:
                    self.assertNotIn(line, code)

    def test_every_direction_returns_the_pair_count_the_recast_preserved(self):
        """Only the spelling changed: no row was added, dropped or re-paired."""
        expected = {
            ("Human", "Mouse"): 16819,
            ("Mouse", "Human"): 16819,
            ("Human", "Rat"): 16752,
            ("Rat", "Human"): 16752,
            ("Mouse", "Rat"): 16369,
            ("Rat", "Mouse"): 16369,
        }
        mapper = self._map()
        for (src, tgt), count in expected.items():
            with self.subTest(pair=f"{src}->{tgt}"):
                self.assertEqual(len(mapper.pairs(src, tgt)), count)


class TestTheEngineStillParses(unittest.TestCase):
    def test_the_r_arm_parses(self):
        if not RSCRIPT:
            self.skipTest("Rscript not available on this host")
        for path in (PREP_R, SCORING_R):
            with self.subTest(path=path.name):
                result = subprocess.run(
                    [RSCRIPT, "-e", f"invisible(parse({str(path)!r}))"],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                self.assertEqual(
                    result.returncode,
                    0,
                    f"{path.name} does not parse: {result.stderr}",
                )


if __name__ == "__main__":
    unittest.main()
