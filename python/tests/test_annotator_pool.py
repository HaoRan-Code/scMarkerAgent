"""The evidence the annotating agent opens with, and the tools it can call for more."""

from __future__ import annotations

from scmarkeragent import annotator_pool as pool_api


def _markers(pool, cluster, candidate):
    return pool_api.find_candidate(pool, cluster, candidate)["markers"]


def test_a_panel_arrives_whole_and_best_corroborated_first(pool):
    """The ordering is the design: the gene the resource says defines the cell is first.

    Publication support is the curated database's own record of which gene identifies
    the identity. Sorting the panel on it means an agent reading from the top meets that
    gene before any incidental marker, whatever the measurement next to it turns out to
    be.
    """
    rows = _markers(pool, "1", "enterocyte")
    assert [row["gene"] for row in rows] == ["VIL1", "SI", "FABP2", "RBP2", "PTPRC"]
    assert [row["n_pub"] for row in rows] == [53, 34, 24, 3, 1]
    assert rows[0]["polarity"] == "positive"
    assert rows[-1]["polarity"] == "negative"
    # The measured contrast, the curated identity of the gene, and the two verdicts
    # retrieval already reached on this row: `sig_pass` admitted it and `m_g` weighed it.
    # `avg_log2FC` still says whether the gene is up; `significant` says whether that
    # rise clears the gate, which is a different question and a different answer on 18%
    # of rows.
    assert set(rows[0]) == {
        "gene",
        "polarity",
        "n_pub",
        "tier",
        "pct_in",
        "pct_out",
        "avg_log2FC",
        "significant",
        "specificity",
    }


def test_detection_is_carried_as_percent_on_the_de_tables_own_scale(pool):
    row = _markers(pool, "1", "enterocyte")[2]
    assert row["gene"] == "FABP2"
    assert row["pct_in"] == 92.4
    assert row["pct_out"] == 68.5



def test_an_unraised_marker_keeps_the_two_readings_distinguishable(pool):
    """Filtering the panel to raised markers would erase both of these cases.

    A dominant population's defining markers go unraised because the comparison is
    against the rest of the same tissue -- FABP2 at 92.4 in versus 68.5 out is the
    enterocyte evidence, and it has a negative fold change. An identity the cells do not
    carry also goes unraised, but near zero: CHGA at 2.3 against 4.1. Showing only what
    rose hides the first as evidence and the second as refutation.
    """
    enterocyte = {row["gene"]: row for row in _markers(pool, "1", "enterocyte")}
    endocrine = {row["gene"]: row for row in _markers(pool, "1", "enteroendocrine cell")}

    assert not pool_api.is_raised(enterocyte["FABP2"])
    assert enterocyte["FABP2"]["pct_in"] > 90.0

    assert not pool_api.is_raised(endocrine["CHGA"])
    assert endocrine["CHGA"]["pct_in"] < 5.0
    assert endocrine["CHGA"]["n_pub"] == 160


def test_raised_is_exactly_a_positive_fold_change(pool):
    rows = {row["gene"]: row for row in _markers(pool, "1", "enterocyte")}
    assert pool_api.is_raised(rows["SI"])  # 0.01 is a raise
    assert not pool_api.is_raised(rows["VIL1"])  # -0.22 is not
    assert pool_api.is_raised({"avg_log2FC": 0.0001})
    assert not pool_api.is_raised({"avg_log2FC": 0.0})
    assert not pool_api.is_raised({})


def test_a_curated_gene_the_dataset_never_measured_is_counted_not_shown(pool):
    """It is neither evidence for the identity nor against it, so it carries no row."""
    entry = pool_api.find_candidate(pool, "1", "enteroendocrine cell")
    assert entry["unmeasured_curated_genes"]["count"] == 1
    assert entry["unmeasured_curated_genes"]["genes"] == ["MISSING1"]
    assert "MISSING1" not in {row["gene"] for row in entry["markers"]}


def test_the_opening_packet_carries_sources_only_behind_a_detected_exclusion(pool):
    """Sources are the largest evidence available and the smallest part of it that decides.

    They are fetched by name, which is the whole reason the agent has tools. The one
    exception is a negative marker the cluster detects: its weight turns entirely on what
    the literature establishes, so those sentences arrive with the panel rather than
    waiting for an agent to think of asking.
    """
    packet = pool_api.cluster_packet(pool, "1")
    assert packet["query"]["cluster_id"] == "1"
    assert packet["query"]["cells_in_cluster"] == 288
    assert [entry["cell_type"] for entry in packet["candidates"]] == [
        "enterocyte",
        "enteroendocrine cell",
        "tuft cell",
    ]
    assert set(packet["candidates"][0]) == {
        "cell_type",
        "markers",
        "exclusion_sources",
        "unmeasured_curated_genes",
        "program",
    }
    # The fixture's only negative marker sits at pct_in 0.0, far under the threshold, so
    # no sentence reaches this packet and nothing else grew a source field.
    assert all(entry["exclusion_sources"] == {} for entry in packet["candidates"])
    assert "sentence" not in repr(packet)


def test_a_detected_exclusion_ships_its_sentences(pool_with_detected_exclusion):
    """The provenance follows the measurement, not the candidate or the polarity alone."""
    pool, sources = pool_with_detected_exclusion
    packet = pool_api.cluster_packet(pool, "1")
    by_name = {entry["cell_type"]: entry for entry in packet["candidates"]}

    shipped = by_name["enterocyte"]["exclusion_sources"]
    assert set(shipped) == {"PTPRC"}
    assert shipped["PTPRC"][0]["sentence"] == "PTPRC is absent from enterocyte."
    assert shipped["PTPRC"][0]["pmid"] == "1"
    # Only the detected exclusion was looked up; no positive marker pulled a sentence.
    assert sources.asked == [("enterocyte", "PTPRC")]

    # Cluster 2 leaves PTPRC undetected, so the same candidate ships nothing there.
    assert pool_api.cluster_packet(pool, "2")["candidates"][0]["exclusion_sources"] == {}
    # What the agent is shown of a row, which is what keeps the two arms aligned. `tier`
    # is not in it: the resource derives it from `n_pub` by a fixed rule, so it is the
    # same fact twice, and the packet's field list never described it.
    assert set(by_name["enterocyte"]["markers"][0]) == {
        "gene", "polarity", "n_pub", "specificity",
        "pct_in", "pct_out", "avg_log2FC", "significant",
    }


def test_the_cluster_table_names_every_claimant_of_a_shared_gene(pool):
    """The fact a panel cannot carry: how many of these candidates claim the same gene.

    PROX1 is curated for two of the three candidates, and each of their panels shows only
    its own claim. Read from the panels, PROX1 looks like a marker of whichever identity
    is being read; read here, it is a gene the two share, and the widest contrast in the
    cluster belongs to a gene that separates nothing.
    """
    table = pool_api.cluster_marker_table(pool, "1")

    # Ordered by contrast, widest first -- not by publication support, which would open
    # on SI (34 publications) and put RBP2 (3) last.
    assert [row["gene"] for row in table] == ["PROX1", "RBP2", "SI"]
    assert set(table[0]) == {"gene", "pct_in", "pct_out", "specificity", "claimed_by"}
    assert table[0]["claimed_by"] == [
        {"cell_type": "enteroendocrine cell", "n_pub": 4},
        {"cell_type": "tuft cell", "n_pub": 1},
    ]
    assert table[1]["claimed_by"] == [{"cell_type": "enterocyte", "n_pub": 3}]
    # Unraised and near-zero rows are absent: CHGA is the enteroendocrine panel's whole
    # refutation and belongs on that panel, but it separates nothing and is not evidence
    # this cluster is anything.
    assert {"FABP2", "VIL1", "CHGA", "NEUROD1", "POU2F3"}.isdisjoint(
        row["gene"] for row in table
    )
    # Negative markers have their own channel; this table is the positive comparison.
    assert "PTPRC" not in {row["gene"] for row in table}
    assert pool_api.cluster_packet(pool, "1")["cluster_markers"] == table


def test_gene_across_clusters_exposes_a_program_the_packet_cannot_show(pool):
    """The opening packet holds one cluster, so a shared program is invisible in it."""
    observation = pool_api.run_tool(pool, "1", "gene_across_clusters", {"gene": "fabp2"})
    assert observation["gene"] == "FABP2"
    assert observation["clusters"] == [
        {"cluster_id": "1", "pct_in": 92.4, "pct_out": 68.5},
        {"cluster_id": "2", "pct_in": 10.0, "pct_out": 80.0},
    ]
    missing = pool_api.run_tool(pool, "1", "gene_across_clusters", {"gene": "NOPE"})
    assert missing["error"] == "not measured in this dataset"


def test_candidates_with_gene_answers_only_for_this_clusters_pool(pool):
    observation = pool_api.run_tool(pool, "1", "candidates_with_gene", {"gene": "PROX1"})
    assert [row["cell_type"] for row in observation["candidates"]] == [
        "enteroendocrine cell",
        "tuft cell",
    ]
    assert observation["candidates"][0]["n_pub"] == 4


def test_sources_are_capped_per_call_and_bound_to_a_candidate_of_this_cluster(pool):
    class _Sources:
        def records_for_marker(self, cell_type, gene, k):
            return [
                {"pmid": "1", "pmcid": "PMC1", "sentence": f"{gene} marks {cell_type}."}
            ][:k]

    observation = pool_api.run_tool(
        pool, "1", "sources", {"candidate": "enterocyte", "genes": ["fabp2"]}, _Sources()
    )
    assert observation["sources"]["FABP2"][0]["pmcid"] == "PMC1"
    assert observation["truncated"] is False

    many = pool_api.run_tool(
        pool,
        "1",
        "sources",
        {
            "candidate": "enterocyte",
            "genes": [f"G{i}" for i in range(pool_api.SOURCE_GENES_PER_QUERY + 4)],
        },
        _Sources(),
    )
    assert len(many["sources"]) == pool_api.SOURCE_GENES_PER_QUERY
    assert many["truncated"] is True
    # The budget is what the template promises, and 64% of the calls in the shipped run
    # asked for exactly it -- so a cluster wanting more paid a second turn, which re-sends
    # the whole 26k-token packet to buy a few sentences.
    assert pool_api.SOURCE_GENES_PER_QUERY == 12

    unknown = pool_api.run_tool(
        pool, "1", "sources", {"candidate": "plasma cell", "genes": ["ALB"]}, _Sources()
    )
    assert "not a candidate" in unknown["error"]


def test_an_unusable_call_is_answered_rather_than_raised(pool):
    """The agent is mid-conversation; telling it what went wrong lets it correct itself."""
    observation = pool_api.run_tool(pool, "1", "gene_stats", {"genes": ["ALB"]})
    assert "unknown tool" in observation["error"]
    assert pool_api.TOOLS == (
        "sources",
        "gene_across_clusters",
        "candidates_with_gene",
        "pool_search",
    )


def test_the_agent_can_ask_what_the_retrieved_menu_left_out(pool):
    """Retrieval's cut has to be visible to something, or a miss cannot be reported.

    `candidates_with_gene` answers only about the candidates that were supplied, so before
    this tool nothing downstream could tell "no other identity claims this gene" from "the
    identity that claims it was not retrieved". Those are opposite facts.
    """
    inside = pool_api.run_tool(pool, "1", "candidates_with_gene", {"gene": "PROX1"})
    outside = pool_api.run_tool(pool, "1", "pool_search", {"gene": "PROX1"})
    here = set(pool_api.candidate_names(pool, "1"))

    assert {row["cell_type"] for row in inside["candidates"]} <= here
    named = {row["cell_type"] for row in outside["not_in_candidates"]}
    assert named == {"M cell"} and not (named & here)
    # And it is the better-published claimant, which is the case that matters: the two
    # candidates claiming PROX1 hold 4 and 1 publications for it against M cell's 8.
    assert outside["not_in_candidates"][0]["n_pub"] > max(
        row["n_pub"] for row in inside["candidates"]
    )
    # A name and a publication count, deliberately not a panel: evidence that the menu may
    # be short, not a sixteenth candidate to answer with.
    assert set(outside["not_in_candidates"][0]) == {
        "cell_type",
        "n_pub",
        "tier",
        "polarity",
    }
    assert len(outside["not_in_candidates"]) <= pool_api.POOL_SEARCH_ROWS
    assert outside["total"] >= len(outside["not_in_candidates"])


def test_a_warning_reads_the_whole_panel_and_not_the_quoted_gene(pool):
    """The check has to be able to see what the claim did not mention.

    A claim can be written on one incidental marker with a wide contrast while the genes
    the resource says define that cell sit flat. Auditing only the quoted gene would
    find nothing wrong with exactly the claims that are wrong, so the warning is
    computed over the identity's whole measured positive panel.
    """
    warned = pool_api.claim_warnings(pool, "1", [("cooc", "enteroendocrine cell")])
    assert len(warned) == 1
    name, line = warned[0]
    assert name == "enteroendocrine cell"
    assert line.startswith("cooc enteroendocrine cell: not_raised 2/3")
    # Ordered by publication support, so the most telling absence is named first.
    assert "CHGA(2.3/4.1)" in line
    assert line.index("CHGA") < line.index("NEUROD1")
    assert line.endswith("+0 more")


def test_a_claim_whose_panel_is_fully_raised_carries_no_warning(pool):
    assert pool_api.claim_warnings(pool, "2", [("selected", "tuft cell")]) == []


def test_a_finer_name_is_exclusive_only_where_the_parent_never_curates_the_gene(pool):
    """The comparison is against the parent's WHOLE panel, not the part it raises here.

    PROX1 is raised in cluster 2 and is curated for both identities, so it is the parent
    program read at a narrower name however wide its contrast; POU2F3 is curated for the
    finer identity alone and is the only marker that separates them.
    """
    assert pool_api.subtype_exclusive_markers(
        pool, "2", "enteroendocrine cell", "tuft cell"
    ) == ["POU2F3"]


def test_an_unraised_exclusive_marker_does_not_count(pool):
    """The gene has to be exclusive AND up here; being curated apart establishes nothing."""
    assert pool_api.subtype_exclusive_markers(
        pool, "1", "enteroendocrine cell", "tuft cell"
    ) == []
    assert pool_api.subtype_exclusive_markers(pool, "1", "enterocyte", "plasma cell") == []


def test_the_exclusive_requirement_is_more_than_one_marker(pool):
    """One exclusive gene can be incidental; the requirement is what makes it a refinement."""
    assert pool_api.SUBTYPE_EXCLUSIVE_REQUIRED == 2
    assert len(
        pool_api.subtype_exclusive_markers(pool, "2", "enterocyte", "enteroendocrine cell")
    ) == 3


def test_a_dominant_population_is_flagged_and_not_corrected(pool):
    """The warning is a flag, so a correct call on a dominant population still gets one.

    Reporting it is the honest outcome: three of the enterocyte panel's four positive
    markers really are unraised here, and the reason -- a high `pct_out` beside a high
    `pct_in` -- is in the line itself for a reader to weigh.
    """
    name, line = pool_api.claim_warnings(pool, "1", [("selected", "enterocyte")])[0]
    assert name == "enterocyte"
    assert line.startswith("selected enterocyte: not_raised 2/4")
    assert "VIL1(70.8/72.1)" in line
    assert "FABP2(92.4/68.5)" in line


def test_a_warning_names_the_identity_even_when_the_name_carries_a_colon():
    """212 curated cell-type names hold a colon of their own.

    The identity used to be recovered downstream by cutting the warning line at its
    first colon, which filed this one under `BCR` -- an identity that does not exist.
    The name now travels with its line, so nothing is parsed back out of the string.
    """
    name = "BCR::ABL1 positive primitive cell"
    colon_pool = {
        "clusters": {
            "1": {
                "candidates": [
                    {
                        "cell_type": name,
                        "markers": [
                            {
                                "gene": "ABL1",
                                "polarity": "positive",
                                "avg_log2FC": -1.0,
                                "pct_in": 1.0,
                                "pct_out": 2.0,
                            }
                        ],
                    }
                ]
            }
        }
    }
    warned = pool_api.claim_warnings(colon_pool, "1", [("selected", name)])
    assert warned == [
        (
            name,
            f"selected {name}: not_raised 1/1 | top_n_pub: ABL1(1.0/2.0) | +0 more",
        )
    ]


def test_the_judge_packet_narrows_to_the_label_under_test(pool):
    """What quality control may see, asserted as an exact key set.

    A judge handed `candidates` or the retrieval order could prefer a different label,
    and that disagreement would say nothing about whether THIS one is carried. The same
    assertion is made against the R builder in tests/test_two_arm_packet_parity.R.
    """
    from conftest import FakeSources

    packet = pool_api.judge_packet(pool, "1", "enterocyte", FakeSources())
    assert set(packet) == {
        "query",
        "label_under_test",
        "panel",
        "unmeasured_curated_genes",
        "sources",
    }
    assert packet["label_under_test"] == "enterocyte"
    # The whole definition, ordered by publication support, raised or not: VIL1 and
    # FABP2 are detected in most of the cluster AND in most cells outside it, and that
    # they fail to separate is exactly what a judge has to see. The undetected negative
    # marker is left out -- an exclusion nobody detects excludes nothing.
    assert [row["gene"] for row in packet["panel"]] == ["VIL1", "SI", "FABP2", "RBP2"]
    assert set(packet["sources"]) == {"VIL1", "SI", "FABP2", "RBP2"}


def test_the_judge_reads_the_whole_panel_and_asks_for_the_sentences(pool):
    """No row is cut, because cutting rows is what hid the evidence.

    The judge panel used to stop at the fourteen best-published positives. Measured over
    the shipped run that removed 61 of the genes the annotator's own answers rested on --
    10.3% of them, over 29 of 104 clusters: C1QB under Kupffer cell, VSIG4 under alveolar
    macrophage, PKHD1 under cholangiocyte. A judge cannot weigh a row it was never shown,
    and it could not ask. Rows are cheap and sentences are not, so every row ships and the
    sentences past the head are asked for.
    """
    from conftest import FakeSources

    pushed = pool_api.JUDGE_SOURCES_PUSHED
    entry = {
        "markers": [
            {
                "gene": f"POS{index:02d}",
                "polarity": "positive",
                "n_pub": 100 - index,
                "pct_in": 80.0,
                "pct_out": 5.0,
                "avg_log2FC": 1.0,
                "significant": True,
                "specificity": 0.5,
            }
            for index in range(pushed + 6)
        ]
        + [
            {
                "gene": "NEGDEEP",
                "polarity": "negative",
                "n_pub": 1,
                "pct_in": 70.0,
                "pct_out": 4.0,
                "avg_log2FC": 0.9,
                "significant": True,
                "specificity": 0.5,
            }
        ]
    }
    panel = pool_api._judge_panel(entry)
    assert len(panel) == pushed + 7
    assert panel[-1]["gene"] == "NEGDEEP"

    shipped = pool_api._judge_sources("enterocyte", panel, FakeSources())
    # The best-published positives, and the detected exclusion whatever its rank: an
    # exclusion's whole weight is what its sentence does, so the row without the sentence
    # is the half that cannot be read.
    assert set(shipped) == {f"POS{index:02d}" for index in range(pushed)} | {"NEGDEEP"}
    # The rest are on the panel and answerable, which is the point: visible, and askable.
    assert {row["gene"] for row in panel} - set(shipped) == {
        f"POS{index:02d}" for index in range(pushed, pushed + 6)
    }


def test_a_judge_may_ask_only_about_the_labels_it_is_judging(pool):
    """The pull channel cannot become a candidate list.

    A judge that could fetch a rival's sentences would start choosing between identities,
    and its verdict would no longer be about whether THIS label is carried. So the gate is
    the tool set plus the label, not a sentence in the template asking it not to look.
    """
    from conftest import FakeSources

    assert pool_api.JUDGE_TOOLS == ("sources", "gene_across_clusters")
    assert "candidates_with_gene" not in pool_api.JUDGE_TOOLS
    assert "pool_search" not in pool_api.JUDGE_TOOLS

    labels = ("tuft cell", "enterocyte")
    allowed = pool_api.run_judge_tool(
        pool, "1", labels, "sources", {"label": "enterocyte", "genes": ["VIL1"]},
        FakeSources(),
    )
    assert allowed["sources"]["VIL1"]

    refused = pool_api.run_judge_tool(
        pool, "1", labels, "sources",
        {"label": "enteroendocrine cell", "genes": ["CHGA"]}, FakeSources(),
    )
    assert "not under test here" in refused["error"]

    for tool in ("candidates_with_gene", "pool_search"):
        blocked = pool_api.run_judge_tool(pool, "1", labels, tool, {"gene": "PROX1"})
        assert "unknown tool" in blocked["error"]

    # Cross-cluster expression names no identity, so it is allowed: whether a raised marker
    # is a tissue-wide program is a question about these cells.
    across = pool_api.run_judge_tool(pool, "1", labels, "gene_across_clusters",
                                     {"gene": "PROX1"})
    assert [row["cluster_id"] for row in across["clusters"]]


def test_the_judge_packet_carries_the_parent_when_a_refinement_is_under_test(pool):
    from conftest import FakeSources

    packet = pool_api.judge_packet(
        pool, "1", "tuft cell", FakeSources(), parent="enterocyte"
    )
    assert packet["parent_label"] == "enterocyte"
    assert [row["gene"] for row in packet["parent_panel"]] == [
        "VIL1", "SI", "FABP2", "RBP2"
    ]
    assert packet["raised_and_absent_from_parent_panel"] == ["PROX1"]
    # POU2F3 defines tuft cells and sits at 1.0 here: the row that decides this is not
    # one, and the one a raised-only panel would have hidden.
    assert [row["gene"] for row in packet["panel"]] == ["POU2F3", "PROX1"]
