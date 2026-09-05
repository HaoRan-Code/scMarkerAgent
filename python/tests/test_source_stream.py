"""Source sentences are a stream the agent can draw from, not a fixed shortlist.

A marker's rows share one publication count, so the first few are simply a few of however
many exist. What is asserted here is that asking again reaches the rest, that nothing is
served twice, and that the order is the same wherever it is computed -- the last of which
is what lets the R arm serve the same sentences without asking this one.
"""

from __future__ import annotations

import pytest

from scmarkeragent.marker_sources import SourceServer, order_key


class Context:
    """A pair with more sentences than one batch, spread over two publication counts."""

    def __init__(self):
        self.rows = [
            {"sentence": f"s{i}", "pmid": str(i), "pmcid": f"PMC{i}", "n_pub": n}
            for i, n in enumerate([9, 9, 9, 9, 9, 4, 4], start=1)
        ]

    def all_records(self, cell_type, gene):
        ordered = sorted(
            self.rows,
            key=lambda row: (-row["n_pub"], order_key("seed", row["pmcid"], row["sentence"])),
        )
        return [
            {"sentence": r["sentence"], "pmid": r["pmid"], "pmcid": r["pmcid"]}
            for r in ordered
        ]

    def records_for_marker(self, cell_type, gene, k=3):
        return self.all_records(cell_type, gene)[:k]


def test_asking_again_reaches_the_sentences_the_first_answer_left_behind():
    server = SourceServer(Context(), batch=3, max_batches=0)
    first = server.take("enterocyte", "SI")
    second = server.take("enterocyte", "SI")
    assert [r["sentence"] for r in first["sources"]] != [
        r["sentence"] for r in second["sources"]
    ]
    assert not set(r["sentence"] for r in first["sources"]) & set(
        r["sentence"] for r in second["sources"]
    )
    assert first["remaining"] == 4 and not first["exhausted"]
    assert second["remaining"] == 1


def test_the_stream_ends_and_says_so():
    server = SourceServer(Context(), batch=3, max_batches=0)
    seen = []
    for _ in range(4):
        answer = server.take("enterocyte", "SI")
        seen += [r["sentence"] for r in answer["sources"]]
    assert sorted(seen) == sorted(f"s{i}" for i in range(1, 8))
    assert answer["exhausted"] and answer["remaining"] == 0


def test_higher_publication_counts_are_served_before_lower_ones():
    server = SourceServer(Context(), batch=3, max_batches=0)
    first = [r["sentence"] for r in server.take("enterocyte", "SI")["sources"]]
    # s6 and s7 carry the lower count, so nothing from them reaches the first batch while
    # five better-corroborated sentences are still unread.
    assert "s6" not in first and "s7" not in first


def test_sentences_already_on_the_page_are_not_served_again():
    server = SourceServer(Context(), batch=3, max_batches=0)
    opening = server.opening("enterocyte", "SI")
    answer = server.take("enterocyte", "SI")
    assert not set(r["sentence"] for r in opening) & set(
        r["sentence"] for r in answer["sources"]
    )


def test_a_pair_can_only_be_drawn_from_so_many_times():
    server = SourceServer(Context(), batch=1, max_batches=2)
    server.take("enterocyte", "SI")
    server.take("enterocyte", "SI")
    stopped = server.take("enterocyte", "SI")
    assert stopped["sources"] == []
    assert stopped["limit_reached"] and stopped["remaining"] == 5


def test_one_pair_running_out_does_not_affect_another():
    server = SourceServer(Context(), batch=3, max_batches=1)
    server.take("enterocyte", "SI")
    other = server.take("enterocyte", "VIL1")
    assert len(other["sources"]) == 3


@pytest.mark.parametrize("seed", ["scmarkeragent-v1", "another-draw"])
def test_the_order_is_fixed_by_the_seed_rather_than_by_the_sort(seed):
    keys = [order_key(seed, "PMC1", "a"), order_key(seed, "PMC1", "a")]
    assert keys[0] == keys[1]
    assert order_key(seed, "PMC1", "a") != order_key(seed, "PMC2", "a")


def test_a_different_seed_is_a_different_draw():
    a = order_key("one", "PMC1", "sentence")
    b = order_key("two", "PMC1", "sentence")
    assert a != b
