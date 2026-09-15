"""Matrix chemistry (book chapter 3, 12.5.2, 13.2; Banzhaf 1993): published tables, counts, closure, ODE."""

import numpy as np
import pytest
from odes import rhs

from chemart import generate_network
from chemart.chemistries.matrix_chemistry import bitstring, operation, reaction_table


def parse_table(text):
    rows = [list(map(int, line.split())) for line in text.strip().splitlines()]
    return {row[0]: row[1:] for row in rows}


# Book table 3.3: N = 4, first folding, operators and strings 1..15.
TABLE_3_3 = parse_table("""
 1  1 0 1  4  5  4  5  0  1  0  1  4  5  4  5
 2  0 1 1  0  0  1  1  4  4  5  5  4  4  5  5
 3  1 1 1  4  5  5  5  4  5  5  5  4  5  5  5
 4  2 0 2  8 10  8 10  0  2  0  2  8 10  8 10
 5  3 0 3 12 15 12 15  0  3  0  3 12 15 12 15
 6  2 1 3  8 10  9 11  4  6  5  7 12 14 13 15
 7  3 1 3 12 15 13 15  4  7  5  7 12 15 13 15
 8  0 2 2  0  0  2  2  8  8 10 10  8  8 10 10
 9  1 2 3  4  5  6  7  8  9 10 11 12 13 14 15
10  0 3 3  0  0  3  3 12 12 15 15 12 12 15 15
11  1 3 3  4  5  7  7 12 13 15 15 12 13 15 15
12  2 2 2  8 10 10 10  8 10 10 10  8 10 10 10
13  3 2 3 12 15 14 15  8 11 10 11 12 15 14 15
14  2 3 3  8 10 11 11 12 14 15 15 12 14 15 15
15  3 3 3 12 15 15 15 12 15 15 15 12 15 15 15
""")

# Book table 3.5: N = 9, first folding, operators and strings 165..179.
TABLE_3_5 = parse_table("""
165 283 287 287 280 281 284 285 283 283 287 287 312 313 316 317
166 347 351 351 344 344 349 349 347 347 351 351 376 376 381 381
167 347 351 351 344 345 349 349 347 347 351 351 376 377 381 381
168 274 278 278 272 274 276 278 274 274 278 278 304 306 308 310
169 275 278 279 280 283 284 287 282 283 286 287 304 307 308 311
170 338 343 343 336 338 341 343 338 338 343 343 376 378 381 383
171 339 343 343 344 347 349 351 346 347 351 351 376 379 381 383
172 283 287 287 280 282 284 286 283 283 287 287 312 314 316 318
173 283 287 287 280 283 284 287 283 283 287 287 312 315 316 319
174 347 351 351 344 346 349 351 347 347 351 351 376 378 381 383
175 347 351 351 344 347 349 351 347 347 351 351 376 379 381 383
176 402 406 406 400 400 406 406 402 402 406 406 432 432 438 438
177 403 406 407 408 409 414 415 410 411 414 415 432 433 438 439
178 466 471 471 464 464 471 471 466 466 471 471 504 504 511 511
179 467 471 471 472 473 479 479 474 475 479 479 504 505 511 511
""")

# Book table 3.7 (= table 13.1): N = 9, first folding, operators and strings 1..27; 0 = elastic.
TABLE_3_7 = parse_table("""
1  1 0 1 0 1 0 1 8 9 8 9 8 9 8 9 0 1 0 1 0 1 0 1 8 9 8 9
2  0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 8 8 9 9 8 8 9 9 8 8 9 9
3  1 1 1 0 1 1 1 8 9 9 9 8 9 9 9 8 9 9 9 8 9 9 9 8 9 9 9
4  0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0
5  1 0 1 1 1 1 1 8 9 8 9 9 9 9 9 0 1 0 1 1 1 1 1 8 9 8 9
6  0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 8 8 9 9 9 9 9 9 8 8 9 9
7  1 1 1 1 1 1 1 8 9 9 9 9 9 9 9 8 9 9 9 9 9 9 9 8 9 9 9
8  2 0 2 0 2 0 2 16 18 16 18 16 18 16 18 0 2 0 2 0 2 0 2 16 18 16 18
9  3 0 3 0 3 0 3 24 27 24 27 24 27 24 27 0 3 0 3 0 3 0 3 24 27 24 27
10 2 1 3 0 2 1 3 16 18 17 19 16 18 17 19 8 10 9 11 8 10 9 11 24 26 25 27
11 3 1 3 0 3 1 3 24 27 25 27 24 27 25 27 8 11 9 11 8 11 9 11 24 27 25 27
12 2 0 2 1 3 1 3 16 18 16 18 17 19 17 19 0 2 0 2 1 3 1 3 16 18 16 18
13 3 0 3 1 3 1 3 24 27 24 27 25 27 25 27 0 3 0 3 1 3 1 3 24 27 24 27
14 2 1 3 1 3 1 3 16 18 17 19 17 19 17 19 8 10 9 11 9 11 9 11 24 26 25 27
15 3 1 3 1 3 1 3 24 27 25 27 25 27 25 27 8 11 9 11 9 11 9 11 24 27 25 27
16 0 2 2 0 0 2 2 0 0 2 2 0 0 2 2 16 16 18 18 16 16 18 18 16 16 18 18
17 1 2 3 0 1 2 3 8 9 10 11 8 9 10 11 16 17 18 19 16 17 18 19 24 25 26 27
18 0 3 3 0 0 3 3 0 0 3 3 0 0 3 3 24 24 27 27 24 24 27 27 24 24 27 27
19 1 3 3 0 1 3 3 8 9 11 11 8 9 11 11 24 25 27 27 24 25 27 27 24 25 27 27
20 0 2 2 1 1 3 3 0 0 2 2 1 1 3 3 16 16 18 18 17 17 19 19 16 16 18 18
21 1 2 3 1 1 3 3 8 9 10 11 9 9 11 11 16 17 18 19 17 17 19 19 24 25 26 27
22 0 3 3 1 1 3 3 0 0 3 3 1 1 3 3 24 24 27 27 25 25 27 27 24 24 27 27
23 1 3 3 1 1 3 3 8 9 11 11 9 9 11 11 24 25 27 27 25 25 27 27 24 25 27 27
24 2 2 2 0 2 2 2 16 18 18 18 16 18 18 18 16 18 18 18 16 18 18 18 16 18 18 18
25 3 2 3 0 3 2 3 24 27 26 27 24 27 26 27 16 19 18 19 16 19 18 19 24 27 26 27
26 2 3 3 0 2 3 3 16 18 19 19 16 18 19 19 24 26 27 27 24 26 27 27 24 26 27 27
27 3 3 3 0 3 3 3 24 27 27 27 24 27 27 27 24 27 27 27 24 27 27 27 24 27 27 27
""")

# Book table 12.3: the organisations of the 4-bit matrix chemistry, IDs 0..51.
TABLE_12_3 = [
    [], [1], [8], [9], [15], [1, 8], [1, 9], [6, 9], [8, 9], [9, 15],
    [6, 9, 15], [7, 9, 15], [9, 14, 15], [11, 13, 15], [1, 8, 9],
    [1, 2, 4, 8], [1, 3, 5, 15], [7, 11, 13, 15], [8, 10, 12, 15], [9, 11, 13, 15], [11, 13, 14, 15],
    [1, 2, 4, 8, 9], [1, 3, 5, 9, 15], [7, 11, 13, 14, 15], [7, 9, 11, 13, 15], [8, 9, 10, 12, 15],
    [9, 11, 13, 14, 15],
    [1, 2, 4, 6, 8, 9], [1, 3, 5, 7, 9, 15], [1, 3, 5, 11, 13, 15], [7, 9, 11, 13, 14, 15],
    [8, 9, 10, 12, 14, 15], [8, 10, 11, 12, 13, 15],
    [1, 3, 5, 7, 11, 13, 15], [1, 3, 5, 9, 11, 13, 15], [6, 7, 9, 11, 13, 14, 15],
    [8, 9, 10, 11, 12, 13, 15], [8, 10, 11, 12, 13, 14, 15],
    [1, 3, 5, 7, 9, 11, 13, 15], [8, 9, 10, 11, 12, 13, 14, 15],
    [1, 2, 3, 4, 5, 8, 10, 12, 15], [1, 2, 3, 4, 5, 8, 9, 10, 12, 15],
    [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15], [1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 15],
    [1, 2, 3, 4, 5, 8, 9, 10, 12, 14, 15], [1, 2, 3, 4, 5, 8, 10, 11, 12, 13, 15],
    [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 15],
    [1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15], [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 15],
    [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15],
    [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    list(range(1, 16)),
]

ALL_4BIT = list(range(1, 16))


def number(species_id):
    return int(species_id[1:])


def made(reaction):
    (s,) = [s for s in reaction.products if reaction.products[s] > reaction.reactants.get(s, 0)]
    return s


def ordered_products(net):
    """{(operator, string): product} recovered from a network's multiset reactions."""
    product = operation(net.params["N"], net.params["folding"], net.params["Theta"])
    out = {}
    for r in net.reactions:
        pair = [number(s) for s, n in r.reactants.items() for _ in range(n)]
        c = number(made(r))
        hits = {(a, b) for a, b in (pair, pair[::-1]) if product(a, b) == c}
        assert len(hits) == r.rate["k"]
        out.update({h: c for h in hits})
    return out


def replication_counts(table, strings):
    """(self-replications, replications) over the given strings, as in book 3.2 and table 3.4."""
    T = table[np.ix_(strings, strings)]
    ids = np.array(strings)
    self_rep = int((np.diag(T) == ids).sum())
    off = ~np.eye(len(ids), dtype=bool)
    rep = int((((T == ids[:, None]) | (T == ids[None, :])) & off).sum())
    return self_rep, rep


# --- the reaction rule -------------------------------------------------------
def test_table_3_3_is_the_network_for_n4():
    net = generate_network("matrix-chemistry", N=4, seed_species=ALL_4BIT)
    assert net.status == "complete" and [s.id for s in net.species] == [f"s{i}" for i in ALL_4BIT]
    table = ordered_products(net)
    for op, row in TABLE_3_3.items():
        for string, c in zip(ALL_4BIT, row):
            assert table.get((op, string), 0) == c, (op, string)


def test_worked_examples_and_bit_order():
    product = operation(4)
    assert bitstring(5, 4) == "1010"                        # book 3.2: s(5) = (1,0,1,0)
    assert product(8, 8) == 8 and product(1, 11) == 1 and product(1, 6) == 4   # eqs. 3.18-3.20
    # Banzhaf 1993 [63]: P = ((0,0),(0,1)) on s' = (1,0,0,1) gives (0,0,0,1).
    assert product(int("0001"[::-1], 2), int("1001"[::-1], 2)) == int("0001"[::-1], 2)


def test_n4_has_4_self_replications_and_76_replications():
    net = generate_network("matrix-chemistry", N=4, seed_species=ALL_4BIT)
    self_rep = [r for r in net.reactions if len(r.reactants) == 1 and r.products == {next(iter(r.reactants)): 3}]
    assert sorted(number(next(iter(r.reactants))) for r in self_rep) == [1, 8, 9, 15]
    assert net.extras["analysis"]["self_replicators"] == ["s1", "s8", "s9", "s15"]
    replications = sum(r.rate["k"] for r in net.reactions if len(r.reactants) == 2 and made(r) in r.reactants)
    assert replications == 76                                # table 3.3 caption
    assert replication_counts(reaction_table(4), ALL_4BIT) == (4, 76)


@pytest.mark.parametrize("folding, expected", [(1, (14, 12028)), (2, (122, 21310)), (3, (18, 11822)), (4, (94, 16830))])
def test_table_3_4_counts_for_n9(folding, expected):
    assert replication_counts(reaction_table(9, folding), list(range(1, 512))) == expected


def test_table_3_5_and_table_3_7_for_n9():
    table = reaction_table(9)
    for op, row in TABLE_3_5.items():
        assert table[op, 165:180].tolist() == row, op
    for op, row in TABLE_3_7.items():
        assert table[op, 1:28].tolist() == row, op
    product = operation(9)
    assert all(product(op, s) == table[op, s] for op in range(0, 512, 7) for s in range(0, 512, 5))


def test_table_12_2_adds_the_destructor_row_and_column():
    table = reaction_table(4)
    assert table[0].tolist() == [0] * 16 and table[:, 0].tolist() == [0] * 16
    assert all(table[op, 1:].tolist() == row for op, row in TABLE_3_3.items())


# --- closure, organisations, dynamics ------------------------------------------
def test_default_is_the_book_closure_from_s1_to_s15():
    net = generate_network("matrix-chemistry")
    assert net.status == "complete"
    assert {number(s.id) for s in net.species} == set(range(1, 28)) - {20, 21, 22, 23}   # book 3.3
    assert net.initial_state == {f"s{i}": pytest.approx(1 / 15) for i in range(1, 16)}
    assert {s.id: s.structure for s in net.species}["s27"] == "110110000"
    table = ordered_products(net)
    strings = [number(s.id) for s in net.species]
    assert all(table.get((a, b), 0) == TABLE_3_7[a][b - 1] for a in strings for b in strings)


def test_eq_3_24_rate_equations():
    net = generate_network("matrix-chemistry")
    ids, f = rhs(net)
    strings = [number(s) for s in ids]
    x = np.random.default_rng(5).dirichlet(np.ones(len(ids)))
    W = {(TABLE_3_7[j][k - 1], j, k) for j in strings for k in strings if TABLE_3_7[j][k - 1]}
    production = np.zeros(len(ids))
    for i, j, k in W:
        production[strings.index(i)] += x[strings.index(j)] * x[strings.index(k)]
    phi = production.sum()                                   # eq. 3.25
    assert f(0.0, x) == pytest.approx(production - x * phi / x.sum(), abs=1e-12)


def test_organisations_of_the_4bit_system():
    net = generate_network("matrix-chemistry", N=4, seed_species=ALL_4BIT)
    ids = [s.id for s in net.species]
    index = {s: i for i, s in enumerate(ids)}
    sets = np.arange(1 << len(ids), dtype=np.int64)
    closed = np.ones(len(sets), dtype=bool)
    produced = np.zeros(len(sets), dtype=np.int64)
    for r in net.reactions:
        inside = np.bitwise_and.reduce([(sets >> index[s]) & 1 for s in r.reactants]).astype(bool)
        bit = index[made(r)]
        closed &= ~(inside & ((sets >> bit) & 1 == 0))
        produced |= np.where(inside, 1 << bit, 0)
    # No reaction consumes anything, so self-maintenance under the dilution flow
    # means every member is produced inside the set.
    organisations = sets[closed & ((produced & sets) == sets)]
    as_sets = {frozenset(number(ids[i]) for i in range(len(ids)) if m >> i & 1) for m in organisations}
    assert int(closed.sum()) == 122
    assert len(as_sets) == 54                                # prose of 12.5.2
    assert {frozenset(s) for s in TABLE_12_3} <= as_sets     # every entry of table 12.3
    assert as_sets - {frozenset(s) for s in TABLE_12_3} == {
        frozenset({1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 15}),
        frozenset({1, 2, 3, 4, 5, 8, 10, 11, 12, 13, 14, 15}),
    }


def test_13_2_1_strings_2_to_7_feed_s1():
    net = generate_network("matrix-chemistry", seed_species=[2, 3, 4, 5, 6, 7])
    assert {number(s.id) for s in net.species} == {1, 2, 3, 4, 5, 6, 7}
    among_seeds = [r for r in net.reactions if "s1" not in r.reactants]
    assert among_seeds and all(made(r) == "s1" for r in among_seeds)


# --- options -------------------------------------------------------------------
def test_budget_truncates_the_closure():
    net = generate_network("matrix-chemistry", max_species=20)
    assert net.status == "truncated" and len(net.species) <= 20


def test_destructor_can_be_made_reactive():
    net = generate_network("matrix-chemistry", N=4, seed_species=[1, 2], destructor_elastic=False)
    lines = net.to_text().splitlines()
    assert "s1 + s2 -> s1 + s2 + s0  [mass-action k=2.0]" in lines   # s(1) on s(2) and s(2) on s(1) both give s(0)
    assert "s0 + s1 -> 2 s0 + s1  [mass-action k=2.0]" in lines


def test_soup_records_fired_reactions():
    kwargs = dict(N=4, seed_species=ALL_4BIT, method="soup", M=300, steps=4000)
    net = generate_network("matrix-chemistry", seed=3, **kwargs)
    assert net == generate_network("matrix-chemistry", seed=3, **kwargs)
    assert net.status == "observed" and all(r.count > 0 for r in net.reactions)
    assert sum(r.count for r in net.reactions) <= 4000
    assert sum(net.initial_state.values()) == 300
    assert sum(net.extras["analysis"]["final_population"].values()) == 300
    full = generate_network("matrix-chemistry", N=4, seed_species=ALL_4BIT)
    rates = {(frozenset(r.reactants.items()), frozenset(r.products.items())): r.rate for r in full.reactions}
    for r in net.reactions:
        assert rates[(frozenset(r.reactants.items()), frozenset(r.products.items()))] == r.rate


@pytest.mark.parametrize("given, message", [
    ({"N": 8}, "perfect square"),
    ({"N": 4, "seed_species": [16]}, "seed_species"),
    ({"seed_species": [0, 1]}, "destructor"),
    ({"seed_species": []}, "seed_species"),
])
def test_invalid_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network("matrix-chemistry", **given)
