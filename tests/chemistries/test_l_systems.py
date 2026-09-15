"""L-systems: book 9.8; Prusinkiewicz & Lindenmayer, The Algorithmic Beauty of Plants (ABOP), chapter 1."""

import math
from collections import defaultdict
from math import comb

import pytest

from chemart import catalog, generate_network
from chemart.chemistries.l_systems import SYSTEMS, LSystem, symbols, word_text

ID = "l-systems"
FIB = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]


def chain(net):
    """The words of a deterministic derivation, in order from the axiom."""
    word = {s.id: s.structure for s in net.species}
    succ = {}
    for r in net.reactions:
        (a,), (b,) = r.reactants, r.products
        assert a not in succ, "a deterministic word has one successor"
        succ[a] = b
    out, current = [], "w0"
    while current is not None and len(out) <= len(net.species):
        out.append(word[current])
        current = succ.get(current)
    return out


def joined(text):
    return word_text(symbols(text))


def count(net, letters):
    return [sum(1 for s in symbols(w) if s in letters) for w in chain(net)]


def test_catalog_offers_every_published_system():
    entry = next(c for c in catalog.load() if c.id == ID)
    choices = next(p for p in entry.params if p.name == "system").choices
    assert set(choices) == set(SYSTEMS) | {"custom"}


def test_parallel_derivation_of_figure_1_3():
    """ABOP 1.2 / fig. 1.3: b, a, ab, aba, abaab, abaababa, abaababaabaab."""
    net = generate_network(ID, system="custom", axiom="b", productions=["a -> ab", "b -> a"], iterations=6)
    assert chain(net) == ["b", "a", "ab", "aba", "abaab", "abaababa", "abaababaabaab"]
    assert net.status == "complete" and len(net.reactions) == 6
    assert net.extras["probabilities"] == [1.0] * 6


def test_algae_lengths_are_fibonacci():
    """ABOP eq. (1.2), sec. 1.9: letter counts and word lengths follow the Fibonacci series."""
    net = generate_network(ID, iterations=10)
    words = chain(net)
    assert words[:5] == ["a", "ab", "aba", "abaab", "abaababa"]
    assert [len(w) for w in words] == FIB
    assert [w.count("a") for w in words] == [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
    assert net.initial_state == {"w0": 1.0}


def test_anabaena_filament():
    """ABOP eq. (1.1) and the derivation printed in sec. 1.2; cell counts are Fibonacci numbers."""
    net = generate_network(ID, system="anabaena", iterations=8)
    published = ["a_r", "a_l b_r", "b_l a_r a_r", "a_l a_l b_r a_l b_r", "b_l a_r b_l a_r a_r b_l a_r a_r"]
    words = chain(net)
    assert words[:5] == [joined(w) for w in published]
    assert [len(symbols(w)) for w in words] == FIB[:9]


@pytest.mark.parametrize("n", range(5))
def test_koch_curve_has_4_to_the_n_segments(n):
    """Koch construction (ABOP fig. 1.1: generator of 4 sides): 4^n segments from one edge, 3 * 4^n for the snowflake."""
    curve = generate_network(ID, system="custom", axiom="F", productions=["F -> F+F--F+F"], angle=60.0, iterations=n)
    assert count(curve, {"F"})[-1] == 4 ** n
    snowflake = generate_network(ID, system="koch-snowflake", iterations=n)
    assert count(snowflake, {"F"}) == [3 * 4 ** k for k in range(n + 1)]
    assert snowflake.extras["turtle"] == {"angle": 60.0}


def test_quadratic_koch_island_and_snowflake():
    """ABOP fig. 1.6 (derivations 0 to 3, 90 degrees) and fig. 1.7b (5 sides per generator)."""
    island = generate_network(ID, system="koch-island", iterations=3)
    assert chain(island)[0] == "F-F-F-F"
    assert chain(island)[1] == "F-F+F+FF-F-F+F-F-F+F+FF-F-F+F-F-F+F+FF-F-F+F-F-F+F+FF-F-F+F"
    assert count(island, {"F"}) == [4 * 8 ** k for k in range(4)]
    assert island.extras["turtle"] == {"angle": 90.0, "published_iterations": 3}
    snow = generate_network(ID, system="quadratic-snowflake", iterations=4)
    assert count(snow, {"F"}) == [5 ** k for k in range(5)]
    dragon = generate_network(ID, system="dragon-curve", iterations=10)
    assert count(dragon, {"F_l", "F_r"}) == [2 ** k for k in range(11)]
    gasket = generate_network(ID, system="sierpinski-gasket", iterations=6)
    assert count(gasket, {"F_l", "F_r"}) == [3 ** k for k in range(7)]


def test_bracketed_plant_words():
    """ABOP fig. 1.24a: one step of F -> F[+F]F[-F]F; brackets replace themselves."""
    net = generate_network(ID, system="plant-a", iterations=2)
    words = chain(net)
    assert words[1] == "F[+F]F[-F]F"
    assert words[2] == "F[+F]F[-F]F[+F[+F]F[-F]F]F[+F]F[-F]F[-F[+F]F[-F]F]F[+F]F[-F]F"
    assert net.extras["turtle"] == {"angle": 25.7, "published_iterations": 5}


def test_growth_matrix_and_pascal_triangle():
    """ABOP 1.9: reading symbols gives the matrix Q; a_i -> a_i a_{i+1} grows like the Pascal triangle."""
    algae = generate_network(ID, reading="symbols", iterations=5)
    ids, R, P = algae.matrices()
    Q = P.toarray()                               # successor letter counts, one column per production
    assert ids == ["a", "b"] and R.toarray().tolist() == [[1, 0], [0, 1]]
    assert Q[:, 0].tolist() == [1, 1] and Q[:, 1].tolist() == [1, 0]   # a -> ab, b -> a
    assert [h.get("a", 0) for h in algae.extras["analysis"]["counts"]] == [1, 1, 2, 3, 5, 8]
    assert any(r.catalysts for r in algae.reactions)

    rules = [f"a_{i} -> a_{i} a_{i + 1}" for i in range(8)]
    pascal = generate_network(ID, system="custom", reading="symbols", axiom="a_0", productions=rules, iterations=7)
    counts = pascal.extras["analysis"]["counts"]
    for k in range(8):
        assert [counts[k].get(f"a_{i}", 0) for i in range(k + 1)] == [comb(k, i) for i in range(k + 1)]

    elongation = generate_network(ID, system="custom", reading="symbols", axiom="F", productions=["F -> FF"], iterations=6)
    assert elongation.extras["analysis"]["lengths"] == [2 ** n for n in range(7)]


def test_symbols_reading_matches_word_lengths():
    words = generate_network(ID, system="plant-d", iterations=4)
    letters = generate_network(ID, system="plant-d", reading="symbols", iterations=4)
    assert [len(w) for w in chain(words)] == letters.extras["analysis"]["lengths"]


def test_signal_propagation():
    """ABOP 1.8: b < a -> b, b -> a moves the letter b to the right."""
    net = generate_network(ID, system="signal-propagation", iterations=12)
    words = chain(net)
    assert words[:5] == ["baaaaaaaa", "abaaaaaaa", "aabaaaaaa", "aaabaaaaa", "aaaabaaaa"]
    assert words[-1] == "aaaaaaaaa" and len(words) == 10
    assert net.extras["stationary"] == {"w9": 1.0}     # the signal has left the string


def test_context_matching_skips_branches_as_in_figure_1_29():
    """ABOP 1.8: BC < S > G[H]M applies to S in ABC[DE][SG[HI[JK]L]MNO]."""
    word = "ABC[DE][SG[HI[JK]L]MNO]"
    applies = generate_network(ID, system="custom", axiom=word, productions=["BC < S > G[H]M -> X"], iterations=1)
    assert chain(applies)[1] == "ABC[DE][XG[HI[JK]L]MNO]"
    for rule in ["BC < S > G[H]N -> X", "DE < S > G -> X", "BC < S > G[I] -> X"]:
        other = generate_network(ID, system="custom", axiom=word, productions=[rule], iterations=1)
        assert other.reactions == [], rule


def test_acropetal_and_basipetal_signals():
    """ABOP fig. 1.30 (#ignore: +-); words derived by hand from the matching rules of sec. 1.8."""
    up = chain(generate_network(ID, system="acropetal-signal", iterations=6))
    assert up == [joined(w) for w in [
        "F_b[+F_a]F_a[-F_a]F_a[+F_a]F_a",
        "F_b[+F_b]F_b[-F_a]F_a[+F_a]F_a",
        "F_b[+F_b]F_b[-F_b]F_b[+F_a]F_a",
        "F_b[+F_b]F_b[-F_b]F_b[+F_b]F_b",
    ]]
    down = chain(generate_network(ID, system="basipetal-signal", iterations=6))
    assert down == [joined(w) for w in [
        "F_a[+F_a]F_a[-F_a]F_a[+F_a]F_b",
        "F_a[+F_a]F_a[-F_a]F_b[+F_a]F_b",
        "F_a[+F_a]F_b[-F_a]F_b[+F_a]F_b",
        "F_b[+F_a]F_b[-F_a]F_b[+F_a]F_b",
    ]]


def test_square_root_growth():
    """ABOP eq. (1.5), fig. 1.33: the 2L-system grows as floor(sqrt(n)) + 4."""
    net = generate_network(ID, system="square-root-growth", iterations=50)
    assert [len(symbols(w)) for w in chain(net)] == [math.isqrt(n) + 4 for n in range(51)]


def test_hogeweg_hesper_structures():
    """ABOP fig. 1.31a: bracketed 2L-system over {0, 1} with #ignore +-F, derivation length 30."""
    lsys = LSystem(SYSTEMS["hogeweg-hesper-a"]["axiom"], SYSTEMS["hogeweg-hesper-a"]["productions"], "+-F")
    words = lsys.derivation(30)
    assert all(set(w) <= set("01F+-[]") for w in words)
    assert len(words[30]) > len(words[0])
    assert any("[" in w for w in words)          # 0 < 0 > 1 -> 1[+F1F1] fires: the structure branches
    net = generate_network(ID, system="hogeweg-hesper-a", iterations=30)
    assert net.status == "complete"
    word = {s.id: s.structure for s in net.species}
    assert word["w0"] == "F1F1F1"
    for r in net.reactions:
        (a,), (b,) = r.reactants, r.products
        assert lsys.derive(symbols(word[a])) == symbols(word[b])
    assert {word[s.id] for s in net.species} <= {"".join(w) for w in words}


def test_stochastic_alternatives_sum_to_one():
    """ABOP 1.7: F -> F[+F]F[-F]F (.33), F[+F]F (.33), F[-F]F (.34), each occurrence chosen independently."""
    net = generate_network(ID, system="stochastic-plant", iterations=2)
    assert net.status == "complete"
    word = {s.id: s.structure for s in net.species}
    first = {word[next(iter(r.products))]: pr for r, pr in zip(net.reactions, net.extras["probabilities"])
             if "w0" in r.reactants}
    assert first == pytest.approx({"F[+F]F[-F]F": 0.33, "F[+F]F": 0.33, "F[-F]F": 0.34})
    total = defaultdict(float, net.extras["stationary"])
    for r, pr in zip(net.reactions, net.extras["probabilities"]):
        total[next(iter(r.reactants))] += pr
    expanded = [s for s, d in zip(net.species, net.extras["analysis"]["steps"]) if d < 2]
    assert len(expanded) == 4
    assert all(total[s.id] == pytest.approx(1.0) for s in expanded)
    # F[+F]F[-F]F has five F: 3^5 combinations, all giving distinct words
    assert sum(1 for r in net.reactions if word[next(iter(r.reactants))] == "F[+F]F[-F]F") == 243
    assert all(r.rate is None for r in net.reactions)   # probabilities are not rates


def test_stochastic_budget_truncates():
    net = generate_network(ID, system="stochastic-plant", iterations=3, max_species=500)
    assert net.status == "truncated"
    assert len(net.species) <= 500 and net.extras["analysis"]["unexpanded"]


def test_stochastic_symbols_reading():
    net = generate_network(ID, system="stochastic-plant", reading="symbols", iterations=3)
    assert len(net.reactions) == 3 and sum(net.extras["probabilities"]) == pytest.approx(1.0)
    # expected number of F: each F becomes 5 F (.33) or 3 F (.67): factor 3.66 per step
    assert net.extras["analysis"]["counts"][3]["F"] == pytest.approx(3.66 ** 3)


def test_length_budget_truncates():
    net = generate_network(ID, system="plant-a", iterations=10, max_length=5000)
    assert net.status == "truncated" and max(net.extras["analysis"]["lengths"]) <= 5000


@pytest.mark.parametrize(
    "given, message",
    [
        (dict(system="custom"), "axiom"),
        (dict(system="custom", axiom="a"), "productions"),
        (dict(productions=["a -> b"]), "only used with system='custom'"),
        (dict(system="custom", axiom="F", productions=["F -> F[+F : 0.5", "F -> F : 0.5"]), "unbalanced"),
        (dict(system="custom", axiom="F", productions=["F -> FF : 0.5", "F -> F : 0.4"]), "sum to"),
        (dict(system="custom", axiom="F", productions=["F -> FF", "F -> F"]), "probability"),
        (dict(system="custom", axiom="ab", productions=["ab -> b"]), "one symbol"),
        (dict(system="custom", axiom="a", productions=["[ -> a"]), "brackets"),
        (dict(system="signal-propagation", reading="symbols"), "context-sensitive"),
    ],
)
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
