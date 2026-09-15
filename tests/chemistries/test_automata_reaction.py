"""The automata reaction reproduces Dittrich & Banzhaf (1998), Artificial Life 4(2):203-220.

Reaction tables and examples are from the preprint (tables 1-2, figs. 3, 4, 9, 10);
the machine itself is a port of the authors' C source autoreac-1.0.
"""

import random
from collections import Counter

import numpy as np
import pytest
from odes import rhs

from chemart import generate_network
from chemart.chemistries.automata_reaction import automata, disassemble, make_react, word_id


# --- the machine ---------------------------------------------------------------
# Paper table 2: s1 (operator) + s2 (operand) => s3.
TABLE_2 = [
    (0x101D662B, 0xCB5D5C2D, 0xCB5D5606),  # typical random reaction
    (0xF82710FB, 0x73FBC890, 0x73FBC891),
    (0x98AC35F1, 0x4C22BC78, 0x4C22BC7C),
    (0xC1126A24, 0xF9C8A25E, 0xC1C8A25E),  # with TDIR
    (0x9A716B98, 0x9878FA02, 0x9878FA02),  # passive replication
    (0x9A716B98, 0x00000000, 0x00000000),
    (0x1E1CA260, 0x078F21FF, 0x1E1CA260),  # active replication
    (0x1E1CA260, 0x00000000, 0x1E1CA260),
]


@pytest.mark.parametrize("table", [1, 2])
@pytest.mark.parametrize("s1, s2, s3", TABLE_2)
def test_paper_table_2(table, s1, s2, s3):
    assert automata(s1, s2, table) == s3


def test_programs_as_listed_in_the_paper():
    # table 1 (the worked example) and the program listings of figs. 3, 4 and 10
    assert disassemble(0x101D662B) == ["EXOR", "SETP 0110", "CPON", "NOP", "MOV", "ID", "MOV"]
    assert disassemble(0x7240A7EF) == ["AND", "ID", "CPOFF", "NOT", "ID", "TDIR", "SETP 0111"]
    assert disassemble(0x1E1CA260) == ["ID", "CPON", "SETP 1010", "EXOR", "MOV", "ID", "MOV"]
    assert disassemble(0x1E64A24E, 2) == ["ID", "TDIR", "SETP 1010", "TDIR", "CPON", "ID", "MOV"]
    assert disassemble(0x0000000A, 1)[0] == "NOT" and disassemble(0x0000000A, 2)[0] == "EQ"


def reaction_table(words, table, forbid=False):
    """Paper section 4: RT(i, j) = k if s_i + s_j => s_k, None if elastic, '*' if outside."""
    react = make_react("automata", table, forbid)
    index = {w: i for i, w in enumerate(words)}
    out = []
    for a in words:
        row = []
        for b in words:
            s3 = react(a, b)
            row.append(None if s3 is None else index.get(s3, "*"))
        out.append(row)
    return out


def test_fig3_organization_reaction_table():
    words = [0x7240A7EF, 0x7240A7EA, 0x7240A7EB, 0x7240A7EE]
    assert reaction_table(words, 1) == [[2, 1, 2, 1], [0, 0, 0, 0], [3, 0, 3, 0], [1, 1, 1, 1]]


def test_fig4_active_replicators():
    words = [0x1E1CA260, 0x1E1CA261, 0x1E1CA264, 0x1011A261]
    assert reaction_table(words, 1) == [[i] * 4 for i in range(4)]
    rng = random.Random(0)
    for _ in range(50):
        assert automata(0x1E1CA260, rng.getrandbits(32), 1) == 0x1E1CA260


def test_fig9_converged_soup_with_filter_f1():
    # generation 3400 of the table-2 run with exact replication forbidden; '-' is elastic
    words = [0x104E264B, 0x104E264A, 0x104F264B, 0x104F264A, 0x140E264B, 0x140E264A, 0x140F264B, 0x140F264A]
    paper = """
        01 -  03 02 05 04 07 06
        -  00 03 02 05 04 07 06
        01 00 03 -  05 04 07 06
        01 00 -  02 05 04 07 06
        -  00 -  02 -  -  -  06
        01 -  03 -  -  -  07 -
        -  00 -  02 -  04 -  -
        01 -  03 -  05 -  -  -
    """
    expected = [[None if c == "-" else int(c) for c in line.split()] for line in paper.strip().splitlines()]
    assert reaction_table(words, 2, forbid=True) == expected


def test_fig10_cross_over():
    assert automata(0x1E64A24E, 0x00000000, 2) == 0x0000024E
    assert automata(0x1064B246, 0x1E64A24E, 2) == 0x1E64A246   # No. 7 + No. 0 => No. 4
    assert automata(0x1E64A24E, 0x1064B246, 2) == 0x1064B24E   # No. 0 + No. 7 => No. 2


def test_passive_and_active_replicators_are_common_and_rare():
    # paper 3.2: about 30% passive replicators, about 0.004% active ones
    rng = random.Random(3)
    strings = [rng.getrandbits(32) for _ in range(3000)]
    operands = [rng.getrandbits(32) for _ in range(3)]
    passive = sum(all(automata(s, o) == o for o in operands) for s in strings) / len(strings)
    active = sum(all(automata(s, o) == s for o in operands) for s in strings) / len(strings)
    assert 0.15 < passive < 0.45
    assert active < 0.002


# --- organizations as closures -------------------------------------------------
@pytest.mark.parametrize("words, table, forbid", [
    (["7240a7ef", "7240a7ea", "7240a7eb", "7240a7ee"], 1, False),                     # fig. 3
    (["1e1ca260", "1e1ca261", "1e1ca264", "1011a261"], 1, False),                     # fig. 4
])
def test_published_organizations_are_closed(words, table, forbid):
    net = generate_network("automata-reaction", method="closure", words=words, code_table=table,
                           forbid_exact_replication=forbid)
    assert net.status == "complete"
    assert {s.id for s in net.species} == {f"w{w}" for w in words}
    assert all(r.catalysts == r.reactants for r in net.reactions), "both reactants are catalysts"


def test_rates_give_the_catalytic_network_equation():
    # paper eqs. 1-2 on the fig. 3 organization: dx_k = sum_ij k_ij^k x_i x_j - x_k sum_ijk k_ij^k x_i x_j
    words = [0x7240A7EF, 0x7240A7EA, 0x7240A7EB, 0x7240A7EE]
    net = generate_network("automata-reaction", method="closure", words=[f"{w:08x}" for w in words])
    ids, f = rhs(net)
    x = np.random.default_rng(0).uniform(0.1, 1.0, len(ids))
    x /= x.sum()
    conc = dict(zip(ids, x))
    production = Counter()
    for a in words:
        for b in words:
            production[word_id(automata(a, b))] += conc[word_id(a)] * conc[word_id(b)]
    flux = sum(production.values())
    ours = dict(zip(ids, f(0.0, x)))
    for s in ids:
        assert ours[s] == pytest.approx(production[s] - conc[s] * flux, abs=1e-12), s


# --- the soup ------------------------------------------------------------------
def test_and_reaction_goes_extinct():
    # fig. 2: the soup is exploited by the lethal string 00000000
    net = generate_network("automata-reaction", seed=0, mechanism="and", M=1000, generations=10)
    final = net.extras["final_state"]
    assert final.get("w00000000", 0) >= 990
    assert net.extras["analysis"]["diversity"][-1] <= 0.01


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_small_soup_converges_to_an_organization(seed):
    # fig. 3, M = 100: diversity and innovativity drop until a few strings dominate
    net = generate_network("automata-reaction", seed=seed, M=100, generations=150)
    a = net.extras["analysis"]
    assert a["diversity"][0] > 0.9 and a["diversity"][-1] <= 0.1
    assert sum(a["innovativity"][:10]) > 1.0 and sum(a["innovativity"][-10:]) == 0
    top = [int(w[1:], 16) for w in list(net.extras["final_state"])[:2]]
    if len(top) == 2:
        assert bin(top[0] ^ top[1]).count("1") <= 8, "surviving strings are syntactically similar"


def test_filter_f1_makes_replication_elastic():
    net = generate_network("automata-reaction", seed=4, M=300, generations=5, code_table=2,
                           forbid_exact_replication=True)
    for r in net.reactions:
        new =Counter(r.products) - Counter(r.reactants)
        assert sum(new.values()) == 1 and not set(new) & set(r.reactants)
    assert all(p < 1.0 for p in net.extras["analysis"]["productivity"])
    free = generate_network("automata-reaction", seed=4, M=300, generations=5, code_table=2)
    assert all(p == 1.0 for p in free.extras["analysis"]["productivity"])


def test_observed_network_rates_count_ordered_pairs():
    net = generate_network("automata-reaction", seed=2, M=200, generations=3)
    assert net.status == "observed" and net.outflow == "constant-total"
    assert sum(net.initial_state.values()) == 200
    for r in net.reactions:
        a, b = [int(s[1:], 16) for s in Counter(r.reactants).elements()]
        (s3,) = (Counter(r.products) - Counter(r.reactants)).elements()
        k = sum(1 for x, y in {(a, b), (b, a)} if word_id(automata(x, y)) == s3)
        assert r.rate == {"law": "mass-action", "k": float(k)} and r.count >= 1
    species = {s.id: s.structure for s in net.species}
    assert all(int(bits, 2) == int(sid[1:], 16) and len(bits) == 32 for sid, bits in species.items())


def test_bad_parameters():
    with pytest.raises(ValueError, match="32-bit"):
        generate_network("automata-reaction", words=["xyz"])
    with pytest.raises(ValueError, match="code_table"):
        generate_network("automata-reaction", mechanism="and", code_table=2)
