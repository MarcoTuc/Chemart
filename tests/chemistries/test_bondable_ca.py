"""Bondable cellular automata reproduce book figure 10.16 (Banzhaf & Yamamoto 2015, 10.7.3).

The figure is a "periodic table" of the 256 one-dimensional binary CA rules of
width 12, ordered by the value at which their mean polarity settles. FIG_10_16
below is that figure read off the book's own page: 202 rules sit in the columns
-12..+12 and the other 54 are drawn in a block under the table, which is the
(otherwise empty) polarity-0 column, so they are recorded here as 0.
"""

import random

import pytest

from chemart import generate_network
from chemart.chemistries.bondable_ca import (
    Molecule, best_bond, bits, bond_strength, initial_config, periodic_table, polarity,
    run_at, settle, settle_molecule, step, step_molecule,
)

W = 12

# Book figure 10.16, rule -> settled mean polarity (scale -12..+12).
FIG_10_16 = {
    0: -12, 1: -2, 2: -10, 3: -1, 4: -12, 5: -2, 6: -9, 7: -1, 8: -12, 9: -2, 10: -8, 11: 0, 12: -10, 13: -2, 14: -8, 15: 0,
    16: -10, 17: -1, 18: -5, 19: 0, 20: -9, 21: -1, 22: -2, 23: 0, 24: -8, 25: -1, 26: -2, 27: 1, 28: 0, 29: 0, 30: 0, 31: 1,
    32: -12, 33: -2, 34: -10, 35: -1, 36: -12, 37: -2, 38: -9, 39: -1, 40: -12, 41: -2, 42: -8, 43: 0, 44: -10, 45: 0, 46: -8, 47: 0,
    48: -10, 49: -1, 50: 0, 51: 0, 52: -9, 53: -1, 54: -2, 55: 0, 56: -8, 57: 0, 58: 2, 59: 1, 60: 0, 61: 1, 62: 3, 63: 1,
    64: -12, 65: -2, 66: -8, 67: -1, 68: -10, 69: -2, 70: 0, 71: 0, 72: -12, 73: -1, 74: -4, 75: 0, 76: -8, 77: 0, 78: 2, 79: 2,
    80: -8, 81: 0, 82: -2, 83: 1, 84: -8, 85: 0, 86: 0, 87: 1, 88: -4, 89: 0, 90: 4, 91: 2, 92: 2, 93: 2, 94: 4, 95: 2,
    96: -12, 97: -2, 98: -8, 99: 0, 100: -10, 101: 0, 102: 0, 103: 1, 104: -12, 105: 0, 106: -12, 107: 2, 108: -8, 109: 1, 110: 1, 111: 2,
    112: -8, 113: 0, 114: 2, 115: 1, 116: -8, 117: 0, 118: 3, 119: 1, 120: -12, 121: 2, 122: 4, 123: 2, 124: 1, 125: 2, 126: 4, 127: 2,
    128: -12, 129: -4, 130: -8, 131: -3, 132: -12, 133: -4, 134: -6, 135: 0, 136: -12, 137: -1, 138: 0, 139: 8, 140: -10, 141: -2, 142: 0, 143: 8,
    144: -8, 145: -3, 146: -5, 147: 2, 148: -6, 149: 0, 150: 0, 151: 2, 152: -8, 153: 0, 154: 3, 155: 9, 156: 0, 157: 0, 158: 6, 159: 9,
    160: -12, 161: -4, 162: -6, 163: -2, 164: -12, 165: -4, 166: -3, 167: 2, 168: -12, 169: 12, 170: 0, 171: 8, 172: -10, 173: 4, 174: 0, 175: 8,
    176: -6, 177: -2, 178: 0, 179: 0, 180: -3, 181: 2, 182: 5, 183: 5, 184: 0, 185: 8, 186: 6, 187: 10, 188: 8, 189: 8, 190: 8, 191: 10,
    192: -12, 193: -1, 194: -8, 195: 0, 196: -10, 197: -2, 198: 0, 199: 0, 200: 0, 201: 8, 202: 10, 203: 10, 204: 0, 205: 8, 206: 10, 207: 10,
    208: 0, 209: 8, 210: 3, 211: 9, 212: 0, 213: 8, 214: 6, 215: 9, 216: 10, 217: 10, 218: 12, 219: 12, 220: 10, 221: 10, 222: 12, 223: 12,
    224: -12, 225: 12, 226: 0, 227: 8, 228: -10, 229: 4, 230: 8, 231: 8, 232: 0, 233: 12, 234: 12, 235: 12, 236: 0, 237: 12, 238: 12, 239: 12,
    240: 0, 241: 8, 242: 6, 243: 10, 244: 0, 245: 8, 246: 8, 247: 10, 248: 12, 249: 12, 250: 12, 251: 12, 252: 12, 253: 12, 254: 12, 255: 12,
}

# The eight rules the figure draws at an exact half-integer, four rounded toward
# zero and four away from it, so the figure's own tie-breaking is inconsistent.
HALF_INTEGER = {26, 82, 154, 166, 167, 180, 181, 210}


@pytest.fixture(scope="module")
def table():
    return {row["rule"]: row for row in periodic_table(W)}


# --- figure 10.16 ----------------------------------------------------------------
def test_periodic_table_reproduces_figure_10_16(table):
    """Book fig. 10.16: the settled mean polarity of all 256 rules at width 12."""
    assert len(FIG_10_16) == 256
    exact, off = [], {}
    for rule, published in FIG_10_16.items():
        got = table[rule]["mean_polarity"]
        assert abs(got - published) <= 0.5, f"rule {rule}: figure {published}, got {got}"
        if got == published:
            exact.append(rule)
        else:
            off[rule] = got
    assert len(exact) == 226, "the figure's integer columns are hit exactly by 226 rules"
    # Every rule the figure misses by more than nothing has a fractional cycle mean
    # and is drawn at the nearest column; none of them is an integer.
    assert all(v != int(v) for v in off.values())
    assert {r for r, v in off.items() if abs(v - round(v)) == 0.5} == HALF_INTEGER


def test_the_figure_is_ordered_by_settled_mean_polarity(table):
    """'these 256 types can now be ordered into the analogue of a periodic table' (book 10.7.3)."""
    rows = periodic_table(W)
    assert [r["rule"] for r in rows] != list(range(256))
    assert [r["mean_polarity"] for r in rows] == sorted(r["mean_polarity"] for r in rows)
    assert len(rows) == 256
    # The extremes of the -12..+12 scale, and the empty-looking centre column.
    assert {r for r, v in FIG_10_16.items() if v == -12} == {
        0, 4, 8, 32, 36, 40, 64, 72, 96, 104, 106, 120, 128, 132, 136, 160, 164, 168, 192, 224}
    assert len({r for r, v in FIG_10_16.items() if v == 12}) == 20
    assert len({r for r, v in FIG_10_16.items() if v == 0}) == 54


def test_settling_times_fit_the_figure_shading(table):
    """Fig. 10.16 shades each rule by how long its mean polarity took to settle, up to 2048."""
    times = [row["settling_iteration"] for row in table.values()]
    assert max(times) <= 2048
    assert max(times) == 28
    assert table[0]["settling_iteration"] <= 2 and table[255]["settling_iteration"] <= 2


def test_known_wolfram_rules(table):
    """Rules whose behaviour is known independently of the figure."""
    assert initial_config(W) == 0b111111 and bits(initial_config(W), W) == "111111000000"
    # Rule 0 empties the ring, rule 255 fills it.
    assert table[0]["mean_polarity"] == -12 and table[0]["settled_config"] == "0" * 12
    assert table[255]["mean_polarity"] == 12 and table[255]["settled_config"] == "1" * 12
    # Density-preserving rules hold the balanced start, so they settle at 0.
    for rule in (204, 170, 240, 184, 232, 51):        # identity, shifts, traffic, majority, inversion
        assert table[rule]["mean_polarity"] == 0, rule
        assert FIG_10_16[rule] == 0
    # Rule 90 (XOR of the neighbours) and rule 110 (universal).
    assert table[90]["mean_polarity"] == 4 == FIG_10_16[90]
    assert table[110]["mean_polarity"] == pytest.approx(2 / 3) and FIG_10_16[110] == 1
    # Rule 2 leaves a single cell circling the ring forever: polarity 1 - 11 = -10.
    assert table[2]["mean_polarity"] == -10 and table[2]["cycle_length"] == 12
    # Rule 6 alternates between one and two live cells: (-10 - 8) / 2 = -9.
    assert table[6]["mean_polarity"] == -9


@pytest.mark.slow
def test_the_balanced_block_start_is_what_matches_the_figure():
    """The figure pins the initial condition: a contiguous half-ring of ones."""
    def score(config):
        return sum(1 for rule, published in FIG_10_16.items()
                   if abs(settle(rule, W, config)["mean_polarity"] - published) <= 0.5)

    assert score(initial_config(W)) == 256
    rng = random.Random(0)
    balanced = [sum(1 << i for i in rng.sample(range(W), W // 2)) for _ in range(48)]
    others = sorted(score(c) for c in balanced)
    assert max(others) < 200, "no other balanced start comes close to reproducing the figure"
    assert others[len(others) // 2] < 160, "a typical balanced start matches barely half of it"
    assert score(0b101010101010) < 100, "the alternating start does not reproduce it either"


# --- bonding strength (figure 10.15) ---------------------------------------------
def test_bond_strength_is_the_longest_complementary_run():
    """Fig. 10.15: 'bonding strength based on largest contiguous sequence of complementary binary numbers'."""
    ones = (1 << W) - 1
    assert bond_strength(0, ones, W) == W                     # complementary everywhere
    assert bond_strength(0, 1, W) == 1                        # one cell differs
    assert bond_strength(0, 0, W) == 0                        # nothing complementary at any offset
    # A ring equals its own complement under a half-ring rotation.
    assert best_bond(0b000000111111, 0b000000111111, W) == (12, 6, ones)
    assert best_bond(0b101010101010, 0b101010101010, W)[:2] == (12, 1)
    # Sliding a half-block against itself grows the run one cell at a time.
    a = 0b000000111111
    assert [run_at(a, a, W, d)[0] for d in range(W)] == [0, 1, 2, 3, 4, 5, 12, 5, 4, 3, 2, 1]


def test_bond_strength_is_symmetric_and_bounded():
    rng = random.Random(1)
    for _ in range(400):
        a, b = rng.randrange(1 << W), rng.randrange(1 << W)
        s = bond_strength(a, b, W)
        assert s == bond_strength(b, a, W)
        assert 0 <= s <= W
        length, offset, mask = best_bond(a, b, W)
        assert mask.bit_count() == length
        # Every cell of the run really is complementary at that offset.
        for i in range(W):
            if (mask >> i) & 1:
                assert ((a >> i) & 1) != ((b >> ((i + offset) % W)) & 1)


# --- the coupling actually acts ---------------------------------------------------
def _pair(ra, rb):
    a, b = settle(ra, W)["config"], settle(rb, W)["config"]
    return Molecule([ra, rb], [a, b], [best_bond(a, b, W)[1]], W)


def test_uncoupled_atoms_run_independently():
    m = _pair(0, 175)
    assert step_molecule(m, "none") == tuple(
        step(c, rule, W) for c, rule in zip(m.configs, m.rules))


def test_bonding_changes_the_dynamics():
    """Book 10.7.3: once bonded, each CA's neighbourhood is influenced by the other."""
    for ra, rb in ((0, 175), (7, 5)):
        m = _pair(ra, rb)
        assert bond_strength(*m.configs, W) >= 8
        bonded = settle_molecule(m, "replace", 256)[0]
        free = settle_molecule(m, "none", 256)[0]
        assert bonded.configs != free.configs, (ra, rb)
    differ = 0
    for ra in range(0, 256, 17):
        for rb in range(0, 256, 13):
            m = _pair(ra, rb)
            if bond_strength(*m.configs, W) < 8:
                continue
            if settle_molecule(m, "replace", 256)[0].configs != settle_molecule(m, "none", 256)[0].configs:
                differ += 1
    assert differ > 20, "the coupling changes the outcome for many bondable pairs"


def test_polarity_is_the_figure_scale():
    assert polarity(0, W) == -W and polarity((1 << W) - 1, W) == W
    assert polarity(initial_config(W), W) == 0


# --- the reaction network ---------------------------------------------------------
def test_default_network_is_a_closure_of_associations_and_dissociations():
    net = generate_network("bondable-ca", seed=1)
    assert net.status == "complete"
    assert len(net.species) == 23 and len(net.reactions) == 11
    kinds = {}
    for r in net.reactions:
        assert r.rate is None, "the book gives no kinetics for BCA"
        n_in, n_out = sum(r.reactants.values()), sum(r.products.values())
        assert (n_in, n_out) in ((2, 1), (1, 2))
        kinds.setdefault("association" if n_in == 2 else "dissociation", []).append(r)
    assert len(kinds["association"]) == 7 and len(kinds["dissociation"]) == 4
    law = net.extras["interaction_law"]
    assert law["coupling"]["mode"] == "replace" and law["settle_iterations"] == 256
    assert len(net.extras["analysis"]["periodic_table"]) == 256


def test_atoms_are_conserved():
    """Association and dissociation only rearrange atoms, so every atom type is conserved."""
    import numpy as np

    net = generate_network("bondable-ca", seed=1)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    laws = net.extras["conservation"]
    assert laws
    for law in laws:
        m = np.array([law["vector"][s] for s in ids])
        assert np.all(m @ S == 0), law["name"]
    # The seeded atom types are species in their own right.
    assert {"r0", "r90", "r255"} <= {s.id for s in net.species}
    for s in net.species:
        assert s.structure.count("~") % 2 == 0


def test_generation_is_deterministic():
    """The chemistry is deterministic: the seed cannot change the network."""
    a = generate_network("bondable-ca", seed=7)
    b = generate_network("bondable-ca", seed=7)
    c = generate_network("bondable-ca", seed=99)
    assert a.to_dict() == b.to_dict()
    assert [r.to_text() for r in a.reactions] == [r.to_text() for r in c.reactions]


def test_coupling_none_gives_a_different_chemistry():
    """Switching the book's interaction off changes which molecules survive."""
    on = generate_network("bondable-ca", seed=1)
    off = generate_network("bondable-ca", seed=1, coupling="none")
    assert {s.id for s in on.species} != {s.id for s in off.species}


def test_bad_parameters():
    with pytest.raises(ValueError, match="bond_threshold must be at most"):
        generate_network("bondable-ca", w=8, bond_threshold=10)
    with pytest.raises(ValueError, match="Wolfram rule numbers"):
        generate_network("bondable-ca", atom_types=[0, 300])
    with pytest.raises(ValueError, match="non-empty list of integers"):
        generate_network("bondable-ca", atom_types=[])
