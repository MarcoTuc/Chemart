"""Conrad's lock-and-key enzymatic processor, as realised in the wet lab by
Zauner & Conrad, *Enzymatic Computing*, Biotechnol. Prog. 17(3):553-559 (2001).

What is reproduced here: the assay composition of the Materials and Methods
section, Table 1 (which milieu states must be active for each operation and the
signal strength each one needs), the Results statement that AND, OR, XOR, NAND
and NOR are realisable while NXOR is not, the Figure 7 argument that a
monotonic element cannot do the XOR, and the Figure 9 result that the tabletop
device classified all 135 presented patterns correctly.

No absorbance is asserted: the paper publishes its response surface as a
picture, so the response levels used below are ordinal.
"""

import math

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.conrad_enzymatic import (
    KEYS,
    LOCKS,
    OPERATIONS,
    PATTERNS,
    SHAPE_BITS,
    STATES,
    classify,
    implementable,
    match,
    milieu_state,
    signal_strength,
    truth_table,
)

ID = "conrad-enzymatic"

#: Table 1 of the paper. Output states are listed for the input states in the
#: table's own column order (I1 I2 = 11, 01, 10, 00), followed by Delta s.
TABLE_1 = {
    "AND": ("1000", lambda r: r["c"] - max(r["a"], r["b"])),
    "OR": ("1110", lambda r: min(r["b"], r["c"]) - r["a"]),
    "XOR": ("0110", lambda r: r["b"] - max(r["a"], r["c"])),
    "NAND": ("0111", lambda r: min(r["a"], r["b"]) - r["c"]),
    "NOR": ("0001", lambda r: r["a"] - max(r["b"], r["c"])),
    "NXOR": ("1001", lambda r: min(r["a"], r["c"]) - r["b"]),
}
TABLE_1_COLUMNS = ("11", "01", "10", "00")


def texts(net):
    return set(net.to_text().splitlines())


def bump(x, a=4.0, b=1.0):
    """Response levels at 0, 1 and 2 encoding units of signal.

    A stimulation x suppression product (1 + a u) exp(-b u): it rises to a
    maximum and then falls below its zero-signal value, the shape the paper
    reports for MDH along each ion axis ("possibly due to opposition between
    the stimulatory and suppressive effects"). Ordinal, not absorbance.
    """
    f = lambda u: (1.0 + a * u) * math.exp(-b * u)
    return {"a": f(0.0), "b": f(x), "c": f(2.0 * x)}


# --- the published assay -----------------------------------------------------
def test_default_network_is_the_mgcl2_device_of_figures_8_and_9():
    net = generate_network(ID)
    assert net.status == "complete"
    assert [s.id for s in net.species] == [
        "MDH", "MDH_NAD", "MDH_NAD_MAL",
        "MDH_Mg", "MDH_Mg_NAD", "MDH_Mg_NAD_MAL",
        "NAD", "MAL", "OAA", "NADH", "H", "Mg",
    ]
    assert len(net.reactions) == 12
    # two 1-signals put two encoding units of Mg2+ into the cuvette
    assert net.initial_state == {"NAD": 1.29, "MAL": 5.41, "Mg": pytest.approx(144.76)}
    assert net.extras["buffered"] == ["H"]


def test_the_default_composition_is_the_published_protocol():
    """Materials and Methods, Signal Processor: 2 x 0.8 mL signal + 0.5 mL enzyme."""
    net = generate_network(ID)
    d = net.extras["experiment"]["signal_processor"]
    total = 2 * d["signal_volume_mL"] + d["enzyme_volume_mL"]
    assert total == pytest.approx(2.1)
    assert net.params["nad_mM"] == pytest.approx(
        d["nad_in_enzyme_solution_mM"] * d["enzyme_volume_mL"] / total, abs=5e-3)
    assert net.params["malate_mM"] == pytest.approx(
        d["malate_in_signal_mM"] * 2 * d["signal_volume_mL"] / total, abs=5e-3)
    assert net.params["encoding_mM"]["Mg"] == pytest.approx(
        d["mgcl2_in_1_signal_mM"] * d["signal_volume_mL"] / total, abs=5e-3)


def test_the_turnover_is_the_assayed_reaction_and_carries_no_kinetics():
    net = generate_network(ID)
    assert "MDH_NAD_MAL -> MDH + OAA + NADH + H" in texts(net)
    assert "MDH_Mg_NAD_MAL -> MDH_Mg + OAA + NADH + H" in texts(net)
    assert net.extras["experiment"]["overall_reaction"] == \
        "L-malate + NAD+ -> oxaloacetate + NADH + H+"
    # the paper is phenomenological and publishes no rate constants
    assert all(r.rate is None for r in net.reactions)


def test_one_signal_pattern_sets_the_milieu_state():
    for inputs, units in [("00", 0), ("01", 1), ("10", 1), ("11", 2)]:
        net = generate_network(ID, inputs=inputs)
        assert net.extras["analysis"]["milieu_state"] == milieu_state(inputs)
        assert net.initial_state["Mg"] == pytest.approx(units * 72.38)


# --- lock and key ------------------------------------------------------------
def test_only_complementary_shapes_bind():
    """Book 11.2.1: recognition is shape complementarity, nothing else binds."""
    net = generate_network(ID, encoding_mM={"Mg": 20.0, "Ca": 40.0})
    scores = net.extras["interaction_law"]["match_scores"]
    perfect = {"coenzyme/NAD", "substrate/MAL", "modulator/Mg"}
    for pair, score in scores.items():
        site, ligand = pair.split("/")
        assert score == match(KEYS[ligand], LOCKS[site])
        if pair in perfect:
            assert score == SHAPE_BITS
        elif pair == "modulator/Ca":
            assert score == 5
        else:
            assert score <= 4, f"{pair} would bind spuriously"
    # exactly the recognised pairs appear as binding reactions
    assert "MDH + NAD -> MDH_NAD" in texts(net)
    assert "MDH_NAD + MAL -> MDH_NAD_MAL" in texts(net)
    assert "MDH + Mg -> MDH_Mg" in texts(net) and "MDH + Ca -> MDH_Ca" in texts(net)
    for bad in ("MDH + MAL ->", "MDH_NAD + NADH ->", "MDH_NAD + OAA ->", "MDH + NADH ->"):
        assert not any(t.startswith(bad) for t in texts(net))


def test_calcium_is_the_weaker_key_and_drops_out_at_perfect_stringency():
    """Results: the signal strength is stronger for MgCl2 alone than for CaCl2 alone."""
    scores = generate_network(ID).extras["interaction_law"]["match_scores"]
    assert scores["modulator/Ca"] < scores["modulator/Mg"] == SHAPE_BITS
    with pytest.raises(ValueError, match="Ca does not recognise"):
        generate_network(ID, encoding_mM={"Mg": 20.0, "Ca": 40.0}, match_bits=6)
    assert len(generate_network(ID, match_bits=6).reactions) == 12


def test_two_carriers_give_the_lattice_of_conformers():
    """Both ions "either alone or in combination": 4 conformers, 2^2 milieu contexts."""
    net = generate_network(ID, encoding_mM={"Mg": 20.0, "Ca": 40.0}, operation="XOR")
    conformers = [s.id for s in net.species if s.id.startswith("MDH") and "NAD" not in s.id]
    assert conformers == ["MDH", "MDH_Mg", "MDH_Ca", "MDH_Mg_Ca"]
    assert len(net.species) == 19 and len(net.reactions) == 28
    assert {"MDH_Mg + Ca -> MDH_Mg_Ca", "MDH_Ca + Mg -> MDH_Mg_Ca"} <= texts(net)
    # the encoding the paper names for the XOR
    assert {"operation": "XOR", "encoding_mM": {"Mg": 20.0, "Ca": 40.0}} in \
        net.extras["analysis"]["published_encodings"]
    assert net.initial_state["Mg"] == 40.0 and net.initial_state["Ca"] == 80.0


# --- Table 1 -----------------------------------------------------------------
@pytest.mark.parametrize("operation", OPERATIONS)
def test_truth_tables_are_table_1(operation):
    outputs, _ = TABLE_1[operation]
    assert truth_table(operation) == {
        pattern: int(bit) for pattern, bit in zip(TABLE_1_COLUMNS, outputs)
    }
    assert generate_network(ID, operation=operation).extras["analysis"]["truth_table"] == \
        truth_table(operation)


@pytest.mark.parametrize("operation", OPERATIONS)
def test_signal_strength_is_table_1(operation):
    """Delta s = min(r over the active states) - max(r over the rest)."""
    _, delta_s = TABLE_1[operation]
    rng = np.random.default_rng(11)
    for _ in range(200):
        r = dict(zip(STATES, rng.random(3)))
        assert signal_strength(operation, r) == pytest.approx(delta_s(r))
        assert implementable(operation, r) == (delta_s(r) > 0)


def test_the_five_realisable_operations_appear_and_nxor_never_does():
    """Figures 5A-C and 6A-B realise AND, OR, XOR, NAND and NOR; NXOR is not
    implementable with a high-response-is-active readout."""
    found = set()
    for x in np.linspace(0.02, 4.0, 400):
        found |= {op for op in OPERATIONS if implementable(op, bump(float(x)))}
    assert found == {"AND", "OR", "XOR", "NAND", "NOR"}


def test_a_monotonic_element_cannot_do_the_xor():
    """Figure 7: the XOR is linearly inseparable, so a monotonic response would
    need two thresholds; the nonmonotonic response makes it separable."""
    rng = np.random.default_rng(3)
    for _ in range(200):
        levels = np.sort(rng.random(3))
        rising = dict(zip(STATES, map(float, levels)))
        falling = dict(zip(STATES, map(float, levels[::-1])))
        assert not implementable("XOR", rising) and not implementable("XOR", falling)
        assert not implementable("NXOR", rising) and not implementable("NXOR", falling)
        # the linearly separable operations need only a monotonic response
        assert implementable("AND", rising) and implementable("OR", rising)
        assert implementable("NOR", falling) and implementable("NAND", falling)
    peaked = bump(1.5)          # the maximum lies at the single-1-signal state
    assert peaked["b"] > peaked["a"] and peaked["b"] > peaked["c"]
    assert implementable("XOR", peaked)


def test_nxor_would_need_an_interior_minimum():
    """A response that rises then falls never dips at the middle milieu state."""
    rng = np.random.default_rng(5)
    for _ in range(300):
        r = bump(float(rng.uniform(0.01, 6.0)), a=float(rng.uniform(1.5, 8.0)))
        assert r["b"] >= min(r["a"], r["c"])
        assert signal_strength("NXOR", r) <= 0


# --- the device result (Figure 9) -------------------------------------------
def test_the_device_classifies_all_135_presentations():
    net = generate_network(ID)
    published = net.extras["analysis"]["published_classification"]
    assert published["operation"] == "XOR"
    assert published["signalling_substance"] == "MgCl2"
    assert published["by_pattern"] == {"00": 45, "01/10": 46, "11": 44}
    assert published["presentations"] == published["correct"] == 135
    assert published["accuracy"] == 1.0
    assert published["readout_time_s"] == 10.0

    # Figure 9: the single-1 patterns respond high, 00 and 11 low, so one
    # threshold separates the three distributions. Levels are ordinal.
    response = {"a": 0.0, "b": 1.0, "c": 0.0}
    assert implementable("XOR", response)
    outputs = classify("XOR", response, threshold=0.5)
    assert outputs == {"00": 0, "01": 1, "10": 1, "11": 0}
    presented = ["00"] * 45 + ["01"] * 23 + ["10"] * 23 + ["11"] * 44
    assert len(presented) == 135
    assert sum(outputs[p] == p.count("1") % 2 for p in presented) == 135


def test_every_pattern_maps_to_a_milieu_state():
    assert [milieu_state(p) for p in PATTERNS] == ["a", "b", "b", "c"]


# --- bookkeeping -------------------------------------------------------------
def test_matter_is_conserved_through_the_recognition_cycle():
    net = generate_network(ID, encoding_mM={"Mg": 20.0, "Ca": 40.0})
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    laws = {law["name"]: law["vector"] for law in net.extras["conservation"]}
    assert set(laws) == {"enzyme", "nicotinamide dinucleotide", "C4 dicarboxylate",
                         "Mg (free and bound)", "Ca (free and bound)"}
    for name, vector in laws.items():
        m = np.array([vector.get(i, 0) for i in ids], dtype=float)
        assert np.allclose(S.T @ m, 0.0), name
    assert laws["nicotinamide dinucleotide"]["MDH_Mg_NAD_MAL"] == 1
    assert laws["C4 dicarboxylate"] == {
        "MDH_NAD_MAL": 1, "MDH_Mg_NAD_MAL": 1, "MDH_Ca_NAD_MAL": 1,
        "MDH_Mg_Ca_NAD_MAL": 1, "MAL": 1, "OAA": 1,
    }


def test_the_network_does_not_depend_on_the_seed():
    assert generate_network(ID, seed=1).to_dict() == generate_network(ID, seed=2).to_dict() | {"seed": 1}


@pytest.mark.parametrize(
    "given, message",
    [
        (dict(encoding_mM={}), "non-empty"),
        (dict(encoding_mM={"Na": 20.0}), "unknown signal carrier"),
        (dict(encoding_mM={"Mg": -1.0}), ">= 0"),
        (dict(encoding_mM={"Mg": "20"}), "concentration in mM"),
        (dict(encoding_mM={"Mg": 20.0, "Ca": 40.0}, match_bits=6), "Ca does not recognise"),
        (dict(match_bits=4), "match_bits"),
        (dict(malate_mM=-1.0), "malate_mM"),
    ],
)
def test_rejects_impossible_experiments(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
