"""Conrad's lock-and-key enzymatic processor (book 11.2.1, 17.1.2). Catalog id: conrad-enzymatic.

The wet realisation of Conrad's lock-key paradigm is Zauner & Conrad,
*Enzymatic Computing*, Biotechnol. Prog. 17(3):553-559 (2001): mitochondrial
malate dehydrogenase (MDH) reads a two-bit input pattern carried by Mg2+ and
Ca2+ and classifies it by its catalytic activity, measured as NADH absorbance
at 339 nm. This module returns the in-silico network of that experiment.

Molecules carry a 6-bit shape: ligands a *key*, the enzyme three *locks*
(coenzyme, substrate and modulator sites). A ligand binds a site when at least
`match_bits` positions are complementary - the book's lock-and-key criterion
(11.2.1). The ions bound at the modulator site select a conformer, each
conformer runs the same recognition-and-turnover cycle, and the activity that
results fuses the milieu signals nonlinearly; that nonmonotonic fusion is what
lets one enzyme perform the linearly inseparable XOR.

The response levels themselves are measured, not modelled: Figure 1's surface
was published as a picture, so nothing here invents an absorbance. What is
reproduced is the classification law of Table 1 (which milieu states must be
active, and the signal strength each operation needs).
"""

from __future__ import annotations

from itertools import combinations

from chemart.helpers.explicit import network
from chemart.network import Species

SHAPE_BITS = 6

#: Chemart's shape encoding of the molecules of the MDH assay: ligand keys and
#: enzyme-site locks. A perfect key is the bit-inverse of its lock.
KEYS = {
    "Mg": "111000",     # signal carrier, perfect complement of the modulator site
    "Ca": "011000",     # signal carrier, one bit off: the weaker of the two
    "NAD": "110100",    # coenzyme
    "MAL": "101010",    # L-malate, the substrate
    "NADH": "110111",   # reduced coenzyme, released
    "OAA": "011010",    # oxaloacetate, released
}
LOCKS = {"coenzyme": "001011", "substrate": "010101", "modulator": "000111"}

#: The two signal carriers of the paper, in the order used for conformer names.
IONS = ("Mg", "Ca")
LIGANDS = ("NAD", "MAL", "Mg", "Ca", "NADH", "OAA")

NAMES = {
    "Mg": "Mg2+ signal carrier",
    "Ca": "Ca2+ signal carrier",
    "NAD": "NAD+ coenzyme",
    "MAL": "L-malate substrate",
    "NADH": "NADH, read at 339 nm",
    "OAA": "oxaloacetate",
    "H": "proton released on turnover (glycine/NaOH buffer, pH 10)",
}

#: Milieu states: a = no 1-signal, b = one 1-signal, c = two 1-signals.
STATES = ("a", "b", "c")
PATTERNS = ("00", "01", "10", "11")
OPERATIONS = ("AND", "OR", "XOR", "NAND", "NOR", "NXOR")

#: Table 1 of Zauner & Conrad (2001): the milieu states that must give an
#: active (high) output for each of the six commutative two-bit operations.
ACTIVE_STATES = {
    "AND": ("c",),
    "OR": ("b", "c"),
    "XOR": ("b",),
    "NAND": ("a", "b"),
    "NOR": ("a",),
    "NXOR": ("a", "c"),
}
#: The same table's signal strength, Delta s (last column), as written there.
SIGNAL_STRENGTH = {
    "AND": "r(c) - max(r(a), r(b))",
    "OR": "min(r(b), r(c)) - r(a)",
    "XOR": "r(b) - max(r(a), r(c))",
    "NAND": "min(r(a), r(b)) - r(c)",
    "NOR": "r(a) - max(r(b), r(c))",
    "NXOR": "min(r(a), r(c)) - r(b)",
}

#: The published protocol (paper, Materials and Methods / Results).
EXPERIMENT = {
    "source": "K.-P. Zauner and M. Conrad, Enzymatic Computing, "
              "Biotechnol. Prog. 17(3):553-559, 2001 (doi:10.1021/bp010004n)",
    "enzyme": "mitochondrial malate dehydrogenase, pig heart (ICN Biomedicals); "
              "6 uL/mL of the suspension in the enzyme solution - its molar "
              "concentration is not published",
    "overall_reaction": "L-malate + NAD+ -> oxaloacetate + NADH + H+",
    "readout": "NADH absorbance at 339 nm, 1 cm light path",
    "response_surface": {
        "assays": 36,
        "design": "6 x 6 MgCl2 x CaCl2 factor levels, run as 7 blocks of 6 cuvettes",
        "ion_levels_uL": "0; 25 of 2 M; 25, 50, 75, 100 of 4 M, made up to 300 uL",
        "assay_volume_mL": 3.0,
        "readout_time_s": 300.0,
        "buffer": "0.12 M glycine/NaOH, pH 10",
        "temperature_C": 29.0,
        "shape": "rises to a maximum and then falls along both ion axes "
                 "('convex', i.e. strictly nonmonotonic)",
    },
    "signal_processor": {
        "signal_volume_mL": 0.8,
        "enzyme_volume_mL": 0.5,
        "malate_in_signal_mM": 7.1,
        "nad_in_enzyme_solution_mM": 5.4,
        "mgcl2_in_1_signal_mM": 190.0,
        "readout_time_s": 10.0,
    },
}

#: Results the paper states in numbers.
PUBLISHED_CLASSIFICATION = {
    "operation": "XOR",
    "signalling_substance": "MgCl2",
    "presentations": 135,
    "by_pattern": {"00": 45, "01/10": 46, "11": 44},
    "correct": 135,
    "accuracy": 1.0,
    "readout_time_s": 10.0,
    "figure": "Figure 9",
}
#: Two encodings named in the Results section (in-mixture concentrations).
PUBLISHED_ENCODINGS = [
    {"operation": "OR", "encoding_mM": {"Mg": 20.0}},
    {"operation": "XOR", "encoding_mM": {"Mg": 20.0, "Ca": 40.0}},
]


# --- the lock-and-key criterion ---------------------------------------------
def match(key: str, lock: str) -> int:
    """Complementary positions between a key and a lock (perfect = SHAPE_BITS)."""
    return sum(a != b for a, b in zip(key, lock))


def recognised(ligand: str, site: str, match_bits: int) -> bool:
    return match(KEYS[ligand], LOCKS[site]) >= match_bits


# --- the classification law (Table 1) ---------------------------------------
def milieu_state(inputs: str) -> str:
    """0, 1 or 2 one-signals in the cuvette give milieu state a, b or c."""
    return STATES[inputs.count("1")]


def truth_table(operation: str) -> dict[str, int]:
    """Output bit per input pattern, from the milieu states that must be active."""
    return {p: int(milieu_state(p) in ACTIVE_STATES[operation]) for p in PATTERNS}


def signal_strength(operation: str, response: dict[str, float]) -> float:
    """Delta s: the minimum response difference the enzyme has to supply."""
    active = ACTIVE_STATES[operation]
    return (min(response[s] for s in active)
            - max(response[s] for s in STATES if s not in active))


def implementable(operation: str, response: dict[str, float]) -> bool:
    """An operation is realisable with this response exactly when Delta s > 0."""
    return signal_strength(operation, response) > 0


def classify(operation: str, response: dict[str, float], threshold: float) -> dict[str, int]:
    """Threshold the response: above threshold is an active output."""
    return {p: int(response[milieu_state(p)] > threshold) for p in PATTERNS}


# --- the network -------------------------------------------------------------
def _carriers(encoding: object, match_bits: int) -> dict[str, float]:
    if not isinstance(encoding, dict) or not encoding:
        raise ValueError(
            f"encoding_mM must be a non-empty object mapping signal carriers to mM, got {encoding!r}"
        )
    unknown = [k for k in encoding if k not in IONS]
    if unknown:
        raise ValueError(
            f"encoding_mM: unknown signal carrier(s) {unknown}; "
            f"the published carriers are {list(IONS)}"
        )
    for ion, conc in encoding.items():
        if isinstance(conc, bool) or not isinstance(conc, (int, float)) or conc < 0:
            raise ValueError(f"encoding_mM[{ion!r}] must be a concentration in mM >= 0, got {conc!r}")
        score = match(KEYS[ion], LOCKS["modulator"])
        if score < match_bits:
            raise ValueError(
                f"encoding_mM: {ion} does not recognise the modulator site at "
                f"match_bits={match_bits} (its key is complementary in {score} of "
                f"{SHAPE_BITS} positions); lower match_bits or drop {ion}"
            )
    return {ion: float(encoding[ion]) for ion in IONS if ion in encoding}


def _conformer(bound: tuple[str, ...]) -> str:
    return "_".join(("MDH",) + bound)


def _conformer_structure(bound: tuple[str, ...]) -> str:
    sites = " ".join(f"{site}={lock}" for site, lock in LOCKS.items())
    if not bound:
        state = "modulator site free"
    else:
        state = "modulator site: " + ", ".join(
            f"{ion} key={KEYS[ion]} ({match(KEYS[ion], LOCKS['modulator'])}/{SHAPE_BITS})"
            for ion in bound
        )
    return f"malate dehydrogenase; locks {sites}; {state}"


def generate(p, rng):
    carriers = _carriers(p.encoding_mM, p.match_bits)
    if p.operation not in OPERATIONS:
        raise ValueError(f"operation must be one of {list(OPERATIONS)}, got {p.operation!r}")
    if p.inputs not in PATTERNS:
        raise ValueError(f"inputs must be one of {list(PATTERNS)}, got {p.inputs!r}")

    coenzymes = [l for l in LIGANDS if recognised(l, "coenzyme", p.match_bits)]
    substrates = [l for l in LIGANDS if recognised(l, "substrate", p.match_bits)]
    ions = [ion for ion in IONS if ion in carriers]
    conformers = [c for n in range(len(ions) + 1) for c in combinations(ions, n)]

    species: list[Species] = []
    reactions: list[tuple[str, None]] = []
    holds = {}          # species id -> set of ions bound, for the conservation laws
    carries_nad = []
    carries_c4 = []

    for bound in conformers:
        enzyme = _conformer(bound)
        species.append(Species(enzyme, _conformer_structure(bound)))
        holds[enzyme] = set(bound)
        for co in coenzymes:
            binary = f"{enzyme}_{co}"
            species.append(Species(binary, f"{_conformer_structure(bound)}; coenzyme site: "
                                           f"{co} key={KEYS[co]} "
                                           f"({match(KEYS[co], LOCKS['coenzyme'])}/{SHAPE_BITS})"))
            holds[binary] = set(bound)
            carries_nad.append(binary)
            reactions.append((f"{enzyme} + {co} -> {binary}", None))
            reactions.append((f"{binary} -> {enzyme} + {co}", None))
            for sub in substrates:
                ternary = f"{binary}_{sub}"
                species.append(Species(ternary, f"{_conformer_structure(bound)}; coenzyme site: "
                                                f"{co} key={KEYS[co]}; substrate site: {sub} "
                                                f"key={KEYS[sub]} "
                                                f"({match(KEYS[sub], LOCKS['substrate'])}/{SHAPE_BITS})"))
                holds[ternary] = set(bound)
                carries_nad.append(ternary)
                carries_c4.append(ternary)
                reactions.append((f"{binary} + {sub} -> {ternary}", None))
                reactions.append((f"{ternary} -> {binary} + {sub}", None))
                if (co, sub) == ("NAD", "MAL"):     # the assayed turnover
                    reactions.append((f"{ternary} -> {enzyme} + OAA + NADH + H", None))

    # Modulator binding: every edge of the lattice of bound-ion sets.
    for bound in conformers:
        for ion in ions:
            if ion in bound:
                continue
            after = _conformer(tuple(i for i in ions if i in set(bound) | {ion}))
            reactions.append((f"{_conformer(bound)} + {ion} -> {after}", None))
            reactions.append((f"{after} -> {_conformer(bound)} + {ion}", None))

    for ligand in ("NAD", "MAL", "OAA", "NADH"):
        species.append(Species(ligand, f"{NAMES[ligand]}; key={KEYS[ligand]}"))
    species.append(Species("H", NAMES["H"]))
    for ion in ions:
        species.append(Species(ion, f"{NAMES[ion]}; key={KEYS[ion]}"))

    enzymes = list(holds)
    conservation = [
        {"name": "enzyme", "vector": {s: 1 for s in enzymes}},
        {"name": "nicotinamide dinucleotide",
         "vector": {**{s: 1 for s in carries_nad}, "NAD": 1, "NADH": 1}},
        {"name": "C4 dicarboxylate",
         "vector": {**{s: 1 for s in carries_c4}, "MAL": 1, "OAA": 1}},
    ]
    for ion in ions:
        conservation.append({
            "name": f"{ion} (free and bound)",
            "vector": {**{s: 1 for s, b in holds.items() if ion in b}, ion: 1},
        })

    n_ones = p.inputs.count("1")
    initial = {"NAD": float(p.nad_mM), "MAL": float(p.malate_mM)}
    initial.update({ion: n_ones * carriers[ion] for ion in ions})

    extras = {
        "interaction_law": {
            "recognition": "lock and key: a ligand binds a site when at least match_bits "
                           "of the shape_bits positions are complementary (bit-inverted)",
            "shape_bits": SHAPE_BITS,
            "match_bits": int(p.match_bits),
            "keys": dict(KEYS),
            "locks": dict(LOCKS),
            "match_scores": {f"{site}/{lig}": match(KEYS[lig], LOCKS[site])
                             for site in LOCKS for lig in LIGANDS},
            "conformational_fusion": "Mg2+ and Ca2+ bind the modulator site and select a "
                                     "conformer; the conformational dynamics fuses the milieu "
                                     "signals nonlinearly, so the catalytic activity is a "
                                     "nonmonotonic function of each signal concentration",
            "response": "r = NADH absorbance at 339 nm at a fixed time after initiation; "
                        "measured to rise to a maximum and then fall along both ion axes",
            "signal_coding": "a 1-signal is the presence of one encoding unit of the carrier "
                             "ions, a 0-signal their absence; two input lines give milieu "
                             "states a (00), b (01 or 10) and c (11)",
            "signal_strength": "Delta s = min(r over the milieu states that must be active) "
                               "- max(r over the rest); the operation is realisable iff Delta s > 0",
            "output": "threshold the response; above threshold is an active output",
        },
        "analysis": {
            "operation": p.operation,
            "truth_table": truth_table(p.operation),
            "active_milieu_states": list(ACTIVE_STATES[p.operation]),
            "signal_strength": SIGNAL_STRENGTH[p.operation],
            "inputs": p.inputs,
            "milieu_state": milieu_state(p.inputs),
            "encoding_mM": dict(carriers),
            "milieu_mM": {ion: {s: k * carriers[ion] for k, s in enumerate(STATES)}
                          for ion in ions},
            "linearly_separable": p.operation not in ("XOR", "NXOR"),
            "implementable_on_the_published_response_surface": p.operation != "NXOR",
            "not_implementable": ["NXOR"],
            "published_classification": dict(PUBLISHED_CLASSIFICATION),
            "published_encodings": [dict(e) for e in PUBLISHED_ENCODINGS],
        },
        "conservation": conservation,
        "buffered": ["H"],
        "experiment": EXPERIMENT,
    }
    return network(reactions, species=species, initial_state=initial, extras=extras)
