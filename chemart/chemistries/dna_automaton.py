"""DNA automaton (Benenson & Shapiro): a finite automaton made of DNA and FokI.

Catalog id: dna-automaton. Book section 19.3.1 ("DNA Automaton"), whose detailed
source is Benenson, Adar, Paz-Elizur, Livneh & Shapiro, PNAS 100:2191 (2003) -
the book's reference [100] and the ATP-free successor of the original automaton
(Nature 414:430, 2001). This module returns the in-silico reaction network of
that published experiment, with the published sequences.

Molecules (PNAS 2003, fig. 2):

- **hardware**: the type IIS restriction enzyme FokI, which binds ``GGATG`` and
  cleaves 9 nt downstream on the sense strand and 13 nt downstream on the
  antisense strand, leaving a 4-nt 5' sticky end;
- **software**: eight transition molecules T1-T8, each a short duplex carrying a
  4-nt 5' antisense overhang that detects one <state, symbol> pair, the FokI
  site, and a spacer of 0, 1 or 2 bp that sets where FokI will cut the input and
  hence the next state (0 bp: S1 -> S0, 1 bp: keep the state, 2 bp: S0 -> S1);
- **input**: a duplex in which every symbol is 5 bp (a = ``TGGCT``,
  b = ``GCAGG``, terminator t = ``GTCGG``) separated by the 3-bp spacer
  ``GCC``. The exposed 4-nt sticky end reads <state, symbol>: the leftmost four
  bases of a symbol mean S1, the rightmost four mean S0. The input is delivered
  already cleaved, so its sticky end is <S0, first symbol>.

One transition (fig. 2E) is three reactions::

    FokI + Ti            <-> FokI-Ti                    (Kd = 2 nM)
    FokI-Ti + X          <-> FokI-Ti-X                  (sticky-end hybridisation, Kd = 50 uM)
    FokI-Ti-X            ->  FokI-Ti + fragment + X'    (FokI cleavage, irreversible)

The software and the enzyme are recycled; the cleaved symbol (a 7-9 nt duplex
with only 3-5 bp paired) is scattered, and the hydrolysis of its two
phosphodiester bonds is what fuels the computation. The cascade stops when the
sticky end encodes the terminator, which no transition molecule detects: the
remaining molecule is the output, whose antisense strand is 15 nt for S0 and
16 nt for S1 - the two bands of the published gels.

The transition rule of every software molecule here is *derived from its
published oligonucleotides* (sticky end -> detected <state, symbol>, spacer
length -> next state) and reproduces fig. 2C exactly; the tests check this.
"""

from __future__ import annotations

from typing import Any

from chemart.network import Network, Reaction, Species

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

#: FokI: recognition site and the two cut offsets (PNAS 2003, fig. 2B).
FOKI_SITE = "GGATG"
FOKI_SENSE_CUT = 9
FOKI_ANTISENSE_CUT = 13
STICKY = FOKI_ANTISENSE_CUT - FOKI_SENSE_CUT      # 4-nt 5' overhang

#: Symbol encoding (PNAS 2003, fig. 2A): 5 bp per symbol, 3-bp spacers.
SYMBOL_DUPLEX = {"a": "TGGCT", "b": "GCAGG", "t": "GTCGG"}
SYMBOL_SPACER = "GCC"
TERMINATOR_TAIL = "TACCGATTAAGTTGGA"   # the last base is the 1-nt 3' extrusion
ALPHABET = ("a", "b")
TERMINATOR = "t"
STATES = ("S0", "S1")

#: Published software oligonucleotides (PNAS 2003, Materials and Methods).
SOFTWARE_SENSE_OLIGOS = {
    "TN1368": "AAGAGCTAGAGTCGGATGC",
    "TN24": "AAGAGCTAGAGTCGGATGCC",
    "TN57": "AAGAGCTAGAGTCGGATG",
}
SOFTWARE_OLIGOS = {
    "T1": ("TN1368", "AGCCGCATCCGACTCTAGCTCT"),
    "T2": ("TN24", "AGCCGGCATCCGACTCTAGCTCT"),
    "T3": ("TN1368", "CCTGGCATCCGACTCTAGCTCT"),
    "T4": ("TN24", "CCTGGGCATCCGACTCTAGCTCT"),
    "T5": ("TN57", "GCCACATCCGACTCTAGCTCT"),
    "T6": ("TN1368", "GCCAGCATCCGACTCTAGCTCT"),
    "T7": ("TN57", "CTGCCATCCGACTCTAGCTCT"),
    "T8": ("TN1368", "CTGCGCATCCGACTCTAGCTCT"),
}

#: spacer length of the software molecule -> next state, given the detected state.
NEXT_STATE = {("S0", 1): "S0", ("S0", 2): "S1", ("S1", 0): "S0", ("S1", 1): "S1"}

#: The published programs (PNAS 2003, fig. 1A and 1C) and the book's example
#: automaton (Banzhaf & Yamamoto, fig. 19.13), which was not run in the lab.
PROGRAMS = {
    "A1": {"transitions": ["T2", "T3", "T5", "T8"], "accepting": "S0",
           "decides": "the input contains an even number of a symbols"},
    "A2": {"transitions": ["T2", "T4", "T5", "T7"], "accepting": "S0",
           "decides": "the input contains an even number of symbols"},
    "A3": {"transitions": ["T1", "T4", "T5", "T8"], "accepting": "S1",
           "decides": "the input ends with b"},
    "book": {"transitions": ["T1", "T4", "T6", "T7"], "accepting": "S0",
             "decides": "the input contains an even number of b symbols "
                        "(book fig. 19.13; not among the automata run in the lab)"},
}

#: The published inputs (PNAS 2003, fig. 1D).
PUBLISHED_INPUTS = {
    "I1": "abb", "I2": "abba", "I3": "babbabb", "I4": "babbabba",
    "I5": "baaaabb", "I6": "baaaabba", "I7": "abbbbabbabb", "I8": "abbbbaaaabba",
}

#: Published reaction conditions. `seconds_per_transition` is only given where
#: the paper measured it; concentrations are micromolar.
CONDITIONS = {
    "standard": {
        "experiment": "fig. 3B: automata A1-A3 on inputs I1-I8, 10 ul, 8 C, 20 min",
        "input_uM": 1.0, "software_each_uM": 1.0, "hardware_uM": 4.0,
        "temperature_C": 8.0, "seconds_per_transition": None,
    },
    "fast": {
        "experiment": "fastest computation: program A2 (1 uM each) and 4 uM FokI "
                      "with 10 nM input I3 at 8 C, 20 s per transition",
        "input_uM": 0.01, "software_each_uM": 1.0, "hardware_uM": 4.0,
        "temperature_C": 8.0, "seconds_per_transition": 20.0,
    },
    "parallel": {
        "experiment": "maximal parallelism: 10 uM program A2 (2.5 uM each) and "
                      "10 uM FokI with 5 uM input I3, 45 s per transition, "
                      "6.646e10 transitions per second per ul",
        "input_uM": 5.0, "software_each_uM": 2.5, "hardware_uM": 10.0,
        "temperature_C": 8.0, "seconds_per_transition": 45.0,
    },
    "reuse": {
        "experiment": "fig. 3C: software recycling, A1 on I8 for 18 h with each "
                      "software molecule at 0.075 molar ratio to the input",
        "input_uM": 1.0, "software_each_uM": 0.075, "hardware_uM": 1.0,
        "temperature_C": 8.0, "seconds_per_transition": None,
    },
}

#: Numbers the paper reports, carried through to extras["published"].
PUBLISHED = {
    "single_step_fidelity": 0.999,
    "average_fidelity_12_symbol_inputs": 0.995,
    "fidelity_A1_on_I2": 0.95,
    "automata_per_ul": 3.0e12,
    "transitions_per_second_per_ul": 6.646e10,
    "seconds_per_transition_fastest": 20.0,
    "seconds_per_transition_parallel": 45.0,
    "seconds_per_transition_2001_automaton": 1000.0,
    "heat_dissipation_W_per_ul": 5.3e-9,
    "free_energy_per_transition_kT": -33.9,
    "dG_two_phosphodiester_bonds_kJ_per_mol_298K": -44.31,
    "dG_two_phosphodiester_bonds_kJ_per_mol_281K": -41.78,
    "Kd_hardware_software_nM": 2.0,
    "Kd_software_input_uM": 50.0,
    "transitions_per_software_molecule_fig3C": {"T2": 29, "T5": 21, "T8": 54},
    "output_antisense_nt": {"S0": 15, "S1": 16},
}


# --- sequences ---------------------------------------------------------------
def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT[base] for base in reversed(seq))


def sticky_end(state: str, symbol: str) -> str:
    """The 4-nt 5' sense overhang that encodes <state, symbol> (fig. 2A)."""
    duplex = SYMBOL_DUPLEX[symbol]
    return duplex[:STICKY] if state == "S1" else duplex[1:1 + STICKY]


def software_duplex(name: str) -> tuple[str, str]:
    """The (sense, antisense) oligonucleotide pair of a transition molecule."""
    sense_oligo, antisense = SOFTWARE_OLIGOS[name]
    return SOFTWARE_SENSE_OLIGOS[sense_oligo], antisense


def _rule(name: str) -> dict[str, Any]:
    """Read one transition rule off the published oligonucleotides."""
    sense, antisense = software_duplex(name)
    spacer = len(sense) - (sense.index(FOKI_SITE) + len(FOKI_SITE))
    end = antisense[:STICKY]
    detected = [(state, symbol) for state in STATES for symbol in ALPHABET
                if revcomp(end) == sticky_end(state, symbol)]
    if len(detected) != 1 or (detected[0][0], spacer) not in NEXT_STATE:
        raise ValueError(f"software molecule {name} does not encode a legal transition")
    state, symbol = detected[0]
    return {
        "state": state, "symbol": symbol, "next": NEXT_STATE[(state, spacer)],
        "spacer": spacer, "sticky_end": end, "sense": sense, "antisense": antisense,
    }


#: T1-T8 with the transition each one implements (reproduces fig. 2C).
TRANSITIONS = {name: _rule(name) for name in SOFTWARE_OLIGOS}


def input_molecule(state: str, remaining: str) -> tuple[str, str]:
    """The (sense, antisense) duplex exposing <state, first symbol of `remaining`>.

    `remaining` is the not-yet-processed symbol string; the terminator is always
    appended, so an empty `remaining` gives the output molecule.
    """
    symbols = remaining + TERMINATOR
    first = SYMBOL_DUPLEX[symbols[0]]
    sense = first if state == "S1" else first[1:]
    for symbol in symbols[1:]:
        sense += SYMBOL_SPACER + SYMBOL_DUPLEX[symbol]
    sense += TERMINATOR_TAIL
    return sense, revcomp(sense[STICKY:-1])


def cleave(sense: str, spacer: int) -> tuple[str, str, str]:
    """FokI cleavage of a software/input hybrid whose software has `spacer` bp.

    Counting from the recognition site, the software contributes `spacer` bases
    before the input, so FokI cuts the input's sense strand after 9 - spacer
    bases and its antisense strand after 13 - spacer. Returns the scattered
    symbol (sense, antisense) and the sense strand of the remaining input.
    """
    sense_cut = FOKI_SENSE_CUT - spacer
    antisense_cut = FOKI_ANTISENSE_CUT - spacer
    if len(sense) < antisense_cut:
        raise ValueError(f"molecule of {len(sense)} nt is too short for a FokI transition")
    return sense[:sense_cut], revcomp(sense[STICKY:antisense_cut]), sense[sense_cut:]


def molecule_id(state: str, remaining: str) -> str:
    return f"{state}_{remaining}{TERMINATOR}"


def run(rules: dict[str, dict], symbols: str) -> dict[str, Any]:
    """Run a program (name -> rule) on a symbol string; returns the trace."""
    by_sticky = {rule["sticky_end"]: name for name, rule in rules.items()}
    state, remaining, steps = "S0", symbols, []
    while True:
        sense, _ = input_molecule(state, remaining)
        name = by_sticky.get(revcomp(sense[:STICKY]))
        if name is None:
            break
        rule = rules[name]
        fragment_sense, fragment_antisense, rest = cleave(sense, rule["spacer"])
        steps.append({
            "transition": name, "from": state, "symbol": remaining[0],
            "to": rule["next"], "molecule": molecule_id(state, remaining),
            "fragment": (fragment_sense, fragment_antisense),
            "product": molecule_id(rule["next"], remaining[1:]),
        })
        state, remaining = rule["next"], remaining[1:]
        if rest != input_molecule(state, remaining)[0]:   # pragma: no cover
            raise ValueError("cleavage did not produce the expected molecule")
    output = input_molecule(state, remaining)
    return {
        "steps": steps, "final_state": state, "suspended": remaining != "",
        "unprocessed": remaining, "output": molecule_id(state, remaining),
        "output_antisense_nt": len(output[1]),
    }


# --- parameters --------------------------------------------------------------
def _program(p) -> tuple[dict[str, dict], str, str]:
    """Resolve the program to (rules by name, accepting state, description)."""
    if p.transitions:
        names = []
        for raw in p.transitions:
            if raw not in TRANSITIONS:
                raise ValueError(
                    f"unknown transition molecule {raw!r}; the published software "
                    f"molecules are {', '.join(TRANSITIONS)}"
                )
            if raw not in names:
                names.append(raw)
        accepting, description = p.accepting_state, "custom program"
    else:
        program = PROGRAMS[p.program]
        names = list(program["transitions"])
        accepting, description = program["accepting"], program["decides"]
    seen: dict[tuple[str, str], str] = {}
    for name in names:
        rule = TRANSITIONS[name]
        key = (rule["state"], rule["symbol"])
        if key in seen:
            raise ValueError(
                f"{seen[key]} and {name} both detect <{key[0]}, {key[1]}>: the automaton "
                "would be nondeterministic, which this deterministic model does not cover"
            )
        seen[key] = name
    return {name: TRANSITIONS[name] for name in names}, accepting, description


def _inputs(values) -> list[tuple[str, str]]:
    """Resolve input names or literal symbol strings to (name, symbols) pairs."""
    if not values:
        raise ValueError("inputs must list at least one input molecule")
    out: list[tuple[str, str]] = []
    for raw in values:
        if not isinstance(raw, str):
            raise ValueError(f"inputs entries must be strings, got {raw!r}")
        if raw in PUBLISHED_INPUTS:
            item = (raw, PUBLISHED_INPUTS[raw])
        elif raw and set(raw) <= set(ALPHABET):
            item = (raw, raw)
        else:
            raise ValueError(
                f"input {raw!r} is neither a published input "
                f"({', '.join(PUBLISHED_INPUTS)}) nor a non-empty string over "
                f"{{{', '.join(ALPHABET)}}}"
            )
        if item not in out:
            out.append(item)
    return out


# --- the network -------------------------------------------------------------
class _Builder:
    def __init__(self) -> None:
        self.species: dict[str, Species] = {}
        self.nucleotides: dict[str, int] = {}
        self.reactions: dict[tuple, Reaction] = {}

    def add(self, sid: str, structure: str | None, nucleotides: int) -> str:
        if sid not in self.species:
            self.species[sid] = Species(sid, structure)
            self.nucleotides[sid] = nucleotides
        return sid

    def react(self, left: list[str], right: list[str], rate: dict | None = None) -> None:
        key = (tuple(sorted(left)), tuple(sorted(right)))
        if key not in self.reactions:
            self.reactions[key] = Reaction.of(left, right, rate=rate)


def _duplex(sense: str, antisense: str) -> str:
    return f"{sense}/{antisense}"


def generate(p, rng):
    rules, accepting, description = _program(p)
    inputs = _inputs(p.inputs)
    conditions = CONDITIONS[p.conditions]
    seconds = conditions["seconds_per_transition"]
    cleavage_rate = None if seconds is None else {
        "law": "mass-action", "k": round(1.0 / seconds, 6), "units": "per second",
        "conditions": conditions["experiment"],
    }

    build = _Builder()
    hardware = build.add(
        "FokI",
        "protein: FokI type IIS endonuclease; binds GGATG and cleaves 9/13 nt downstream",
        0,
    )
    complexes = {}
    for name, rule in rules.items():
        sense, antisense = rule["sense"], rule["antisense"]
        software = build.add(name, _duplex(sense, antisense), len(sense) + len(antisense))
        loaded = build.add(f"FokI-{name}", f"FokI + {_duplex(sense, antisense)}",
                           len(sense) + len(antisense))
        complexes[name] = loaded
        build.react([hardware, software], [loaded])
        build.react([loaded], [hardware, software])

    initial: dict[str, float] = {hardware: conditions["hardware_uM"]}
    for name in rules:
        initial[name] = conditions["software_each_uM"]

    runs = []
    for input_name, symbols in inputs:
        trace = run(rules, symbols)
        start = molecule_id("S0", symbols)
        state, remaining = "S0", symbols
        for step in trace["steps"]:
            sense, antisense = input_molecule(state, remaining)
            molecule = build.add(step["molecule"], _duplex(sense, antisense),
                                 len(sense) + len(antisense))
            loaded = complexes[step["transition"]]
            hybrid = build.add(
                f"{loaded}-{step['molecule']}",
                f"{build.species[loaded].structure} + {_duplex(sense, antisense)}",
                build.nucleotides[loaded] + len(sense) + len(antisense),
            )
            fragment_sense, fragment_antisense = step["fragment"]
            fragment = build.add(
                f"frag_{fragment_sense}_{fragment_antisense}",
                _duplex(fragment_sense, fragment_antisense),
                len(fragment_sense) + len(fragment_antisense),
            )
            state, remaining = step["to"], remaining[1:]
            product_sense, product_antisense = input_molecule(state, remaining)
            product = build.add(step["product"], _duplex(product_sense, product_antisense),
                                len(product_sense) + len(product_antisense))
            build.react([loaded, molecule], [hybrid])
            build.react([hybrid], [loaded, molecule])
            build.react([hybrid], [loaded, fragment, product], rate=cleavage_rate)
        output_sense, output_antisense = input_molecule(trace["final_state"], trace["unprocessed"])
        build.add(trace["output"], _duplex(output_sense, output_antisense),
                  len(output_sense) + len(output_antisense))
        initial[start] = initial.get(start, 0.0) + conditions["input_uM"]
        runs.append({
            "input": input_name,
            "symbols": symbols,
            "molecule": start,
            "final_state": trace["final_state"],
            "accepted": (not trace["suspended"]) and trace["final_state"] == accepting,
            "suspended": trace["suspended"],
            "unprocessed": trace["unprocessed"],
            "output": trace["output"],
            "output_antisense_nt": trace["output_antisense_nt"],
            "transitions": [
                {k: step[k] for k in ("transition", "from", "symbol", "to", "molecule", "product")}
                for step in trace["steps"]
            ],
        })

    species = list(build.species.values())
    return Network(
        species=species,
        reactions=list(build.reactions.values()),
        status="complete",
        initial_state=initial,
        extras={
            "conservation": [{
                "name": "nucleotides",
                "vector": [build.nucleotides[s.id] for s in species],
            }],
            "conditions": {"name": p.conditions, **conditions},
            "analysis": {
                "program": "custom" if p.transitions else p.program,
                "decides": description,
                "accepting_state": accepting,
                "transition_table": {
                    name: f"{rule['state']} -{rule['symbol']}-> {rule['next']}"
                    for name, rule in rules.items()
                },
                "runs": runs,
                "transitions_per_software_molecule": {
                    name: sum(step["transition"] == name
                              for r in runs for step in r["transitions"])
                    for name in rules
                },
            },
            "published": PUBLISHED,
        },
    )
