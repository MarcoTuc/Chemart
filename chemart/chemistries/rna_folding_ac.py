"""RNA-folding ribozyme artificial chemistry (Ullrich & Flamm). Catalog id: rna-folding-ac.

Book section 18.1.1. Molecules are RNA sequences. Every sequence is folded into
its minimum-free-energy secondary structure with the Vienna RNA package (the
genotype -> phenotype map), and the *longest loop* of that structure is read as
the imaginary transition structure (ITS) of the reaction the molecule catalyses
(the phenotype -> function map of Flamm et al. 2010):

- level 1, the size of the ITS, is the length of the loop cycle: the unpaired
  bases of the loop plus the two bases of each stem that delimits it;
- level 3, the arrangement of the electron re-ordering, comes from the stems
  delimiting the loop;
- level 4, the atom types, comes from the sequence inside the loop; here it is
  the recognition word the ribozyme binds.

A mono-cyclic ITS of size 2n breaks n bonds and forms n bonds alternately
around the cycle (Fujita 1986), so only loops of even size are catalytic.
Hairpin loops (one stem) give cleavage, interior and bulge loops (two stems)
give ligation, and multiloops (composite ITS) are not catalytic:

    E + X      -> E + Y + Z     cleavage of X = YZ at the recognition site
    E + Y + Z  -> E + X         ligation of a junction that reads the site

The ribozyme binds the Watson-Crick reverse complement of its loop sequence, so
both reactions conserve the number of nucleotides (extras["conservation"]).

Two faces. `generate` returns every reaction reachable from the initial pool
(the closure), with one substrate for a cleavage and two for a ligation.
`evolve` runs a well-stirred multiset (chemart.soup.stir), a frame every
`pool` collisions, and returns the reactions that fired, each with its firing
count; a collision there always draws three molecules, so a cleavage recorded
in a run carries the third molecule through unchanged as a spectator.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Callable

import RNA

from chemart.expand import expand
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame

ALPHABET = "ACGU"
COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


def fold(seq: str) -> tuple[str, float]:
    """MFE secondary structure in dot-bracket notation and its free energy."""
    structure, energy = RNA.fold(seq)
    return structure, float(energy)


def loops(structure: str) -> list[tuple[int, int, tuple[int, ...]]]:
    """Every loop of the secondary-structure graph.

    Returns (cycle size, number of delimiting stems, unpaired positions) per
    loop, where the cycle size counts the unpaired bases plus the two bases of
    each delimiting pair. The exterior loop is not a cycle and is excluded.
    """
    pairs = list(RNA.ptable(structure))
    n = pairs[0]
    out: list[tuple[int, int, tuple[int, ...]]] = []
    for i in range(1, n + 1):
        j = pairs[i]
        if j <= i:
            continue
        unpaired: list[int] = []
        stems = 1
        k = i + 1
        while k < j:
            if pairs[k] > k:          # an inner stem delimits this loop
                stems += 1
                k = pairs[k] + 1
            else:
                unpaired.append(k)
                k += 1
        out.append((len(unpaired) + 2 * stems, stems, tuple(unpaired)))
    return out


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT[c] for c in reversed(seq))


def function(seq: str, structure: str, its_min: int, its_max: int,
             recognition: int) -> dict[str, Any] | None:
    """The catalytic function of a folded sequence, or None if it is inert."""
    cycles = loops(structure)
    if not cycles:
        return None
    size, stems, unpaired = max(cycles)
    if size % 2 or not (its_min <= size <= its_max) or stems > 2:
        return None
    motif = "".join(seq[i - 1] for i in unpaired)
    if len(motif) < recognition:
        return None
    return {
        "its_size": int(size),
        "stems": int(stems),
        "loop_sequence": motif,
        "site": revcomp(motif[:recognition]),
        "reaction": "cleavage" if stems == 1 else "ligation",
    }


def cleave(substrate: str, site: str, min_fragment: int) -> tuple[str, str] | None:
    """Cut immediately 3' of the recognition site, or None if it is absent."""
    pos = substrate.find(site)
    if pos < 0:
        return None
    cut = pos + len(site)
    if cut < min_fragment or len(substrate) - cut < min_fragment:
        return None
    return substrate[:cut], substrate[cut:]


def joins(left: str, right: str, site: str) -> bool:
    """True if the junction left|right spells the recognition site."""
    return any(left.endswith(site[:a]) and right.startswith(site[a:])
               for a in range(1, len(site)))


def reaction_rule(function_of: Callable[[str], dict | None], min_fragment: int):
    """Closure rule: one substrate is cleaved, two substrates are ligated.

    The first molecule is the catalyst. Keeping cleavage bimolecular and
    ligation termolecular leaves no spectator species in the closure.
    """

    def react(*molecules: str):
        catalyst, substrates = molecules[0], molecules[1:]
        f = function_of(catalyst)
        if f is None:
            return None
        site, kind = f["site"], f["reaction"]
        out = []
        if kind == "cleavage" and len(substrates) == 1:
            pieces = cleave(substrates[0], site, min_fragment)
            if pieces:
                out.append((catalyst, *pieces))
        elif kind == "ligation" and len(substrates) == 2:
            for left, right in (substrates, substrates[::-1]):
                if joins(left, right, site):
                    out.append((catalyst, left + right))
        return out or None

    return react


def collision_rule(function_of: Callable[[str], dict | None], min_fragment: int, rng):
    """Well-stirred rule: three molecules meet and at most one reaction happens.

    A cleavase cuts one of the two substrates and the other is carried through
    unchanged, so nothing is lost from the population.
    """

    def collide(*molecules: str):
        catalyst, substrates = molecules[0], molecules[1:]
        f = function_of(catalyst)
        if f is None:
            return None
        site, kind = f["site"], f["reaction"]
        out = []
        if kind == "cleavage":
            for k, substrate in enumerate(substrates):
                pieces = cleave(substrate, site, min_fragment)
                if pieces:
                    spectators = substrates[:k] + substrates[k + 1:]
                    out.append((catalyst, *pieces, *spectators))
        else:
            for left, right in (substrates, substrates[::-1]):
                if joins(left, right, site):
                    out.append((catalyst, left + right))
        if not out:
            return None
        return out[int(rng.integers(len(out)))]

    return collide


def _check(p) -> None:
    if p.its_min > p.its_max:
        raise ValueError(f"its_min ({p.its_min}) cannot exceed its_max ({p.its_max})")
    if p.min_recognition > p.seq_length:
        raise ValueError(
            f"min_recognition ({p.min_recognition}) cannot exceed seq_length ({p.seq_length})"
        )


def _maps(p):
    """Cached sequence -> (structure, energy) and sequence -> function maps."""
    folded: dict[str, tuple[str, float]] = {}
    functions: dict[str, dict | None] = {}

    def fold_of(seq: str) -> tuple[str, float]:
        if seq not in folded:
            folded[seq] = fold(seq)
        return folded[seq]

    def function_of(seq: str) -> dict | None:
        if seq not in functions:
            structure, _ = fold_of(seq)
            functions[seq] = function(
                seq, structure, p.its_min, p.its_max, p.min_recognition
            )
        return functions[seq]

    return fold_of, function_of


def _seeds(p, rng) -> list[str]:
    return list(dict.fromkeys(
        "".join(rng.choice(list(ALPHABET), p.seq_length)) for _ in range(p.pool)
    ))


def generate(p, rng):
    """The closure of the initial pool, cut off by max_species."""
    _check(p)
    fold_of, function_of = _maps(p)
    seeds = _seeds(p, rng)
    react = reaction_rule(function_of, p.min_fragment)
    found, pairs, status = expand(
        react, seeds, arity=(2, 3), max_species=p.max_species, alternatives=True
    )
    reactions = [Reaction.of(list(lhs), list(rhs)) for lhs, rhs in pairs]
    return _network(found, reactions, status, seeds, fold_of, function_of)


def evolve(p, rng):
    """A well-stirred run of three-molecule collisions: a frame every `pool` collisions."""
    _check(p)
    if p.pool < 3:
        raise ValueError(
            "a well-stirred collision draws 3 molecules, so pool must be >= 3 "
            f"to evolve, got {p.pool}"
        )
    fold_of, function_of = _maps(p)
    seeds = _seeds(p, rng)
    collide = collision_rule(function_of, p.min_fragment, rng)
    tally = Tally()
    for step, pop, tally in stir(collide, seeds, p.steps, rng, arity=3, dilution=p.dilution, tally=tally):
        yield Frame(t=float(step), state={s: float(n) for s, n in Counter(pop).items()},
                    fired=tally.flush())
    events = tally.reactions()
    found = list(dict.fromkeys(
        [*seeds, *(m for lhs, rhs, _ in events for m in (*lhs, *rhs))]
    ))
    reactions = [Reaction.of(list(lhs), list(rhs), count=int(n)) for lhs, rhs, n in events]
    return _network(found, reactions, "observed", seeds, fold_of, function_of)


def _network(found, reactions, status, seeds, fold_of, function_of) -> Network:
    species, energies, catalysts = [], {}, {}
    for seq in found:
        structure, energy = fold_of(seq)
        species.append(Species(seq, f"{seq} {structure}"))
        energies[seq] = round(energy, 2)
        f = function_of(seq)
        if f is not None:
            catalysts[seq] = f

    return Network(
        species=species,
        reactions=reactions,
        status=status,
        initial_state={s: float(n) for s, n in Counter(seeds).items()},
        extras={
            "energies": energies,
            "conservation": [{
                "name": "nucleotides",
                "vector": [len(s.id) for s in species],
            }],
            "functions": catalysts,
            "analysis": {
                "catalytic_species": len(catalysts),
                "cleavases": sum(f["reaction"] == "cleavage" for f in catalysts.values()),
                "ligases": sum(f["reaction"] == "ligation" for f in catalysts.values()),
                "sequence_to_function": (
                    "MFE fold (Vienna RNA); the longest loop cycle is the ITS; "
                    "even size in [its_min, its_max] with one stem (hairpin) gives "
                    "cleavage, two stems (interior/bulge) give ligation; the "
                    "recognition site is the reverse complement of the loop sequence"
                ),
            },
        },
    )
