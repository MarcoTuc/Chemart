"""HBCB / PSD: the hierarchical biomolecular covalent bond model and
programmed self-decomposition. Catalog id: hbcb-psd. Book 18.3.2, ref [643]
(Oohashi, Ueno, Maekawa, Kawai, Nishina & Honda, Artificial Life 15(1), 2009).

The HBCB model assumes that the biomolecules of terrestrial life are arranged
in a hierarchy of classes ordered by the energy of the covalent bonds that
hold each class together. The bonds that link biomonomers (BM) into a
biopolymer (BP) are the weakest of the hierarchy, and hydrolysis -- the main
process of the programmed self-decomposition (PSD) mechanism -- is exactly the
cleavage of those bonds; the bonds inside a monomer, and inside every class
below it, are stronger.

The PSD model assumes that an organism decomposes its own biomolecules by a
genetically regulated, endergonic process, and that evolution selected *how
deep* that decomposition goes: cleaving a BP back to BMs leaves matter that
can be re-used directly, while cleaving it into classes below BM leaves matter
that has to be rebuilt through every intermediate class. The paper's
simulation compares virtual species that differ only in that depth, and the
BP -> BM species wins.

One network is the molecular hierarchy of one closed ecosystem. For each
family of biomolecules it holds the classes BP, BM, and the lower classes
C(n-3) ... C0 down to the base class; the hydrolysis (cleavage) and synthesis
(re-use) reactions that connect neighbouring classes, derived as the closure
of the cleavage rule from the biopolymers (chemart.expand); and the one
programmed self-decomposition reaction per family, which turns a BP directly
into the class named by `psd_target`. Matter is finite and recycled, so there
is no flow, and the base units of every family are conserved exactly
(extras["conservation"]).

What is *not* modelled: the ecosystem level. The paper's virtual ecosystem is
a population of competing organisms on a lattice, with mutation, a temperature
gradient and finite matter (see the SIVA lineage, Sayama 1998); selection
between decomposition strategies acts across organisms, not inside one
reaction network. The cost of returning to a BP from each class is reported in
extras["analysis"] instead. No rate constants are published, so every rate is
None.
"""

from __future__ import annotations

from chemart.expand import expand
from chemart.helpers.params import vector
from chemart.network import Network, Reaction, Species

#: Depth of the programmed self-decomposition, as a class of the hierarchy.
PSD_TARGETS = ("none", "BM", "sub-BM", "base")


def class_names(n_levels: int) -> list[str]:
    """Class of every level, from the base class (0) up to the biopolymer.

    n_levels = 4 gives ["C0", "C1", "BM", "BP"]: BP is the top class and BM
    the class directly below it, the two the model names.
    """
    names = []
    for level in range(n_levels):
        if level == n_levels - 1:
            names.append("BP")
        elif level == n_levels - 2:
            names.append("BM")
        else:
            names.append(f"C{level}")
    return names


def species_id(name: str, family: int) -> str:
    """Species id of a class in one family of biomolecules: 'BM_0'."""
    return f"{name}_{family}"


def molecular_energy(level: int, units: int, energies: list[float]) -> float:
    """Total covalent bond energy inside one molecule of `level`.

    A class-l molecule is `units` class-(l-1) molecules joined by units - 1
    covalent bonds of energy energies[l - 1], so the energies accumulate down
    the hierarchy.
    """
    total = 0.0
    for l in range(1, level + 1):
        total = units * total + (units - 1) * energies[l - 1]
    return total


def rebuild_cost(level: int, n_levels: int, units: int,
                 energies: list[float]) -> tuple[int, float]:
    """(syntheses, bond energy) needed to rebuild one BP from class `level`.

    The re-use cost of the PSD model: decomposing to BM leaves one synthesis
    to do, decomposing further multiplies both the number of synthesis steps
    and the covalent bond energy that has to be put back in.
    """
    syntheses, energy = 0, 0.0
    for l in range(level + 1, n_levels):
        count = units ** (n_levels - 1 - l)
        syntheses += count
        energy += count * (units - 1) * energies[l - 1]
    return syntheses, energy


def _bond_energies(p) -> list[float]:
    values = vector("bond_energies", p.bond_energies, p.n_levels - 1)
    if any(e <= 0 for e in values):
        raise ValueError(f"bond_energies entries must be positive, got {p.bond_energies!r}")
    if any(a <= b for a, b in zip(values, values[1:])):
        raise ValueError(
            "bond_energies must decrease strictly from the base class upwards: the HBCB "
            "model orders the hierarchy by bond strength, so the bonds inside the lower "
            f"classes are the stronger ones and the BM-BM links of a BP the weakest, got {values!r}"
        )
    return values


def _psd_level(target: str, n_levels: int) -> int | None:
    """Level that programmed self-decomposition cleaves a BP down to."""
    if target == "none":
        return None
    if target == "BM":
        return n_levels - 2
    if target == "base":
        return 0
    level = n_levels - 3
    if level < 0:
        raise ValueError(
            f"psd_target='sub-BM' needs n_levels >= 3: with n_levels={n_levels} the only "
            "classes are BP and BM, so there is no class below BM to decompose into"
        )
    return level


def generate(p, rng):
    n, units = p.n_levels, p.units_per_class
    energies = _bond_energies(p)
    target = _psd_level(p.psd_target, n)
    names = class_names(n)
    families = list(range(p.n_families))

    def cleave(mol):
        """Hydrolysis: one class-l molecule falls apart into its subunits."""
        family, level = mol
        if level == 0:
            return None
        return tuple((family, level - 1) for _ in range(units))

    # The hierarchy is the closure of the cleavage rule from the biopolymers.
    found, cleavages, status = expand(
        cleave, [(f, n - 1) for f in families], arity=1, ordered=False,
        max_species=p.n_families * n,
    )
    sid = {(f, level): species_id(names[level], f) for f in families for level in range(n)}
    if set(found) != set(sid):
        raise ValueError(f"the cleavage closure reached {len(found)} classes, expected {len(sid)}")

    species = [
        Species(
            sid[(f, level)],
            structure=(
                f"{names[level]} of family {f}: "
                + (f"{units} x {sid[(f, level - 1)]}" if level else "base class, indivisible")
                + f"; {units ** level} base units"
            ),
        )
        for f in families for level in range(n - 1, -1, -1)
    ]

    reactions: list[Reaction] = []
    index: dict[tuple, int] = {}

    def add(lhs: list[str], rhs: list[str]) -> int:
        """Append a reaction, or return the index of an identical one."""
        key = (tuple(sorted(lhs)), tuple(sorted(rhs)))
        if key not in index:
            index[key] = len(reactions)
            reactions.append(Reaction.of(lhs, rhs))     # no published kinetics
        return index[key]

    # Family by family, from the BP down: hydrolysis and its reverse, synthesis.
    for lhs, rhs in sorted(cleavages, key=lambda r: (r[0][0][0], -r[0][0][1])):
        left, right = [sid[m] for m in lhs], [sid[m] for m in rhs]
        add(left, right)
        add(right, left)

    psd_reactions = []
    if target is not None:
        for f in families:
            copies = units ** (n - 1 - target)
            # With psd_target='BM' this *is* the single hydrolysis step of the
            # model; deeper targets are one programmed process of their own.
            psd_reactions.append(reactions[add([sid[(f, n - 1)]], [sid[(f, target)]] * copies)].to_text())

    cost = {}
    for level in range(n - 1):
        syntheses, energy = rebuild_cost(level, n, units, energies)
        cost[names[level]] = {"syntheses_per_bp": syntheses, "bond_energy_per_bp": energy}

    return Network(
        species=species,
        reactions=reactions,
        status=status,
        initial_state={sid[(f, n - 1)]: float(p.initial_bp) for f in families},
        extras={
            "energies": {
                sid[(f, level)]: molecular_energy(level, units, energies)
                for f in families for level in range(n)
            },
            "conservation": [
                {
                    "name": f"base units of family {f}",
                    "vector": {sid[(f, level)]: units ** level for level in range(n)},
                }
                for f in families
            ],
            "analysis": {
                "psd_target_class": names[target] if target is not None else None,
                "rebuild_cost_by_class": cost,
                "recyclable": target is not None and target == n - 2,
            },
            "hbcb": {
                "classes": names,
                "units_per_class": units,
                "bond_energy_per_link": {
                    sid[(f, level + 1)]: energies[level]
                    for f in families for level in range(n - 1)
                },
            },
            "psd": {"target": p.psd_target, "reactions": psd_reactions},
        },
    )
