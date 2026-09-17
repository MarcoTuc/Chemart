"""ToyChem: the graph-based toy model of chemistry of Benko, Flamm & Stadler.

Catalog id: toychem. Book 11.3.1 (one paragraph and figure 11.12); reconstructed
from Benko, Flamm & Stadler, "A graph-based toy model of chemistry", J. Chem.
Inf. Comput. Sci. 43:1085-1093 (2003) [103], from "Generic properties of chemical
networks: artificial chemistry based on graph rewriting", ECAL/LNCS 2801:10-19
(2003) [102], and from "Explicit collision simulation of chemical reactions in a
graph based artificial chemistry", ECAL 2005:725-733 [104].

Molecules. A molecule is its structural formula: a vertex- and edge-labelled
graph of atoms and bonds. The VSEPR rules [333] turn that graph, unambiguously,
into the *orbital graph* of Polansky: one vertex per outer (valence) orbital -
1s for hydrogen and sp3/sp2/sp hybrids plus the leftover p orbitals for C, N and
O - and one edge per overlap of orbitals on adjacent atoms. The steric number
(sigma bonds + lone pairs) fixes the hybridisation, so the structure formula
already carries the whole model.

Energy. A caricature of Extended Hueckel Theory that never embeds the molecule in
space: the overlap integrals S_ij are read off a table of orbital types rather
than integrated, and H is parametrised from them by Wolfsberg-Helmholtz,
H_ij = kappa (H_ii + H_jj) S_ij / 2 with H_ii = -I_i. Solving H c = E S c and
filling the lowest orbitals with two electrons each gives the total energy
E = sum_alpha n_alpha E_alpha, and the total atomization energy (TAE) is that
energy minus the valence-state energy of the separated atoms. All parameters are
Table 1 of the 2003 paper (see TABLE_1_* below).

Reactions. Graph rewriting rules on the structure formula. A rule has a left
graph, a right graph and a context, and must conserve vertex labels (atoms) and
total bond order (valence electrons) - the two conservation laws the papers
impose on every chemical rule. Bimolecular mechanisms are split into one half
rule per educt plus a join rule, as in the papers' figures 4 and 5. Implemented:
the Diels-Alder rule transcribed verbatim from appendix B of the 2003 paper, the
aldol condensation of figure 5, and keto-enol tautomerism.

Reactivity. Frontier Molecular Orbital theory: the reactivity of a channel is
inversely proportional to the frontier gap E_LUMO - E_HOMO of the two partners,
the paper's simplification Delta E = xi / (E_zeta - E_alpha) of the
Klopman-Salem formula. The gap is the network's activation energy, so reactions
carry Arrhenius rates and `barrier_cutoff` is the papers' reactivity threshold,
which is what decides how large the network grows.

Networks. The closure of a seed set under the rules ("orderly generation"),
computed by chemart.expand.expand, with molecules identified up to isomorphism
by their canonical SMILES - exactly the isomorphism test the 2003 paper uses.
The two published example networks are offered: repetitive Diels-Alder and the
formose reaction.
"""

from __future__ import annotations

import math
from collections import Counter, deque

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import BondType, RWMol
from scipy.linalg import eigh

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

RDLogger.DisableLog("rdApp.*")

NETWORKS = ("diels-alder", "formose", "custom")
REWRITE_MODES = ("all", "random", "priority")
RULE_NAMES = ("diels-alder", "aldol", "keto-enol")

#: kcal/mol per electronvolt
EV = 23.060548

SINGLE, DOUBLE = BondType.SINGLE, BondType.DOUBLE

# ============================================================================
# Table 1 of Benko, Flamm & Stadler (2003): the whole parametrisation
# ============================================================================
#: Coulomb integrals H_ii = -I_i (eV) of the hybrid (sigma) orbitals.
TABLE_1_I_SIGMA = {
    ("H", "s"): -13.6,
    ("C", "sp3"): -13.9, ("C", "sp2"): -14.5, ("C", "sp"): -15.4,
    ("N", "sp3"): -16.6, ("N", "sp2"): -17.6, ("N", "sp"): -19.7,
    ("O", "sp3"): -19.2, ("O", "sp2"): -20.6,
}
#: Coulomb integrals of the p orbitals that carry the pi system.
TABLE_1_I_PI = {"C": -11.4, "N": -13.4, "O": -14.8}

_SIGMA_ORDER = ("H_s", "C_sp3", "C_sp2", "C_sp", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2")
_SIGMA_ROWS = (
    (0.75, 0.69, 0.65, 0.66, 0.62, 0.63, 0.63, 0.55, 0.57),
    (0.69, 0.65, 0.67, 0.71, 0.60, 0.63, 0.65, 0.54, 0.57),
    (0.65, 0.67, 0.77, 0.80, 0.70, 0.73, 0.77, 0.64, 0.68),
    (0.66, 0.71, 0.80, 0.87, 0.77, 0.80, 0.84, None, None),
    (0.62, 0.60, 0.70, 0.77, 0.58, 0.61, 0.65, 0.63, 0.67),
    (0.63, 0.63, 0.73, 0.80, 0.61, 0.70, 0.73, 0.63, 0.67),
    (0.63, 0.65, 0.77, 0.84, 0.65, 0.73, 0.82, None, None),
    (0.55, 0.54, 0.64, None, 0.63, 0.63, None, None, None),
    (0.57, 0.57, 0.68, None, 0.67, 0.67, None, None, None),
)
#: sigma overlap integrals of hybrid orbitals pointing at each other along a bond.
TABLE_1_S_SIGMA = {
    (_SIGMA_ORDER[i], _SIGMA_ORDER[j]): _SIGMA_ROWS[i][j]
    for i in range(len(_SIGMA_ORDER)) for j in range(len(_SIGMA_ORDER))
    if _SIGMA_ROWS[i][j] is not None
}
#: pi overlap integrals between p orbitals on adjacent atoms.
TABLE_1_S_PI = {
    ("C", "C"): 0.38, ("C", "N"): 0.31, ("N", "C"): 0.31, ("N", "N"): 0.31,
    ("C", "O"): 0.26, ("O", "C"): 0.26, ("N", "O"): 0.26, ("O", "N"): 0.26,
    ("O", "O"): 0.26,
}
#: scaling factors of Table 1: banana bonds in strained rings, and sigma-pi coupling.
TABLE_1_RING_SCALE = {3: 0.7, 4: 0.8}
TABLE_1_HYPERCONJUGATION = 0.8
TABLE_1_KAPPA = 1.75
TABLE_1_SEMI_DIRECT = 0.1

#: valence electrons; the model is parametrised for these elements only.
VALENCE = {"H": 1, "C": 4, "N": 5, "O": 6}
_HYBRID = {2: "sp", 3: "sp2", 4: "sp3"}


# ============================================================================
# Molecules: structure formula -> orbital graph -> energy
# ============================================================================
def canonical(smiles: str) -> str:
    """The canonical SMILES of a molecule: its identity up to graph isomorphism.

    This is the isomorphism test of the 2003 paper, which transforms every new
    molecular graph into canonical SMILES so that "the isomorphism test then
    reduces to simple string comparison".
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot read the structure formula {smiles!r} as SMILES")
    return Chem.MolToSmiles(mol)


def _molecule(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot read the structure formula {smiles!r} as SMILES")
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in VALENCE:
            raise ValueError(
                f"{smiles!r} contains {atom.GetSymbol()}, which Table 1 of the toy model does "
                f"not parametrise; the implemented elements are {sorted(VALENCE)}")
        if atom.GetFormalCharge():
            raise ValueError(
                f"{smiles!r} is charged; the toy model's energy calculation is limited to "
                "neutral molecules (2003 paper, section 2)")
        if atom.GetNumRadicalElectrons():
            raise ValueError(f"{smiles!r} is a radical; the toy model treats closed shells only")
    return mol


def orbital_graph(smiles: str) -> dict:
    """The Polansky orbital graph of a structure formula.

    Returns the orbital list ``[(atom, element, kind)]`` with kind in
    {s, sp3, sp2, sp, p}, the hybridisation of every atom, and the bookkeeping
    the energy calculation needs. Hybridisation follows VSEPR: the steric number
    (sigma bonds + lone pairs) is 4, 3 or 2 for sp3, sp2 or sp.
    """
    mol = Chem.AddHs(_molecule(smiles))
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)

    orbitals: list[tuple[int, str, str]] = []
    hybrid: dict[int, str] = {}
    hybrids_of: dict[int, list[int]] = {}
    p_of: dict[int, list[int]] = {}
    for atom in kekule.GetAtoms():
        index, element = atom.GetIdx(), atom.GetSymbol()
        if element == "H":
            hybrid[index] = "s"
            hybrids_of.setdefault(index, []).append(len(orbitals))
            orbitals.append((index, "H", "s"))
            continue
        order = sum(int(b.GetBondTypeAsDouble()) for b in atom.GetBonds())
        lone_pairs = (VALENCE[element] - order) // 2
        steric = atom.GetDegree() + lone_pairs
        if steric not in _HYBRID:
            raise ValueError(
                f"{smiles!r}: atom {index} ({element}) has steric number {steric}, which VSEPR "
                "does not map to an sp3, sp2 or sp hybridisation in this model")
        hybrid[index] = _HYBRID[steric]
        for _ in range(steric):
            hybrids_of.setdefault(index, []).append(len(orbitals))
            orbitals.append((index, element, _HYBRID[steric]))
        for _ in range(4 - steric):
            p_of.setdefault(index, []).append(len(orbitals))
            orbitals.append((index, element, "p"))
    return {"mol": kekule, "orbitals": orbitals, "hybrid": hybrid,
            "hybrids_of": hybrids_of, "p_of": p_of}


def _overlap_key(element: str, hyb: str) -> str:
    return "H_s" if element == "H" else f"{element}_{hyb}"


def eht(smiles: str, kappa: float = TABLE_1_KAPPA,
        semi_direct: float = TABLE_1_SEMI_DIRECT) -> dict:
    """Solve the toy model's Extended Hueckel problem for one molecule.

    Returns the total energy, the total atomization energy (TAE), the frontier
    orbital energies and the orbital spectrum, all in eV except ``tae_kcal``.
    """
    graph = orbital_graph(smiles)
    mol, orbitals = graph["mol"], graph["orbitals"]
    hybrid, hybrids_of, p_of = graph["hybrid"], graph["hybrids_of"], graph["p_of"]
    n = len(orbitals)

    # one hybrid orbital per neighbour; the rest hold lone pairs
    pointing: dict[tuple[int, int], int] = {}
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        for slot, neighbour in enumerate(atom.GetNeighbors()):
            pointing[(i, neighbour.GetIdx())] = hybrids_of[i][slot]

    rings = mol.GetRingInfo().AtomRings()
    S = np.eye(n)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ei = mol.GetAtomWithIdx(i).GetSymbol()
        ej = mol.GetAtomWithIdx(j).GetSymbol()
        smallest = min((len(r) for r in rings if i in r and j in r), default=0)
        scale = TABLE_1_RING_SCALE.get(smallest, 1.0)       # banana bonds
        key = (_overlap_key(ei, hybrid[i]), _overlap_key(ej, hybrid[j]))
        base = TABLE_1_S_SIGMA.get(key)
        if base is None:
            raise ValueError(
                f"{smiles!r}: Table 1 gives no sigma overlap for {key[0]}/{key[1]}")
        oi, oj = pointing[(i, j)], pointing[(j, i)]
        S[oi, oj] = S[oj, oi] = base * scale
        # "semi-direct": exactly one of the two orbitals points along the bond.
        # Table 1 sets the semi-direct overlap with hydrogen, and every indirect
        # overlap, to zero.
        if ei != "H" and ej != "H":
            for other in hybrids_of[i]:
                if other != oi:
                    S[other, oj] = S[oj, other] = semi_direct * base * scale
            for other in hybrids_of[j]:
                if other != oj:
                    S[oi, other] = S[other, oi] = semi_direct * base * scale
        # pi overlaps: one per bond order above one
        for extra in range(int(bond.GetBondTypeAsDouble()) - 1):
            pi_i, pi_j = p_of.get(i, []), p_of.get(j, [])
            if extra < len(pi_i) and extra < len(pi_j):
                S[pi_i[extra], pi_j[extra]] = S[pi_j[extra], pi_i[extra]] = (
                    TABLE_1_S_PI[(ei, ej)] * scale)
        # hyperconjugation: a p orbital against an sp3 orbital of a neighbour
        for (a, b) in ((i, j), (j, i)):
            ea = mol.GetAtomWithIdx(a).GetSymbol()
            eb = mol.GetAtomWithIdx(b).GetSymbol()
            if ea == "H" or eb == "H" or hybrid[b] != "sp3":
                continue
            for pa in p_of.get(a, []):
                target = pointing[(b, a)]
                if S[pa, target]:
                    continue
                S[pa, target] = S[target, pa] = (
                    TABLE_1_HYPERCONJUGATION * TABLE_1_S_PI[(ea, eb)] * scale)

    diagonal = np.array([
        TABLE_1_I_PI[element] if kind == "p" else TABLE_1_I_SIGMA[(element, kind)]
        for _atom, element, kind in orbitals])
    H = np.diag(diagonal).astype(float)
    for a in range(n):
        for b in range(n):
            if a != b and S[a, b]:
                H[a, b] = kappa * (diagonal[a] + diagonal[b]) * S[a, b] / 2.0

    levels = eigh(H, S, eigvals_only=True)
    electrons = sum(VALENCE[a.GetSymbol()] for a in mol.GetAtoms())
    if electrons % 2:
        raise ValueError(f"{smiles!r} has an odd number of valence electrons (open shell)")
    occupied = electrons // 2
    total = 2.0 * float(np.sum(levels[:occupied]))
    # the separated atoms in their valence state: one electron per valence orbital
    reference = float(np.sum(diagonal))
    homo = float(levels[occupied - 1])
    lumo = float(levels[occupied]) if occupied < n else float("inf")
    return {
        "total_energy": total,
        "tae": total - reference,
        "tae_kcal": (total - reference) * EV,
        "homo": homo,
        "lumo": lumo,
        "gap": lumo - homo,
        "orbitals": len(orbitals),
        "electrons": electrons,
    }


# ============================================================================
# Reactions: graph rewriting on the structure formula
# ============================================================================
#: Appendix B of the 2003 paper, verbatim: context six carbons; left graph
#: 1=2, 2-3, 3=4, 5=6; right graph 1-2, 2=3, 3-4, 4-5, 5-6, 6-1.
_DIENE = Chem.MolFromSmarts("[#6;!a]=[#6;!a]-[#6;!a]=[#6;!a]")
_DIENOPHILE = Chem.MolFromSmarts("[#6;!a]=[#6;!a]")
#: Figure 5: the enol half rule (C=C-O-H) and the carbonyl half rule (C=O).
_ENOL = Chem.MolFromSmarts("[#6;!a]=[#6;!a]-[OX2]-[H]")
_CARBONYL = Chem.MolFromSmarts("[#6;!a]=[OX1]")
_KETO = Chem.MolFromSmarts("[H]-[#6;!a]-[#6;!a]=[OX1]")


def _host(parts: list[str]) -> RWMol:
    """The disjoint union of the educts, with explicit hydrogens, ready to rewrite."""
    combined = Chem.AddHs(_molecule(parts[0]))
    for extra in parts[1:]:
        combined = Chem.CombineMols(combined, Chem.AddHs(_molecule(extra)))
    return RWMol(combined)


def _fragments(rw: RWMol) -> tuple[str, ...] | None:
    """Sanitise the rewritten graph and return the canonical SMILES of its components."""
    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
    except Exception:                                          # noqa: BLE001
        return None
    pieces = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return tuple(sorted(Chem.MolToSmiles(Chem.RemoveHs(piece)) for piece in pieces))


def _carbons(smiles: str) -> int:
    return sum(1 for a in Chem.MolFromSmiles(smiles).GetAtoms() if a.GetSymbol() == "C")


def diels_alder(parts: list[str]) -> list[tuple[str, ...]]:
    """[2+4] cycloaddition, the rule of appendix B.

    With one educt this is the intramolecular rearrangement of figure 4; with two
    it is the intermolecular reaction, the rule being applied to the disjoint
    union of the two graphs.
    """
    host = _host(parts)
    out: set[tuple[str, ...]] = set()
    for a, b, c, d in host.GetSubstructMatches(_DIENE, uniquify=False):
        for e, f in host.GetSubstructMatches(_DIENOPHILE, uniquify=False):
            if {a, b, c, d} & {e, f}:
                continue
            rw = RWMol(host)
            rw.GetBondBetweenAtoms(a, b).SetBondType(SINGLE)
            rw.GetBondBetweenAtoms(b, c).SetBondType(DOUBLE)
            rw.GetBondBetweenAtoms(c, d).SetBondType(SINGLE)
            rw.GetBondBetweenAtoms(e, f).SetBondType(SINGLE)
            if rw.GetBondBetweenAtoms(d, e) or rw.GetBondBetweenAtoms(f, a):
                continue
            rw.AddBond(d, e, SINGLE)
            rw.AddBond(f, a, SINGLE)
            products = _fragments(rw)
            if products:
                out.add(products)
    return sorted(out)


def keto_enol(parts: list[str]) -> list[tuple[str, ...]]:
    """Keto-enol tautomerism, both ways: H-C-C=O <-> C=C-O-H.

    The formose network needs it: an aldose has to enolise before it can act as
    the nucleophile of an aldol condensation, which is the isomerisation the 2003
    ECAL paper calls "toggling the state of a species from non-reactive to
    reactive".
    """
    if len(parts) != 1:
        return []
    host = _host(parts)
    out: set[tuple[str, ...]] = set()
    for h, c2, c3, o in host.GetSubstructMatches(_KETO, uniquify=False):
        rw = RWMol(host)
        rw.RemoveBond(h, c2)
        rw.GetBondBetweenAtoms(c2, c3).SetBondType(DOUBLE)
        rw.GetBondBetweenAtoms(c3, o).SetBondType(SINGLE)
        rw.AddBond(o, h, SINGLE)
        products = _fragments(rw)
        if products:
            out.add(products)
    for c1, c2, o, h in host.GetSubstructMatches(_ENOL, uniquify=False):
        rw = RWMol(host)
        rw.RemoveBond(o, h)
        rw.GetBondBetweenAtoms(c1, c2).SetBondType(SINGLE)
        rw.GetBondBetweenAtoms(c2, o).SetBondType(DOUBLE)
        rw.AddBond(c1, h, SINGLE)
        products = _fragments(rw)
        if products:
            out.add(products)
    return sorted(out)


def aldol(parts: list[str], max_carbons: int = 0) -> list[tuple[str, ...]]:
    """Aldol condensation, the three rules of figure 5.

    Half rule 1 acts on the carbonyl educt (C=O becomes C-O-H and the carbon is
    flagged), half rule 2 on the enol educt (C=C-O-H becomes C-C=O and the
    nucleophilic carbon is flagged), and the join rule closes the new C-C bond
    between the two flagged carbons. The product is a beta-hydroxy carbonyl.

    ``max_carbons`` reproduces the cut of the 2003 paper's figure 8: "in order to
    account for cyclisation, which limits the network, we do not permit carbon
    chains with more than four members to undergo further aldol condensation".
    """
    if len(parts) != 2:
        return []
    host = _host(parts)
    first = Chem.AddHs(_molecule(parts[0])).GetNumAtoms()
    out: set[tuple[str, ...]] = set()
    for c1, c2, o3, h4 in host.GetSubstructMatches(_ENOL, uniquify=False):
        enol_in_first = max(c1, c2, o3, h4) < first
        for c5, o6 in host.GetSubstructMatches(_CARBONYL, uniquify=False):
            carbonyl_in_first = c5 < first
            if enol_in_first == carbonyl_in_first:          # one partner each
                continue
            if len({c1, c2, o3, h4} & {c5, o6}):
                continue
            rw = RWMol(host)
            rw.GetBondBetweenAtoms(c1, c2).SetBondType(SINGLE)
            rw.GetBondBetweenAtoms(c2, o3).SetBondType(DOUBLE)
            rw.RemoveBond(o3, h4)
            rw.GetBondBetweenAtoms(c5, o6).SetBondType(SINGLE)
            rw.AddBond(o6, h4, SINGLE)
            if rw.GetBondBetweenAtoms(c1, c5):
                continue
            rw.AddBond(c1, c5, SINGLE)
            products = _fragments(rw)
            if not products:
                continue
            if max_carbons and any(_carbons(s) > max_carbons for s in products):
                continue
            out.add(products)
    return sorted(out)


def apply_rule(name: str, parts: list[str], max_carbons: int = 0) -> list[tuple[str, ...]]:
    """Every outcome of one named rule on one ordered list of educts."""
    if name == "diels-alder":
        return diels_alder(parts)
    if name == "keto-enol":
        return keto_enol(parts)
    if name == "aldol":
        return aldol(parts, max_carbons)
    raise ValueError(f"unknown rule {name!r}; the implemented rules are {list(RULE_NAMES)}")


def conserves(educts: list[str], products: tuple[str, ...]) -> bool:
    """The two conservation laws every chemical rewrite rule must satisfy.

    Conservation of vertex labels (the atoms) and of total bond order (the
    valence electrons), checked by comparing the label lists of the left and the
    right graph, exactly as section 3 of the 2003 paper prescribes.
    """
    def census(smiles_list):
        atoms: Counter = Counter()
        order = 0
        for smiles in smiles_list:
            mol = Chem.AddHs(_molecule(smiles))
            kekule = Chem.Mol(mol)
            Chem.Kekulize(kekule, clearAromaticFlags=True)
            atoms.update(a.GetSymbol() for a in kekule.GetAtoms())
            order += sum(int(b.GetBondTypeAsDouble()) for b in kekule.GetBonds())
        return atoms, order

    return census(educts) == census(list(products))


# ============================================================================
# The published example networks
# ============================================================================
#: Figure 7 of the 2003 paper: "The initial mixture consists of cyclobutadiene,
#: ethenol, phtalic anhydride, methylbutadiene, and cyclohexa-1,3-diene."
DIELS_ALDER_SEED = ("C1=CC=C1", "C=CO", "O=C1OC(=O)c2ccccc21", "C=CC(C)=C", "C1=CC=CCC1")
#: Figure 8: "The initial mixture consists of formaldehyde H2CO and glycol
#: aldehyde CH2OH-CHO and reacts via aldol condensations and dismutations."
FORMOSE_SEED = ("C=O", "OCC=O")

EXAMPLES = {
    "diels-alder": {"seed": DIELS_ALDER_SEED, "rules": ("diels-alder",)},
    "formose": {"seed": FORMOSE_SEED, "rules": ("keto-enol", "aldol")},
}


def build(p) -> dict:
    """The seed set and the rule set of the chosen network."""
    if p.network != "custom":
        if p.seed_molecules or p.rules:
            raise ValueError("seed_molecules and rules are only used with network='custom'")
        spec = EXAMPLES[p.network]
        return {"seed": list(spec["seed"]), "rules": list(spec["rules"])}
    if (not isinstance(p.seed_molecules, list) or not p.seed_molecules
            or not all(isinstance(s, str) and s for s in p.seed_molecules)):
        raise ValueError("network 'custom' needs seed_molecules: a non-empty list of SMILES "
                         f"structure formulae, e.g. ['C=CC=C', 'C=C']; got {p.seed_molecules!r}")
    if (not isinstance(p.rules, list) or not p.rules
            or not all(r in RULE_NAMES for r in p.rules)):
        raise ValueError(f"network 'custom' needs rules: a non-empty list drawn from "
                         f"{list(RULE_NAMES)}; got {p.rules!r}")
    return {"seed": list(p.seed_molecules), "rules": list(p.rules)}


# ============================================================================
# Generator
# ============================================================================
class Engine:
    """Applies the rule set, scores every channel, and gates it by its barrier."""

    def __init__(self, p, rules: list[str], rng):
        self.p, self.rules, self.rng = p, rules, rng
        self._energy: dict[str, dict] = {}
        self.found: dict[tuple, dict] = {}

    def energy(self, smiles: str) -> dict:
        if smiles not in self._energy:
            self._energy[smiles] = eht(smiles, self.p.kappa, self.p.semi_direct_scale)
        return self._energy[smiles]

    def frontier_gap(self, educts: list[str]) -> float:
        """E_LUMO - E_HOMO over the partners: the FMO denominator, in eV.

        Frontier Molecular Orbital theory takes the reactivity to be inversely
        proportional to the gap between the HOMO of one system and the LUMO of
        the other; the reactive channel is the one with the smaller gap.
        """
        levels = [self.energy(s) for s in educts]
        if len(levels) == 1:                                 # intramolecular: A = B
            return levels[0]["lumo"] - levels[0]["homo"]
        return min(b["lumo"] - a["homo"]
                   for i, a in enumerate(levels)
                   for j, b in enumerate(levels) if i != j)

    def react(self, *educts: str):
        lhs = list(educts)
        outcomes: dict[tuple, list] = {}
        for rule in self.rules:
            for products in apply_rule(rule, lhs, self.p.max_carbons):
                if not conserves(lhs, products):
                    continue
                entry = outcomes.setdefault(products, [])
                if rule not in entry:
                    entry.append(rule)
        if not outcomes:
            return None

        gap = self.frontier_gap(lhs)
        if gap <= 0:
            return None
        barrier = gap * EV                                   # kcal/mol
        if self.p.barrier_cutoff and barrier > self.p.barrier_cutoff:
            return None
        reactivity = self.p.xi / gap

        channels = sorted(outcomes)
        if self.p.rewrite_mode == "random" and len(channels) > 1:
            channels = [channels[int(self.rng.integers(len(channels)))]]
        elif self.p.rewrite_mode == "priority" and len(channels) > 1:
            # the paper's priority mode: keep the channel the model prefers, i.e.
            # the most exothermic one
            channels = [min(channels, key=lambda c: self._reaction_energy(lhs, c))]

        left = tuple(sorted(lhs))
        for products in channels:
            self.found[(left, tuple(sorted(products)))] = {
                "rules": outcomes[products],
                "barrier": barrier,
                "reactivity": reactivity,
                "reaction_energy": self._reaction_energy(lhs, products),
            }
        return channels

    def _reaction_energy(self, educts: list[str], products: tuple[str, ...]) -> float:
        return (sum(self.energy(s)["tae_kcal"] for s in products)
                - sum(self.energy(s)["tae_kcal"] for s in educts))


def substrate_graph(species: list[str], reactions: list[tuple]) -> dict:
    """The one-mode projection of the reaction hypergraph, and its statistics.

    The substrate graph of the 2003 ECAL paper: the species are the vertices and
    two species are joined when they take part in the same reaction, i.e. every
    hyperedge is replaced by a clique. Reported with the Erdos-Renyi references
    <L_rand> = ln n / ln <k> and <C_rand> = <k> / (n - 1) of its equation (5).
    """
    index = {s: i for i, s in enumerate(species)}
    adjacency: dict[int, set[int]] = {i: set() for i in range(len(species))}
    for left, right in reactions:
        members = {index[m] for m in left} | {index[m] for m in right}
        for u in members:
            adjacency[u] |= members - {u}

    n = len(species)
    m = sum(len(a) for a in adjacency.values()) // 2
    degree = 2 * m / n if n else 0.0

    clustering = []
    for u, neighbours in adjacency.items():
        d = len(neighbours)
        if d < 2:
            continue
        q = sum(1 for x in neighbours for y in neighbours if x < y and y in adjacency[x])
        clustering.append(2 * q / (d * (d - 1)))
    mean_clustering = sum(clustering) / len(clustering) if clustering else 0.0

    total, pairs = 0, 0
    for start in adjacency:
        seen = {start: 0}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if v not in seen:
                    seen[v] = seen[u] + 1
                    queue.append(v)
        for v, distance in seen.items():
            if v != start:
                total += distance
                pairs += 1
    mean_path = total / pairs if pairs else 0.0

    random_path = (math.log(n) / math.log(degree)) if n > 1 and degree > 1 else 0.0
    random_clustering = degree / (n - 1) if n > 1 else 0.0
    return {
        "nodes": n,
        "edges": m,
        "mean_degree": degree,
        "mean_path_length": mean_path,
        "mean_clustering": mean_clustering,
        "random_path_length": random_path,
        "random_clustering": random_clustering,
        "small_world": bool(mean_clustering > random_clustering and mean_path <= random_path
                            and random_path > 0),
        "degree_sequence": sorted((len(a) for a in adjacency.values()), reverse=True),
    }


def generate(p, rng):
    if p.network not in NETWORKS:
        raise ValueError(f"network must be one of {NETWORKS}, got {p.network!r}")
    if p.rewrite_mode not in REWRITE_MODES:
        raise ValueError(f"rewrite_mode must be one of {REWRITE_MODES}, got {p.rewrite_mode!r}")
    spec = build(p)
    seed = [canonical(s) for s in spec["seed"]]
    if len(set(seed)) != len(seed):
        raise ValueError(f"the seed set contains the same molecule twice: {spec['seed']}")

    engine = Engine(p, spec["rules"], rng)
    species, pairs, status = expand(engine.react, seed, arity=(1, 2),
                                    max_species=p.max_species, ordered=False,
                                    alternatives=True)

    reactions, rule_labels = [], []
    for left, right in pairs:
        record = engine.found[(tuple(sorted(left)), tuple(sorted(right)))]
        rate = {
            "law": "arrhenius",
            "A": float(p.arrhenius_prefactor),
            "Ea": float(record["barrier"]),
            "reactivity": float(record["reactivity"]),
            "rules": ",".join(record["rules"]),
        }
        reactions.append(Reaction(dict(Counter(left)), dict(Counter(right)), rate))
        rule_labels.append(list(record["rules"]))

    levels = {s: engine.energy(s) for s in species}
    elements = sorted({a.GetSymbol() for s in species
                       for a in Chem.AddHs(_molecule(s)).GetAtoms()})
    conservation = []
    for element in elements:
        vector = {s: sum(1 for a in Chem.AddHs(_molecule(s)).GetAtoms()
                         if a.GetSymbol() == element) for s in species}
        if any(vector.values()):
            conservation.append({"name": f"atoms of type {element}", "vector": vector})

    stats = substrate_graph(species, pairs)
    return Network(
        species=[Species(s, structure=s) for s in species],
        reactions=reactions,
        status=status,
        extras={
            "network": p.network,
            "seed": list(seed),
            "rules": list(spec["rules"]),
            "reaction_rules": rule_labels,
            "energies": {
                "units": "kcal/mol",
                "model": "caricature Extended Hueckel Theory on the orbital graph: overlap "
                         "integrals read off Table 1 of Benko, Flamm & Stadler (2003) by orbital "
                         "type instead of being integrated in space, H_ij = kappa (H_ii + H_jj) "
                         "S_ij / 2 (Wolfsberg-Helmholtz) with H_ii = -I_i, then H c = E S c and "
                         "E = sum_alpha n_alpha E_alpha over doubly occupied orbitals",
                "reference": "total atomization energy: the molecule's energy minus the "
                             "valence-state energy of its separated atoms",
                "total_atomization_energy": {s: float(v["tae_kcal"]) for s, v in levels.items()},
                "homo_ev": {s: float(v["homo"]) for s, v in levels.items()},
                "lumo_ev": {s: float(v["lumo"]) for s, v in levels.items()},
                "reaction_energy": [
                    float(engine.found[(tuple(sorted(a)), tuple(sorted(b)))]["reaction_energy"])
                    for a, b in pairs],
                "activation_energy": [
                    float(engine.found[(tuple(sorted(a)), tuple(sorted(b)))]["barrier"])
                    for a, b in pairs],
                "activation_model":
                    "Frontier Molecular Orbital theory: the reactivity of a channel is inversely "
                    "proportional to the frontier gap E_LUMO - E_HOMO of the two partners (the "
                    "paper's simplification Delta E = xi / (E_zeta - E_alpha) of the "
                    "Klopman-Salem formula). The gap is used as the Arrhenius activation energy "
                    "and xi / gap is reported as the paper's reactivity index.",
            },
            "conservation": conservation,
            "parameters": {
                "kappa": float(p.kappa),
                "semi_direct_scale": float(p.semi_direct_scale),
                "hyperconjugation": TABLE_1_HYPERCONJUGATION,
                "ring_scale": {str(k): v for k, v in TABLE_1_RING_SCALE.items()},
                "xi": float(p.xi),
                "temperature": float(p.temperature),
                "barrier_cutoff": float(p.barrier_cutoff),
            },
            "species_encoding":
                "id = structure = canonical SMILES of the structure formula, so molecules are "
                "identified up to graph isomorphism by string comparison, as in the 2003 paper",
            "rewrite_rules":
                "graph rewriting on the structure formula, conserving vertex labels (atoms) and "
                "total bond order (valence electrons); bimolecular mechanisms are a half rule per "
                "educt plus a join rule",
            "analysis": {
                "substrate_graph": stats,
                "rewrite_mode": p.rewrite_mode,
                "max_carbons": int(p.max_carbons),
            },
        },
    )
