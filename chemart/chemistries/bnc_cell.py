"""Bond-number chemistry (BNC) cell model, Hintze & Adami (2008). Catalog id: bnc-cell.

Reconstructed from "Evolution of complex modular biological networks" (PLoS
Comput Biol 4(2): e23), Results and Methods, Table S1 and Figures 1 and 11:

- a molecule is a linear string of atoms 1, 2, 3; each atom carries exactly
  as many bonds as its numeral, so the bond orders ('-' single, '=' double,
  '#' triple) follow from the atoms. Strings are directed: 1-3=2 and 2=3-1
  are the distinct molecules M4 and M5. Up to 12 atoms there are 608
  molecules, numbered by length and then lexicographically (M0 = 1-1,
  M607 = 2=3-3=3-3=3-3=3-3=3-3=2); the smallest 53 (all of length <= 7) are
  the precursors;
- a reaction A + B -> A' + B' cuts one bond in each reactant and swaps the
  tails. The products are valid exactly when both cut bonds have the same
  order, and they must not exceed max_length atoms;
- mode "chemistry" returns every such reaction (topology only);
- mode "cell" returns the metabolic network of one cell: enzymes catalyse
  single reactions with the flux of eq. 3,
  v_j = [M_l][M_m] A(j) [P_j] / (kout_l kout_m), import proteins move
  precursors into the cell and export proteins move molecules out. Proteins
  carry an expression level and four 12-site affinity domains (Table S1);
  the affinity follows eq. 4 and the D(M, P) score of the Methods.
"""

from __future__ import annotations

import re
from collections import Counter
from functools import lru_cache
from itertools import product
from math import sqrt

from chemart.network import Network, Reaction, Species

ATOMS = "123"
BOND = {1: "-", 2: "=", 3: "#"}
PAPER_MOLECULES = 608          # molecules up to 12 atoms: the denominator of eq. 5
SITES = 12                     # atoms per specificity/affinity domain
MAX_DRAWS = 1_000_000
_TOKEN = re.compile(r"[123]|[-=#]")


# ----------------------------------------------------------------------------
# The chemistry
# ----------------------------------------------------------------------------
def bonds(atoms: str) -> list[int] | None:
    """Bond orders of a linear atom string, or None if it breaks the bond rule."""
    if len(atoms) < 2 or any(a not in ATOMS for a in atoms):
        return None
    out = [int(atoms[0])]
    for a in atoms[1:-1]:
        order = int(a) - out[-1]
        if not 1 <= order <= 3:
            return None
        out.append(order)
    return out if out[-1] == int(atoms[-1]) else None


def name(atoms: str) -> str:
    """'123321' -> '1-2-3=3-2-1'."""
    orders = bonds(atoms)
    return atoms[0] + "".join(BOND[o] + a for o, a in zip(orders, atoms[1:]))


@lru_cache(maxsize=None)
def molecules(max_length: int) -> tuple[str, ...]:
    """Atom strings of every valid molecule up to max_length, in the paper's order M0, M1, ..."""
    return tuple(
        "".join(s) for n in range(2, max_length + 1)
        for s in product(ATOMS, repeat=n) if bonds("".join(s)) is not None
    )


def parse(text: str) -> str:
    """'1-2-3=3-2-1' -> '123321', checking the bond symbols against the bond rule."""
    tokens = _TOKEN.findall(text)
    atoms = "".join(tokens[0::2])
    if "".join(tokens) != text or not tokens or any(t not in ATOMS for t in tokens[0::2]) \
            or any(t in ATOMS for t in tokens[1::2]) or len(tokens) % 2 == 0:
        raise ValueError(f"cannot read molecule {text!r}; write atoms 1, 2, 3 joined by bonds '-', '=', '#', e.g. '1-2-3=3-2-1'")
    if bonds(atoms) is None or name(atoms) != text:
        raise ValueError(f"{text!r} is not a valid molecule: every atom must carry exactly as many bonds as its numeral")
    return atoms


def _sites(mols: list[str]) -> list[tuple[str, str, str, int]]:
    """Every cut: (molecule, head, tail, order of the cut bond)."""
    return [(m, m[:k], m[k:], bonds(m)[k - 1]) for m in mols for k in range(1, len(m))]


def react(a_head: str, a_tail: str, b_head: str, b_tail: str) -> tuple[str, str]:
    """Swap tails: A = a_head a_tail, B = b_head b_tail -> (a_head b_tail, b_head a_tail)."""
    return a_head + b_tail, b_head + a_tail


def _is_identity(a_head: str, a_tail: str, b_head: str, b_tail: str) -> bool:
    return a_tail == b_tail or a_head == b_head


@lru_cache(maxsize=None)
def count_cleavage_reactions(max_length: int) -> int:
    """Number of ordered cleavage pairs (A cut at i, B cut at j) whose products are
    valid, not longer than max_length and not the reactants themselves: the space
    the paper counts as its 5,020,279 legal reactions (5,020,280 here, see decisions).
    """
    sites = _sites(molecules(max_length))
    shape = Counter((o, len(h), len(t)) for _, h, t, o in sites)
    total = 0
    for (o1, h1, t1), n1 in shape.items():
        for (o2, h2, t2), n2 in shape.items():
            if o1 == o2 and h1 + t2 <= max_length and h2 + t1 <= max_length:
                total += n1 * n2
    heads = Counter(h for _, h, _, _ in sites)
    tails = Counter(t for _, _, t, _ in sites)
    return total - sum(n * n for n in heads.values()) - sum(n * n for n in tails.values()) + len(sites)


def _all_reactions(mols: list[str], max_length: int) -> list[tuple[tuple[str, str], tuple[str, str]]]:
    by_order: dict[int, list] = {}
    for site in _sites(mols):
        by_order.setdefault(site[3], []).append(site)
    seen = {}
    for group in by_order.values():
        for x, (a, ah, at, _) in enumerate(group):
            for b, bh, bt, _ in group[x:]:
                if _is_identity(ah, at, bh, bt):
                    continue
                c, d = react(ah, at, bh, bt)
                if len(c) > max_length or len(d) > max_length:
                    continue
                key = (tuple(sorted((a, b))), tuple(sorted((c, d))))
                seen.setdefault(key, key)
    return list(seen)


# ----------------------------------------------------------------------------
# Proteins (Methods, Table S1)
# ----------------------------------------------------------------------------
def code(atoms: str) -> str:
    """The 12-site specificity of a molecule: atoms padded with zeros ('123321000000')."""
    return atoms.ljust(SITES, "0")


def affinity(atoms: str, domain: list[int]) -> float:
    """D(M, P) = 1 - sqrt(sum_i f(m_i EQUAL p_i)^2 / 108), f the value of the 2-bit bitwise EQUAL.

    1 for a domain that is the bitwise complement of the molecule, 0 for an identical one.
    """
    m = [int(c) for c in code(atoms)]
    return 1.0 - sqrt(sum((3 ^ (mi ^ pi)) ** 2 for mi, pi in zip(m, domain)) / 108.0)


def _parse_cut(text: str, max_length: int) -> tuple[str, str]:
    """'1-2|2-1' -> ('12', '21'): a valid molecule written with '|' in place of the cut bond."""
    bad = ValueError(f"enzyme reactant {text!r} is not a valid molecule cut at one bond (write the cut bond as '|', e.g. '1-2|2-1')")
    if text.count("|") != 1:
        raise bad
    head, tail = text.split("|")
    parts = []
    for frag in (head, tail):
        tokens = _TOKEN.findall(frag)
        if not frag or "".join(tokens) != frag or len(tokens) % 2 == 0 \
                or any(t not in ATOMS for t in tokens[0::2]) or any(t in ATOMS for t in tokens[1::2]):
            raise bad
        parts.append("".join(tokens[0::2]))
    h, t = parts
    orders = bonds(h + t)
    if orders is None or name(h + t) != f"{head}{BOND[orders[len(h) - 1]]}{tail}":
        raise bad
    if len(h + t) > max_length:
        raise ValueError(f"enzyme reactant {text!r} has {len(h + t)} atoms, above max_length = {max_length}")
    return h, t


def _parse_proteins(entries: list, mols: list[str], max_length: int) -> list[dict]:
    known = set(mols)
    out = []
    for entry in entries:
        if not isinstance(entry, str) or " " not in entry.strip():
            raise ValueError(f"proteins: cannot read {entry!r}; use 'import 1-2-1', 'export 1-1' or 'enzyme 1-2|1 + 1|2-2-1'")
        kind, spec = entry.strip().split(" ", 1)
        if kind in ("import", "export"):
            atoms = parse(spec.strip())
            if atoms not in known:
                raise ValueError(f"proteins: {spec.strip()!r} has more than max_length = {max_length} atoms")
            out.append({"type": kind, "target": atoms})
        elif kind == "enzyme":
            sides = [s.strip() for s in spec.split(" + ")]
            if len(sides) != 2:
                raise ValueError(f"proteins: enzyme {spec!r} needs two cut reactants, e.g. 'enzyme 1-2|1 + 1|2-2-1'")
            (ah, at), (bh, bt) = (_parse_cut(s, max_length) for s in sides)
            a, b = ah + at, bh + bt
            if bonds(a)[len(ah) - 1] != bonds(b)[len(bh) - 1]:
                raise ValueError(f"proteins: enzyme {spec!r} cuts bonds of different order, so its products break the bond rule")
            if _is_identity(ah, at, bh, bt):
                raise ValueError(f"proteins: enzyme {spec!r} gives back its own reactants")
            c, d = react(ah, at, bh, bt)
            if len(c) > max_length or len(d) > max_length:
                raise ValueError(f"proteins: enzyme {spec!r} makes a product longer than max_length = {max_length}")
            out.append({"type": "enzyme", "cut": (ah, at, bh, bt)})
        else:
            raise ValueError(f"proteins: unknown protein type {kind!r}; use import, export or enzyme")
    return out


def _random_proteins(p, rng, mols: list[str], precursors: list[str], max_length: int) -> list[dict]:
    out = []
    if p.enzymes:
        sites = _sites(mols)
        if count_cleavage_reactions(max_length) == 0:
            raise ValueError(f"max_length = {max_length} has no reactions to catalyse; use max_length >= 4 or enzymes = 0")
        while len(out) < p.enzymes:
            for _ in range(MAX_DRAWS):
                x, y = (int(i) for i in rng.integers(len(sites), size=2))
                _, ah, at, o1 = sites[x]
                _, bh, bt, o2 = sites[y]
                c, d = react(ah, at, bh, bt)
                if o1 == o2 and not _is_identity(ah, at, bh, bt) and len(c) <= max_length and len(d) <= max_length:
                    out.append({"type": "enzyme", "cut": (ah, at, bh, bt)})
                    break
    if p.importers and not precursors:
        raise ValueError("importers > 0 needs n_precursors > 0: import proteins carry precursors into the cell")
    for _ in range(p.importers):
        out.append({"type": "import", "target": precursors[int(rng.integers(len(precursors)))]})
    for _ in range(p.exporters):
        out.append({"type": "export", "target": mols[int(rng.integers(len(mols)))]})
    return out


# ----------------------------------------------------------------------------
def generate(p, rng):
    L = p.max_length
    mols = molecules(L)
    if p.n_precursors > len(mols):
        raise ValueError(f"n_precursors = {p.n_precursors} exceeds the {len(mols)} molecules up to max_length = {L}; the paper uses 53 of 608 (max_length 12)")
    index = {m: i for i, m in enumerate(mols)}
    precursors = mols[: p.n_precursors]
    fitness = {
        "law": "w = prod over molecules produced (1.1 + phi(M_i) Delta(M_i)); phi(M_i) = 0 for precursors, i^2 / 608^2 otherwise (eqs. 5-6)",
        "phi": {name(m): (i * i / PAPER_MOLECULES ** 2 if i >= p.n_precursors else 0.0) for i, m in enumerate(mols)},
    }
    common = {
        "precursors": [name(m) for m in precursors],
        "fitness": fitness,
    }

    if p.mode == "chemistry":
        reactions = [
            Reaction(dict(Counter(name(m) for m in lhs)), dict(Counter(name(m) for m in rhs)))
            for lhs, rhs in _all_reactions(mols, L)
        ]
        ids = [name(m) for m in mols]
        extras = {
            **common,
            "reaction_space": {
                "cleavage_reactions": count_cleavage_reactions(L),
                "distinct_reactions": len(reactions),
                "note": "cleavage_reactions counts ordered (A cut, B cut) pairs, the paper's count (5,020,279 at max_length 12); distinct_reactions merges pairs giving the same A + B -> A' + B'",
            },
            "conservation": _conservation({i: m for i, m in zip(ids, mols)}),
        }
        return Network(
            species=[Species(i, structure=code(m)) for i, m in zip(ids, mols)],
            reactions=reactions,
            extras=extras,
        )

    proteins = _parse_proteins(p.proteins, mols, L) if p.proteins else _random_proteins(p, rng, mols, precursors, L)
    return _cell(proteins, rng, mols, index, precursors, common, L)


def _conservation(atoms_of: dict[str, str]) -> list[dict]:
    return [
        {"name": f"atom {a}", "vector": {s: m.count(a) for s, m in atoms_of.items() if a in m}}
        for a in ATOMS
    ]


def _cell(proteins, rng, mols, index, precursors, common, L):
    out = lambda m: f"{name(m)}_out"   # noqa: E731
    # Draw expression levels and affinity domains in protein order.
    for k, prot in enumerate(proteins):
        prot["expression"] = float(rng.uniform(0.0, 1.0))
        prot["domains"] = [[int(v) for v in rng.integers(0, 4, SITES)] for _ in range(4)]
        prot["id"] = f"{prot['type']}_{k}"

    kout: Counter = Counter()
    for prot in proteins:
        if prot["type"] == "enzyme":
            ah, at, bh, bt = prot["cut"]
            for m in {ah + at, bh + bt}:
                kout[m] += 1
        elif prot["type"] == "export":
            kout[prot["target"]] += 1

    reactions, record = [], []
    inside, outside = set(), set()
    for prot in proteins:
        pid, dom = prot["id"], prot["domains"]
        if prot["type"] == "enzyme":
            ah, at, bh, bt = prot["cut"]
            a, b = ah + at, bh + bt
            c, d = react(ah, at, bh, bt)
            aff = sum(affinity(m, dm) for m, dm in zip((a, b, c, d), dom)) / 4.0      # eq. 4
            k = aff / (kout[a] * kout[b])                                                # eq. 3
            spec = f"{name(a)[:_cut_pos(a, ah)]}|{name(a)[_cut_pos(a, ah) + 1:]} + {name(b)[:_cut_pos(b, bh)]}|{name(b)[_cut_pos(b, bh) + 1:]}"
            reactions.append(Reaction(
                dict(Counter([name(a), name(b), pid])), dict(Counter([name(c), name(d), pid])),
                {"law": "mass-action", "k": k, "affinity": aff, "kout_product": kout[a] * kout[b]},
            ))
            inside.update((a, b, c, d))
            record.append({"id": pid, "type": "enzyme", "reaction": spec, "expression": prot["expression"], "affinity": aff, "domains": _domains(dom)})
        else:
            m = prot["target"]
            aff = sum(affinity(m, dm) for dm in dom) / 4.0
            if prot["type"] == "import":
                reactions.append(Reaction({out(m): 1, pid: 1}, {name(m): 1, pid: 1}))
            else:
                reactions.append(Reaction({name(m): 1, pid: 1}, {out(m): 1, pid: 1}))
            inside.add(m)
            outside.add(m)
            record.append({"id": pid, "type": prot["type"], "target": name(m), "expression": prot["expression"], "affinity": aff, "domains": _domains(dom)})

    order = lambda ms: sorted(ms, key=index.__getitem__)   # noqa: E731
    internal = [name(m) for m in order(inside)]
    external = [out(m) for m in order(outside)]
    pids = [prot["id"] for prot in proteins]
    species = (
        [Species(name(m), structure=code(m)) for m in order(inside)]
        + [Species(out(m), structure=code(m)) for m in order(outside)]
        + [Species(prot["id"], structure=_protein_structure(prot)) for prot in proteins]
    )
    atoms_of = {name(m): m for m in inside} | {out(m): m for m in outside}
    imported = [out(prot["target"]) for prot in proteins if prot["type"] == "import"]
    exported = [out(prot["target"]) for prot in proteins if prot["type"] == "export"]
    extras = {
        **common,
        "proteins": record,
        "buffered": pids + list(dict.fromkeys(imported)),
        "removed_every_update": list(dict.fromkeys(exported)),
        "compartments": {
            "cell": internal + pids,
            "environment": external,
            "leak": "precursors leak into the cell at 1e-6 of their concentration at the cell's location",
        },
        "space": {
            "dimensions": 2,
            "precursor_profile": "[M](d) = [M](0) exp(-d^2 / 2) / sqrt(2 pi): diffusion with D = 1/2 at t = 1 from each source (eq. 1)",
            "sources": "one randomly placed, constantly replenished source per precursor",
            "environments": {
                "static": "sources fixed",
                "quasi-static": "one random precursor source moves each update",
                "dynamic": "all sources move every update and a periodically changing 25% of precursors is unavailable",
            },
            "population": "chemostat of 1,000 cells anchored at the centre; 1 of 16 cells removed per update",
        },
        "rate_law": "eq. 3: v_j = [M_l][M_m] A(j) [P_j] / (kout_l kout_m), written as mass action with the enzyme P_j on both sides (held at its expression level) and k = A(j) / (kout_l kout_m); transport has no published rate law",
        "conservation": _conservation(atoms_of),
    }
    return Network(
        species=species,
        reactions=reactions,
        initial_state={prot["id"]: prot["expression"] for prot in proteins} or None,
        extras=extras,
    )


def _cut_pos(full: str, head: str) -> int:
    """Index in name(full) of the bond symbol after `head`."""
    return 2 * len(head) - 1


def _domains(dom: list[list[int]]) -> list[str]:
    return ["".join(map(str, d)) for d in dom]


def _protein_structure(prot: dict) -> str:
    if prot["type"] == "enzyme":
        ah, at, bh, bt = prot["cut"]
        spec = f"{ah}|{at}+{bh}|{bt}"
    else:
        spec = code(prot["target"])
    return f"{prot['type']}:{spec}:" + "-".join(_domains(prot["domains"]))
