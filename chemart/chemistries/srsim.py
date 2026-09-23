"""SRSim: spatial rule-based simulation of diffusing, geometrically constrained molecules.

Catalog id: srsim. Book 18.3.3 (one paragraph); reconstructed from Gruenert,
Ibrahim, Lenser, Lohel, Hinze & Dittrich, "Rule-based spatial modeling with
diffusing, geometrically constrained molecules", BMC Bioinformatics 11:307
(2010) [351], from its additional file 2 ("Calculation of Kinetic Parameters
from Macroscopic Values"), and from the released SRSim sources and example
models (additional files 1 and 3 of the same paper).

The model. An elementary molecule (EM) is a sphere with a mass, a radius, a
diffusion coefficient and named components (binding sites), each given in
polar coordinates (dist, theta, phi) from the sphere centre and each
optionally carrying an internal state (u/p, nP/P, ...). A bond is the straight
connection of two component vectors, so its ideal length is d_ij + d_kl, and
the ideal angle between two components of one EM follows from their polar
coordinates. A species is a complex: a connected site graph, identified up to
isomorphism by a canonical text (`canonical_text`).

Dynamics. EMs move in a cuboid periodic box by Langevin dynamics (the paper's
equation for F_i) or in its overdamped Brownian limit, under three potentials:
harmonic bonds E_d = K_d (d - d_ij - d_kl)^2, harmonic angles
E_a = K_a (alpha - alpha_ijk)^2 and the soft-sphere repulsion
E_s = A [1 + cos(pi r / r_c)] with r_c = r_i + r_j.

Reactions. After every position update the rule system looks for reactions.
A bimolecular rule may fire between two EMs only if they are geometrically
compatible: their distance deviates by less than `distance_tolerance` from the
ideal bond length and - this is SRSim's own test, the routine testSiteGeo of
its LAMMPS module - for every bond either molecule already carries, the angle
between the direction to the partner and that bond deviates by less than
`angle_tolerance` from the ideal angle. An EM without bonds has no rotational
orientation and always passes (with `orientation="sampled"` it instead passes
with the probability (1 - cos t_ang)/2 that a random orientation offers the
component, which makes the reactive volume of the paper's supplement exact).
A compatible pair reacts with probability 1 - exp(-k dt); when a bond breaks,
both molecules stay refractory for `refractory_time`, which is how SRSim
suppresses geminate recombination.

The chemistry is a gas with one face, ``evolve``: a frame per
molecular-dynamics step, at the simulated time, whose state counts the
complexes. The network is the observed set of complex-level reactions with
firing counts (status "observed"); the geometry, the reactor and the final configuration are
in extras["space"], the force law in extras["interaction_law"], the rule set
and its kinetics in extras["rules"] / extras["kinetics"], and the assembly
statistics in extras["analysis"].
"""

from __future__ import annotations

import hashlib
import math
import re
from collections import Counter

import numpy as np

from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame, ticks

MODELS = ("polymer", "dimerization", "scaffold", "custom")
INTEGRATORS = ("langevin", "brownian")
ORIENTATIONS = ("bound-only", "sampled")

#: Rates of the paper's scaffold example (additional file 3, 001.scaffold/spass_w.bngl).
SCAFFOLD_RATES = {"d1_on": 1.0e-4, "d1_off": 0.1, "d2_on": 0.05, "d2_off": 1.0e-3,
                  "d3_on": 1.0e-5, "d3_off": 0.9}
#: ... and the rate prefactors its LAMMPS script gives the fix srsim
#: (preFactBindR, preFactBreakR, preFactModifyR_1 monomolecular, preFactModifyR_2 bimolecular).
SCAFFOLD_PREFACTORS = {"bind": 4.7, "break": 1.0e-3, "modify1": 1.0e-3, "modify2": 4.7}

_SITE = re.compile(r"^([A-Za-z]\w*)(?:~(\w+))?(?:!(\d+|\+|\?))?$")
_EM = re.compile(r"^([A-Za-z]\w*)\((.*)\)$")
_LABEL = re.compile(r"^'([^']*)'\s*(.*)$", re.S)


# ============================================================================
# Molecule types
# ============================================================================
class Site:
    """One component of an elementary molecule, in polar coordinates."""

    __slots__ = ("name", "dist", "theta", "phi", "states", "unit")

    def __init__(self, name: str, dist: float, theta: float, phi: float, states: tuple):
        self.name, self.dist = name, float(dist)
        self.theta, self.phi = float(theta), float(phi)
        self.states = states
        t, f = math.radians(self.theta), math.radians(self.phi)
        self.unit = (math.sin(t) * math.cos(f), math.sin(t) * math.sin(f), math.cos(t))


class MolType:
    """An elementary species: a sphere plus its components."""

    __slots__ = ("name", "radius", "mass", "diffusion", "sites", "ideal", "index")

    def __init__(self, name, radius, mass, diffusion, sites, index):
        self.name, self.radius, self.mass = name, float(radius), float(mass)
        self.diffusion, self.sites, self.index = float(diffusion), sites, index
        n = len(sites)
        self.ideal = [[0.0] * n for _ in range(n)]
        for a in range(n):
            for b in range(n):
                dot = sum(x * y for x, y in zip(sites[a].unit, sites[b].unit))
                self.ideal[a][b] = math.degrees(math.acos(max(-1.0, min(1.0, dot))))

    def sites_named(self, name: str) -> list[int]:
        return [i for i, s in enumerate(self.sites) if s.name == name]


def _number(value, what: str, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{what} must be a number, got {value!r}")
    if positive and value <= 0:
        raise ValueError(f"{what} must be positive, got {value!r}")
    return float(value)


def parse_types(specs) -> list[MolType]:
    """[{'name': 'M', 'radius': 0.5, 'sites': [{'name': 'a', 'dist': 1, 'theta': 90, 'phi': 0}]}]."""
    if not isinstance(specs, list) or not specs:
        raise ValueError("molecule_types must be a non-empty list of elementary species, e.g. "
                         "[{'name': 'M', 'radius': 0.5, 'mass': 1.0, 'diffusion': 1.0, "
                         "'sites': [{'name': 'a', 'dist': 1.0, 'theta': 90, 'phi': 0}]}]")
    out: list[MolType] = []
    seen: set[str] = set()
    for spec in specs:
        if not isinstance(spec, dict) or not isinstance(spec.get("name"), str) or not spec["name"]:
            raise ValueError(f"every molecule type needs a name, got {spec!r}")
        name = spec["name"]
        if not re.fullmatch(r"[A-Za-z]\w*", name):
            raise ValueError(f"molecule type name {name!r} must be a word starting with a letter")
        if name in seen:
            raise ValueError(f"molecule type {name!r} is declared twice")
        seen.add(name)
        raw_sites = spec.get("sites")
        if not isinstance(raw_sites, list) or not raw_sites:
            raise ValueError(f"molecule type {name!r} needs a non-empty list of sites")
        sites = []
        for s in raw_sites:
            if not isinstance(s, dict) or not isinstance(s.get("name"), str) or not s["name"]:
                raise ValueError(f"every site of {name!r} needs a name, got {s!r}")
            if not re.fullmatch(r"[A-Za-z]\w*", s["name"]):
                raise ValueError(f"site name {s['name']!r} of {name!r} must be a word")
            states = s.get("states", [])
            if not isinstance(states, list) or not all(isinstance(v, str) and v for v in states):
                raise ValueError(f"states of site {name}.{s['name']} must be a list of strings")
            if len(set(states)) != len(states):
                raise ValueError(f"site {name}.{s['name']} lists an internal state twice")
            sites.append(Site(s["name"], _number(s.get("dist", 1.0), f"dist of {name}.{s['name']}"),
                              _number(s.get("theta", 0.0), f"theta of {name}.{s['name']}", False),
                              _number(s.get("phi", 0.0), f"phi of {name}.{s['name']}", False),
                              tuple(states)))
        out.append(MolType(name, _number(spec.get("radius", 1.0), f"radius of {name!r}"),
                           _number(spec.get("mass", 1.0), f"mass of {name!r}"),
                           _number(spec.get("diffusion", 1.0), f"diffusion of {name!r}", False),
                           sites, len(out)))
    return out


# ============================================================================
# Rules: a BNGL-flavoured mini-language
# ============================================================================
class Pattern:
    """One elementary-molecule pattern of a rule: a type and its site constraints."""

    __slots__ = ("type", "constraints")

    def __init__(self, type_name: str, constraints: list):
        self.type, self.constraints = type_name, constraints    # (site, state|None, link)


class Rule:
    """A one-way rule: bind, unbind or modify, with a microscopic rate."""

    __slots__ = ("label", "kind", "k", "patterns", "contact", "change", "text",
                 "options", "d0", "v_react")

    def __init__(self, label, kind, k, patterns, contact, change, text):
        self.label, self.kind, self.k = label, kind, float(k)
        self.patterns, self.contact, self.change, self.text = patterns, contact, change, text
        self.options: list[list[tuple]] = []
        self.d0 = 0.0
        self.v_react = 0.0

    @property
    def order(self) -> int:
        return len(self.patterns)


def _parse_em(text: str, where: str) -> Pattern:
    m = _EM.match(text.strip())
    if not m:
        raise ValueError(f"cannot read the molecule pattern {text.strip()!r} in rule {where!r}: "
                         "expected Name(site, site~state, site!1, ...)")
    constraints = []
    for tok in (t.strip() for t in m.group(2).split(",")):
        if not tok:
            continue
        sm = _SITE.match(tok)
        if not sm:
            raise ValueError(f"cannot read the site {tok!r} of {text.strip()!r} in rule {where!r}")
        link = sm.group(3)
        if link is not None and link.isdigit():
            link = ("bond", int(link))
        constraints.append((sm.group(1), sm.group(2), link))
    return Pattern(m.group(1), constraints)


def _parse_side(text: str, where: str) -> list[Pattern]:
    out = []
    for component in text.split(" + "):
        for em in component.split("."):
            if em.strip():
                out.append(_parse_em(em, where))
    return out


def _one_way(label, lhs, rhs, k, text, types) -> Rule:
    by_name = {t.name: t for t in types}
    if len(lhs) != len(rhs):
        raise ValueError(f"rule {text!r}: both sides must list the same elementary molecules "
                         "(SRSim rules never create or delete molecules)")
    if not 1 <= len(lhs) <= 2:
        raise ValueError(f"rule {text!r}: SRSim executes mono- and bimolecular reactions, so a rule "
                         f"has one or two elementary molecules per side, found {len(lhs)}")
    binds, unbinds, changes = [], [], []
    for pos, (left, right) in enumerate(zip(lhs, rhs)):
        if left.type != right.type:
            raise ValueError(f"rule {text!r}: slot {pos + 1} holds {left.type} on the left and "
                             f"{right.type} on the right")
        if left.type not in by_name:
            raise ValueError(f"rule {text!r}: {left.type!r} is not a declared molecule type "
                             f"({sorted(by_name)})")
        if len(left.constraints) != len(right.constraints):
            raise ValueError(f"rule {text!r}: {left.type} must mention the same sites on both sides")
        mol = by_name[left.type]
        for ci, ((ls, lstate, llink), (rs, rstate, rlink)) in enumerate(
                zip(left.constraints, right.constraints)):
            if ls != rs:
                raise ValueError(f"rule {text!r}: {left.type} site {ls!r} becomes {rs!r}")
            if not mol.sites_named(ls):
                raise ValueError(f"rule {text!r}: molecule {mol.name} has no site {ls!r}")
            allowed: set[str] = set()
            for i in mol.sites_named(ls):
                allowed |= set(mol.sites[i].states)
            for st in (lstate, rstate):
                if st is not None and st not in allowed:
                    raise ValueError(f"rule {text!r}: site {mol.name}.{ls} has no internal state "
                                     f"{st!r} ({sorted(allowed) or 'none declared'})")
            if lstate != rstate:
                if rstate is None or lstate is None:
                    raise ValueError(f"rule {text!r}: site {mol.name}.{ls} must name its internal "
                                     "state on both sides to change it")
                changes.append((pos, ci, rstate))
            if llink != rlink:
                if llink is None and isinstance(rlink, tuple):
                    binds.append((rlink[1], pos, ci))
                elif isinstance(llink, tuple) and rlink is None:
                    unbinds.append((llink[1], pos, ci))
                else:
                    raise ValueError(f"rule {text!r}: site {mol.name}.{ls} cannot change from "
                                     f"{llink!r} to {rlink!r}")
    if bool(binds) + bool(unbinds) + bool(changes) != 1:
        raise ValueError(f"rule {text!r}: a rule must form exactly one bond, break exactly one bond "
                         "or change exactly one internal state")
    if changes:
        if len(changes) != 1:
            raise ValueError(f"rule {text!r}: only one internal state may change per rule")
        pos, ci, _new = changes[0]
        for pat in lhs:
            if not pat.constraints:
                raise ValueError(f"rule {text!r}: {pat.type} must name the site through which it "
                                 "meets its partner")
        contact = [ci if q == pos else 0 for q in range(len(lhs))]
        return Rule(label, "modify", k, lhs, contact, changes[0], text)
    half = binds or unbinds
    if len(half) != 2 or half[0][0] != half[1][0] or half[0][1] == half[1][1]:
        raise ValueError(f"rule {text!r}: a bond must be written with the same label on two "
                         "different molecules, e.g. A(a!1).B(b!1)")
    order = sorted(half, key=lambda h: h[1])
    return Rule(label, "bind" if binds else "unbind", k, lhs, [h[2] for h in order], None, text)


def parse_rules(texts, types) -> list[Rule]:
    """["'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ 0.5, 0.01"] -> one-way rules."""
    if not isinstance(texts, list) or not texts or not all(isinstance(t, str) for t in texts):
        raise ValueError("rules must be a non-empty list of rule strings such as "
                         "\"'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ 0.5, 0.01\"")
    out: list[Rule] = []
    labels: set[str] = set()
    for i, raw in enumerate(texts):
        text = " ".join(raw.split())
        m = _LABEL.match(text)
        label, body = (m.group(1), m.group(2)) if m else (f"r{i + 1}", text)
        if label in labels:
            raise ValueError(f"rule label {label!r} is used twice")
        labels.add(label)
        if "@" not in body:
            raise ValueError(f"rule {raw!r} needs its rate after '@'")
        expr, rate_text = body.rsplit("@", 1)
        rates = []
        for part in rate_text.split(","):
            try:
                rates.append(float(part))
            except ValueError:
                raise ValueError(f"rule {raw!r}: rates must be numbers, got {part.strip()!r}") from None
        if any(k < 0 for k in rates):
            raise ValueError(f"rule {raw!r}: rates must be non-negative")
        if "<->" in expr:
            if expr.count("<->") != 1 or len(rates) != 2:
                raise ValueError(f"reversible rule {raw!r} needs one '<->' and two rates 'k, k_rev'")
            left, right = expr.split("<->")
            lhs, rhs = _parse_side(left, raw), _parse_side(right, raw)
            out.append(_one_way(label, lhs, rhs, rates[0], raw, types))
            out.append(_one_way(label + "_rev", rhs, lhs, rates[1], raw, types))
            continue
        if expr.count("->") != 1 or len(rates) != 1:
            raise ValueError(f"rule {raw!r} needs exactly one '->' and one rate")
        left, right = expr.split("->")
        out.append(_one_way(label, _parse_side(left, raw), _parse_side(right, raw),
                            rates[0], raw, types))
    return out


# ============================================================================
# Reactive volume (paper, additional file 2)
# ============================================================================
def reactive_volume(d0: float, t_dist: float, t_ang: float,
                    orientation: str = "bound-only") -> float:
    """The volume in which a partner is geometrically compatible.

    Additional file 2: the spherical shell between d0 - t_dist and d0 + t_dist,
    (8 pi / 3)(3 d0^2 t_dist + t_dist^3), narrowed by one orientation factor
    (1 - cos t_ang)/2 per partner whose orientation is not fixed by other bonds.
    With ``orientation="sampled"`` both factors apply, giving the supplement's
    P * V_reactor = (2 pi / 3)(3 d0^2 t_dist + t_dist^3)(1 - cos t_ang)^2.
    """
    shell = (8.0 * math.pi / 3.0) * (3.0 * d0 ** 2 * t_dist + t_dist ** 3)
    if orientation == "sampled":
        shell *= ((1.0 - math.cos(math.radians(t_ang))) / 2.0) ** 2
    return shell


# ============================================================================
# Complexes: canonical form up to isomorphism
# ============================================================================
def _refine(members, kind, states, links, rounds):
    """Colour refinement: an isomorphism-invariant colour per elementary molecule."""
    colour = {}
    for i in members:
        t = kind[i]
        colour[i] = t.name + "|" + ",".join(sorted(
            f"{s.name}{'~' + s.states[states[i][si]] if s.states else ''}"
            f"{'+' if (i, si) in links else '-'}" for si, s in enumerate(t.sites)))
    for _ in range(rounds):
        new = {}
        for i in members:
            around = sorted(f"{kind[i].sites[si].name}:{colour[links[(i, si)][0]]}"
                            for si in range(len(kind[i].sites)) if (i, si) in links)
            new[i] = hashlib.blake2s((colour[i] + "|" + "|".join(around)).encode(),
                                     digest_size=8).hexdigest()
        if new == colour:
            break
        colour = new
    return colour


def canonical_text(members, kind, states, links) -> str:
    """Canonical text of the connected complex `members`; equal for isomorphic complexes.

    `kind[i]` is the MolType of elementary molecule i, `states[i][site]` the index
    of the site's internal state, and `links[(i, site)] = (j, site)` are the bonds.
    """
    members = sorted(members)
    colour = _refine(members, kind, states, links, len(members) + 1)
    lowest = min(colour[i] for i in members)
    best = None
    for root in (i for i in members if colour[i] == lowest):
        order, seen, site_order = [root], {root}, {}
        queue = [root]
        while queue:
            i = queue.pop(0)
            sites = sorted(range(len(kind[i].sites)),
                           key=lambda si, i=i: (kind[i].sites[si].name,
                                                0 if (i, si) in links else 1,
                                                colour[links[(i, si)][0]] if (i, si) in links else "",
                                                si))
            site_order[i] = sites
            for si in sites:
                partner = links.get((i, si))
                if partner is not None and partner[0] not in seen:
                    seen.add(partner[0])
                    order.append(partner[0])
                    queue.append(partner[0])
        bonds: dict[tuple, int] = {}
        parts = []
        for em in order:
            written = []
            for si in site_order[em]:
                site = kind[em].sites[si]
                state = f"{{{site.states[states[em][si]]}}}" if site.states else ""
                partner = links.get((em, si))
                if partner is None:
                    link = "."
                else:
                    if (em, si) not in bonds:
                        bonds[(em, si)] = bonds[partner] = len(bonds) // 2 + 1
                    link = str(bonds[(em, si)])
                written.append(f"{site.name}{state}[{link}]")
            parts.append(f"{kind[em].name}({','.join(written)})")
        text = ",".join(parts)
        if best is None or text < best:
            best = text
    return best


def species_id(text: str, kind, members) -> str:
    """The canonical text itself, or formula plus a hash of it for large complexes."""
    if len(text) <= 120:
        return text
    counts = Counter(kind[i].name for i in members)
    formula = "".join(f"{name}{n if n > 1 else ''}" for name, n in sorted(counts.items()))
    return f"{formula}#{hashlib.blake2s(text.encode(), digest_size=3).hexdigest()}"


# ============================================================================
# Built-in models
# ============================================================================
def build_model(p) -> dict:
    """The molecule types, rules and initial molecules of the chosen model."""
    if p.model == "custom":
        if not isinstance(p.init, dict) or not p.init or not all(
                isinstance(k, str) and isinstance(v, int) and not isinstance(v, bool) and v >= 0
                for k, v in p.init.items()):
            raise ValueError("model 'custom' needs init: molecule type -> non-negative count, "
                             f"e.g. {{'M': 50}}; got {p.init!r}")
        return {"molecule_types": p.molecule_types, "rules": p.rules, "init": dict(p.init)}
    if p.molecule_types or p.rules or p.init:
        raise ValueError("molecule_types, rules and init are only used with model='custom'")
    if p.model == "polymer":
        return {
            "molecule_types": [{
                "name": "M", "radius": 0.5, "mass": 1.0, "diffusion": p.diffusion,
                "sites": [{"name": "a", "dist": 1.0, "theta": 90.0, "phi": 0.0},
                          {"name": "b", "dist": 1.0, "theta": 90.0, "phi": p.bond_angle}],
            }],
            "rules": [f"'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ {p.k_on:g}, {p.k_off:g}"],
            "init": {"M": p.n_molecules},
        }
    if p.model == "dimerization":
        return {
            "molecule_types": [
                {"name": "A", "radius": 0.5, "mass": 1.0, "diffusion": p.diffusion,
                 "sites": [{"name": "a", "dist": 1.0, "theta": 90.0, "phi": 0.0}]},
                {"name": "B", "radius": 0.5, "mass": 1.0, "diffusion": p.diffusion,
                 "sites": [{"name": "b", "dist": 1.0, "theta": 90.0, "phi": 0.0}]},
            ],
            "rules": [f"'bind' A(a) + B(b) <-> A(a!1).B(b!1) @ {p.k_on:g}, {p.k_off:g}"],
            "init": {"A": p.n_molecules, "B": p.n_molecules},
        }
    # geometry of spass.geo, rules and rates of spass_w.bngl, prefactors of spass_w.in
    r, q = SCAFFOLD_RATES, SCAFFOLD_PREFACTORS
    phos, dephos = r["d1_on"] * q["modify2"], r["d1_off"] * q["modify1"]
    on2, off2 = r["d2_on"] * q["bind"] * p.scaffold_binding, r["d2_off"] * q["break"]
    on3, off3 = r["d3_on"] * q["bind"] * p.scaffold_binding, r["d3_off"] * q["break"]
    return {
        "molecule_types": [
            {"name": "A", "radius": 1.0, "mass": 1.0, "diffusion": p.diffusion,
             "sites": [{"name": "r", "dist": 1.0, "theta": 0.0, "phi": 0.0, "states": ["nP", "P"]},
                       {"name": "s", "dist": 1.0, "theta": 100.0, "phi": 0.0}]},
            {"name": "S", "radius": 8.0, "mass": 5.0, "diffusion": p.diffusion / 5.0,
             "sites": [{"name": "t", "dist": 8.0, "theta": 0.0, "phi": phi}
                       for phi in (0.0, 72.0, 144.0, 216.0)]},
        ],
        "rules": [
            f"'phosphorylate' A(r) + A(r~nP) -> A(r) + A(r~P) @ {phos:g}",
            f"'dephosphorylate' A(r~P) -> A(r~nP) @ {dephos:g}",
            f"'bind_nP' A(r~nP,s) + S(t) <-> A(r~nP,s!1).S(t!1) @ {on2:g}, {off2:g}",
            f"'bind_P' A(r~P,s) + S(t) <-> A(r~P,s!1).S(t!1) @ {on3:g}, {off3:g}",
        ],
        "init": {"A": p.n_molecules, "S": p.n_scaffolds},
    }


def bngl(spec: dict) -> str:
    """The unflattened model as a BioNetGen-style text (species and reaction rules)."""
    lines = ["begin species"]
    lines += [f"  {name}()  {n}" for name, n in spec["init"].items()]
    lines += ["end species", "", "begin reaction rules"]
    lines += [f"  {text}" for text in spec["rules"]]
    lines += ["end reaction rules"]
    return "\n".join(lines) + "\n"


# ============================================================================
# The reactor
# ============================================================================
class World:
    """Elementary molecules diffusing in a periodic box under a rule system."""

    def __init__(self, types, rules, init, p, rng):
        self.types, self.rules, self.p, self.rng = types, rules, p, rng
        self.L, self.dt = float(p.box), float(p.dt)
        self.well_mixed = bool(p.well_mixed)
        self.volume = self.L ** 3
        self.n_sites = max(len(t.sites) for t in types)
        by_name = {t.name: t for t in types}
        for name in init:
            if name not in by_name:
                raise ValueError(f"init names the unknown molecule type {name!r} "
                                 f"({sorted(by_name)})")
        kinds = [by_name[name] for name, n in init.items() for _ in range(n)]
        if not kinds:
            raise ValueError("the reactor is empty: init must place at least one molecule")
        self.n = len(kinds)
        self.kind = kinds
        self.tid = np.array([t.index for t in kinds], dtype=np.int32)
        self.pos = rng.random((self.n, 3)) * self.L
        self.vel = np.zeros((self.n, 3))
        self.start = self.pos.copy()
        self.unwrapped = self.pos.copy()
        self.state = np.zeros((self.n, self.n_sites), dtype=np.int8)
        self.partner = np.full((self.n, self.n_sites), -1, dtype=np.int32)
        self.partner_site = np.zeros((self.n, self.n_sites), dtype=np.int32)
        self.refract = np.full(self.n, -1.0)
        self.mass = np.array([t.mass for t in kinds])
        self.diff = np.array([t.diffusion for t in kinds])
        self.radius = np.array([t.radius for t in kinds])
        self.t_dist = float(p.distance_tolerance)
        self.t_ang = float(p.angle_tolerance)
        self.orientation = p.orientation
        self.p_orient = ((1.0 - math.cos(math.radians(self.t_ang))) / 2.0
                         if self.orientation == "sampled" else 1.0)
        self.time = 0.0
        self._dirty = True
        self._pairs = None
        self.events: dict[tuple, list] = {}
        self.tally = Tally()
        self.counts: Counter = Counter()
        self.species: dict[str, str] = {}
        self.mean_counts: Counter = Counter()
        self.mean_states: Counter = Counter()
        self.mean_refractory = 0.0
        self.mean_samples = 0
        self._prepare_rules()
        for i in range(self.n):
            self.counts[self._register([i])] += 1
        self.initial = Counter(self.counts)

    # -- rules -----------------------------------------------------------
    def _prepare_rules(self):
        reach = [0.0]
        for rule in self.rules:
            rule.options = [self._options(pat, rule.text) for pat in rule.patterns]
            if rule.order == 2:
                d0 = 0.0
                for q, pat in enumerate(rule.patterns):
                    mol = next(t for t in self.types if t.name == pat.type)
                    name = pat.constraints[rule.contact[q]][0]
                    d0 += mol.sites[mol.sites_named(name)[0]].dist
                rule.d0 = d0
                rule.v_react = reactive_volume(d0, self.t_dist, self.t_ang, self.orientation)
                if rule.kind != "unbind":
                    reach.append(d0 + self.t_dist)
        self.max_reaction_distance = max(reach)

    def _options(self, pat: Pattern, where: str) -> list[tuple]:
        """Every injective assignment of the pattern's constraints to site indices."""
        mol = next(t for t in self.types if t.name == pat.type)
        partial: list[tuple] = [()]
        for site, state, _link in pat.constraints:
            grown = []
            for chosen in partial:
                for si in mol.sites_named(site):
                    if si in chosen:
                        continue
                    if state is not None and state not in mol.sites[si].states:
                        continue
                    grown.append(chosen + (si,))
            partial = grown
        if not partial:
            raise ValueError(f"rule {where!r}: {pat.type} has no set of distinct sites matching "
                             "its pattern")
        return [(mol.index, assign) for assign in partial]

    def _mask(self, pat: Pattern, option: tuple) -> np.ndarray:
        tindex, assign = option
        mask = self.tid == tindex
        mol = self.types[tindex]
        for (_site, state, link), si in zip(pat.constraints, assign):
            if state is not None:
                mask = mask & (self.state[:, si] == mol.sites[si].states.index(state))
            if link is None:
                mask = mask & (self.partner[:, si] < 0)
            elif link == "+" or isinstance(link, tuple):
                mask = mask & (self.partner[:, si] >= 0)
        return mask

    # -- geometry --------------------------------------------------------
    def _delta(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        d = a - b
        return d - self.L * np.round(d / self.L)

    def _links(self) -> dict:
        out = {}
        rows, cols = np.nonzero(self.partner >= 0)
        for i, si in zip(rows.tolist(), cols.tolist()):
            out[(i, si)] = (int(self.partner[i, si]), int(self.partner_site[i, si]))
        return out

    def component(self, seed: int) -> list[int]:
        seen, stack = {seed}, [seed]
        while stack:
            i = stack.pop()
            for si in range(len(self.kind[i].sites)):
                j = int(self.partner[i, si])
                if j >= 0 and j not in seen:
                    seen.add(j)
                    stack.append(j)
        return sorted(seen)

    def _register(self, members: list[int]) -> str:
        text = canonical_text(members, self.kind, self.state, self._links())
        sid = species_id(text, self.kind, members)
        self.species.setdefault(sid, text)
        return sid

    def geometry_ok(self, i: int, si: int, j: int) -> bool:
        """SRSim's testSiteGeo: every bond of i must keep its ideal angle to the direction i->j."""
        bound = [s for s in range(len(self.kind[i].sites))
                 if s != si and self.partner[i, s] >= 0]
        if not bound:
            return self.p_orient >= 1.0 or float(self.rng.random()) < self.p_orient
        direction = self._delta(self.pos[j], self.pos[i])
        ideal = self.kind[i].ideal[si]
        for s in bound:
            other = self._delta(self.pos[int(self.partner[i, s])], self.pos[i])
            norm = float(np.linalg.norm(direction) * np.linalg.norm(other)) + 1e-12
            cos = float(np.dot(direction, other)) / norm
            angle = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
            if abs(angle - ideal[s]) > self.t_ang:
                return False
        return True

    # -- forces and motion -----------------------------------------------
    def _topology(self):
        if not self._dirty:
            return
        links = self._links()
        bonds, angles = [], []
        for (i, si), (j, sj) in links.items():
            if i < j:
                bonds.append((i, j, self.kind[i].sites[si].dist + self.kind[j].sites[sj].dist))
        for i in range(self.n):
            bound = [s for s in range(len(self.kind[i].sites)) if self.partner[i, s] >= 0]
            for a in range(len(bound)):
                for b in range(a + 1, len(bound)):
                    angles.append((i, int(self.partner[i, bound[a]]),
                                   int(self.partner[i, bound[b]]),
                                   self.kind[i].ideal[bound[a]][bound[b]]))
        self._bonds = (np.array([b[0] for b in bonds], dtype=np.int64),
                       np.array([b[1] for b in bonds], dtype=np.int64),
                       np.array([b[2] for b in bonds], dtype=float))
        self._angles = (np.array([a[0] for a in angles], dtype=np.int64),
                        np.array([a[1] for a in angles], dtype=np.int64),
                        np.array([a[2] for a in angles], dtype=np.int64),
                        np.radians(np.array([a[3] for a in angles], dtype=float)))
        self._bonded_pairs = [(min(i, j), max(i, j)) for i, j, _ in bonds]
        self._dirty = False

    def forces(self) -> np.ndarray:
        self._topology()
        f = np.zeros((self.n, 3))
        bi, bj, b0 = self._bonds
        if bi.size:                                     # E_d = K_d (d - d_ij - d_kl)^2
            d = self._delta(self.pos[bj], self.pos[bi])
            r = np.linalg.norm(d, axis=1) + 1e-12
            pull = (2.0 * self.p.k_bond * (r - b0) / r)[:, None] * d
            np.add.at(f, bi, pull)
            np.add.at(f, bj, -pull)
        ci, cj, ck, a0 = self._angles
        if ci.size:                                     # E_a = K_a (alpha - alpha_ijk)^2
            u = self._delta(self.pos[cj], self.pos[ci])
            v = self._delta(self.pos[ck], self.pos[ci])
            ru = np.linalg.norm(u, axis=1) + 1e-12
            rv = np.linalg.norm(v, axis=1) + 1e-12
            cos = np.clip(np.sum(u * v, axis=1) / (ru * rv), -1.0, 1.0)
            alpha = np.arccos(cos)
            sin = np.sqrt(np.maximum(1.0 - cos ** 2, 1e-6))
            scale = 2.0 * self.p.k_angle * (alpha - a0) / sin
            fj = scale[:, None] * (v / (ru * rv)[:, None] - (cos / ru ** 2)[:, None] * u)
            fk = scale[:, None] * (u / (ru * rv)[:, None] - (cos / rv ** 2)[:, None] * v)
            np.add.at(f, cj, fj)
            np.add.at(f, ck, fk)
            np.add.at(f, ci, -(fj + fk))
        if self.p.k_repulsion > 0:                      # E_s = A [1 + cos(pi r / r_c)]
            d = self._delta(self.pos[:, None, :], self.pos[None, :, :])
            r = np.linalg.norm(d, axis=2)
            rc = self.radius[:, None] + self.radius[None, :]
            np.fill_diagonal(r, np.inf)
            for i, j in self._bonded_pairs:
                r[i, j] = r[j, i] = np.inf
            close = r < rc
            if close.any():
                mag = np.zeros_like(r)
                mag[close] = (self.p.k_repulsion * math.pi / rc[close]
                              * np.sin(math.pi * r[close] / rc[close]))
                safe = np.where(close, r, 1.0)
                f += np.sum((mag / safe)[:, :, None] * d, axis=1)
        return f

    def integrate(self):
        if self.well_mixed:              # the non-spatial control: positions play no role
            return
        f = self.forces()
        if self.p.integrator == "brownian":     # overdamped limit of the Langevin equation
            noise = self.rng.standard_normal((self.n, 3))
            step = ((self.diff / self.p.temperature)[:, None] * f * self.dt
                    + np.sqrt(2.0 * self.diff * self.dt)[:, None] * noise)
        else:                                   # F = F^S - gamma v + sqrt(2 kB T gamma) xi
            gamma = self.p.temperature / np.maximum(self.diff, 1e-12)
            noise = self.rng.standard_normal((self.n, 3))
            total = (f - gamma[:, None] * self.vel
                     + np.sqrt(2.0 * self.p.temperature * gamma / self.dt)[:, None] * noise)
            self.vel += self.dt * total / self.mass[:, None]
            step = self.dt * self.vel
        self.pos = self.pos + step
        self.unwrapped = self.unwrapped + step
        self.pos %= self.L

    # -- the run ----------------------------------------------------------
    def run(self):
        """Run the steps, yielding after each one.

        `steps=0` runs on until the caller stops reading; the rate samples
        then start at once, there being no half-way point to wait for.
        """
        every = max(1, self.p.steps // 200)
        half = self.p.steps // 2
        for step in (t - 1 for t in ticks(self.p.steps)):
            self.integrate()
            self.time = (step + 1) * self.dt
            self.react()
            if step >= half and step % every == 0:
                self._sample()
            yield step + 1
        if not self.mean_samples:
            self._sample()

    def frame(self) -> Frame:
        return Frame(t=float(self.time), state={s: float(n) for s, n in self.counts.items()},
                     fired=self.tally.flush())

    def _sample(self):
        self.mean_samples += 1
        self.mean_counts.update(self.counts)
        self.mean_states.update(self.state_counts())
        self.mean_refractory += float(np.count_nonzero(self.refract > self.time))

    def state_counts(self) -> Counter:
        out: Counter = Counter()
        for t in self.types:
            rows = np.flatnonzero(self.tid == t.index)
            for si, site in enumerate(t.sites):
                for k, name in enumerate(site.states):
                    n = int(np.count_nonzero(self.state[rows, si] == k))
                    if n:
                        out[f"{t.name}.{site.name}~{name}"] += n
        return out

    # -- reactions --------------------------------------------------------
    def react(self):
        self._pairs = None
        used = np.zeros(self.n, dtype=bool)
        masks = {}
        for rule in self.rules:
            for q, options in enumerate(rule.options):
                for option in options:
                    masks[(id(rule), q, option)] = self._mask(rule.patterns[q], option)
        for rule in self.rules:
            if rule.k <= 0:
                continue
            probability = 1.0 - math.exp(-rule.k * self.dt)
            if rule.kind == "unbind":
                self._fire_unbind(rule, masks, used, probability)
            elif rule.order == 1:
                self._fire_unimolecular(rule, masks, used, probability)
            elif self.well_mixed:
                self._fire_well_mixed(rule, masks, used, probability)
            else:
                self._fire_bimolecular(rule, masks, used, probability)

    def _fire_unimolecular(self, rule, masks, used, probability):
        for option in rule.options[0]:
            index = np.flatnonzero(masks[(id(rule), 0, option)] & ~used)
            if not index.size:
                continue
            for i in index[self.rng.random(index.size) < probability]:
                i = int(i)
                if used[i]:
                    continue
                self.apply(rule, [(i, option[1][rule.contact[0]])])
                used[i] = True

    def _fire_unbind(self, rule, masks, used, probability):
        links = [(i, si, j, sj) for (i, si), (j, sj) in self._links().items() if i < j]
        for n in self.rng.permutation(len(links)) if links else ():
            i, si, j, sj = links[int(n)]
            if used[i] or used[j] or float(self.rng.random()) >= probability:
                continue
            for order in ((0, 1), (1, 0)):
                first, fs = (i, si) if order[0] == 0 else (j, sj)
                second, ss = (j, sj) if order[0] == 0 else (i, si)
                hit = False
                for oa in rule.options[order[0]]:
                    if not masks[(id(rule), order[0], oa)][first]:
                        continue
                    if oa[1][rule.contact[order[0]]] != fs:
                        continue
                    for ob in rule.options[order[1]]:
                        if not masks[(id(rule), order[1], ob)][second]:
                            continue
                        if ob[1][rule.contact[order[1]]] != ss:
                            continue
                        sites = [None, None]
                        sites[order[0]] = (first, fs)
                        sites[order[1]] = (second, ss)
                        self.apply(rule, sites)
                        used[i] = used[j] = hit = True
                        break
                    if hit:
                        break
                if hit:
                    break

    def _candidate_pairs(self):
        d = self._delta(self.pos[:, None, :], self.pos[None, :, :])
        r2 = np.sum(d * d, axis=2)
        upper = np.triu(np.ones((self.n, self.n), dtype=bool), 1)
        pi, pj = np.nonzero(upper & (r2 < self.max_reaction_distance ** 2))
        return pi, pj, r2[pi, pj]

    def _fire_bimolecular(self, rule, masks, used, probability):
        if self._pairs is None:
            self._pairs = self._candidate_pairs()
        pi, pj, r2 = self._pairs
        if not pi.size:
            return
        low = max(rule.d0 - self.t_dist, 0.0) ** 2
        high = (rule.d0 + self.t_dist) ** 2
        window = (r2 >= low) & (r2 <= high)
        if not window.any():
            return
        for order in ((0, 1), (1, 0)):
            for oa in rule.options[order[0]]:
                ma = masks[(id(rule), order[0], oa)]
                for ob in rule.options[order[1]]:
                    mb = masks[(id(rule), order[1], ob)]
                    sel = np.flatnonzero(window & ma[pi] & mb[pj])
                    if not sel.size:
                        continue
                    for n in sel[self.rng.random(sel.size) < probability]:
                        i, j = int(pi[n]), int(pj[n])
                        if used[i] or used[j]:
                            continue
                        sa = oa[1][rule.contact[order[0]]]
                        sb = ob[1][rule.contact[order[1]]]
                        if rule.kind == "bind" and (self.refract[i] > self.time
                                                    or self.refract[j] > self.time):
                            continue
                        if not (self.geometry_ok(i, sa, j) and self.geometry_ok(j, sb, i)):
                            continue
                        sites = [None, None]
                        sites[order[0]] = (i, sa)
                        sites[order[1]] = (j, sb)
                        self.apply(rule, sites)
                        used[i] = used[j] = True

    def _fire_well_mixed(self, rule, masks, used, probability):
        """The non-spatial control: any pair is compatible with probability V_react / V_reactor."""
        chance = min(probability * rule.v_react / self.volume, 1.0)
        for oa in rule.options[0]:
            first = np.flatnonzero(masks[(id(rule), 0, oa)])
            if not first.size:
                continue
            for ob in rule.options[1]:
                second = np.flatnonzero(masks[(id(rule), 1, ob)])
                if not second.size:
                    continue
                draws = int(self.rng.binomial(int(first.size) * int(second.size), chance))
                for _ in range(draws):
                    i = int(first[self.rng.integers(first.size)])
                    j = int(second[self.rng.integers(second.size)])
                    if i == j or used[i] or used[j]:
                        continue
                    if rule.kind == "bind" and (self.refract[i] > self.time
                                                or self.refract[j] > self.time):
                        continue
                    self.apply(rule, [(i, oa[1][rule.contact[0]]), (j, ob[1][rule.contact[1]])])
                    used[i] = used[j] = True

    # -- applying a rule ---------------------------------------------------
    def _complexes_of(self, sites) -> list[list[int]]:
        out, seen = [], set()
        for i, _si in sites:
            if i in seen:
                continue
            members = self.component(i)
            seen.update(members)
            out.append(members)
        return out

    def apply(self, rule, sites):
        before = [self._register(m) for m in self._complexes_of(sites)]
        if rule.kind == "bind":
            (i, si), (j, sj) = sites
            self.partner[i, si], self.partner_site[i, si] = j, sj
            self.partner[j, sj], self.partner_site[j, sj] = i, si
            self._dirty = True
        elif rule.kind == "unbind":
            (i, si), (j, sj) = sites
            self.partner[i, si] = self.partner[j, sj] = -1
            self.refract[i] = self.refract[j] = self.time + self.p.refractory_time
            self._dirty = True
        else:
            pos, _ci, new = rule.change
            i, si = sites[pos]
            self.state[i, si] = self.kind[i].sites[si].states.index(new)
        after = [self._register(m) for m in self._complexes_of(sites)]
        for sid in before:
            self.counts[sid] -= 1
            if self.counts[sid] <= 0:
                del self.counts[sid]
        for sid in after:
            self.counts[sid] += 1
        key = (tuple(sorted(before)), tuple(sorted(after)))
        entry = self.events.setdefault(key, [Counter(before), Counter(after), 0, []])
        entry[2] += 1
        if rule.label not in entry[3]:
            entry[3].append(rule.label)
        self.tally.add(before, after)
        self._pairs = None


# ============================================================================
# Analysis and network assembly
# ============================================================================
def _analysis(world: World) -> dict:
    links = world._links()
    sizes: Counter = Counter()
    closed = 0
    seen: set[int] = set()
    for i in range(world.n):
        if i in seen:
            continue
        members = world.component(i)
        seen.update(members)
        sizes[len(members)] += 1
        bonds = sum(1 for (a, _s) in links if a in members) // 2
        if bonds >= len(members):
            closed += 1
    samples = max(1, world.mean_samples)
    return {
        "steps": int(world.p.steps),
        "simulated_time": float(world.p.steps * world.dt),
        "events": int(sum(e[2] for e in world.events.values())),
        "events_by_rule": {r.label: int(sum(e[2] for e in world.events.values() if r.label in e[3]))
                           for r in world.rules},
        "largest_complex": int(max(sizes) if sizes else 0),
        "size_histogram": {str(k): int(v) for k, v in sorted(sizes.items())},
        "closed_complexes": int(closed),
        "bonds": int(len(links) // 2),
        "final_counts": {s: int(n) for s, n in sorted(world.counts.items())},
        "mean_counts": {s: float(n) / samples for s, n in sorted(world.mean_counts.items())},
        "final_state_counts": {s: int(n) for s, n in sorted(world.state_counts().items())},
        "mean_state_counts": {s: float(n) / samples for s, n in sorted(world.mean_states.items())},
        "mean_refractory": float(world.mean_refractory) / samples,
        "mean_squared_displacement": float(np.mean(np.sum(
            (world.unwrapped - world.start) ** 2, axis=1))),
        "expected_msd_6Dt": float(6.0 * float(np.mean(world.diff)) * world.p.steps * world.dt),
    }


def _space(world: World, types) -> dict:
    belongs: dict[int, str] = {}
    seen: set[int] = set()
    for i in range(world.n):
        if i not in seen:
            members = world.component(i)
            seen.update(members)
            sid = world._register(members)
            for m in members:
                belongs[m] = sid
    return {
        "kind": "continuous",
        "dimensions": 3,
        "units": "arbitrary (as in the paper's examples)",
        "box": [[0.0, float(world.L)]] * 3,
        "boundary": ("none: the well-mixed control ignores positions" if world.well_mixed
                     else "periodic (minimum image)"),
        "timestep": float(world.dt),
        "steps": int(world.p.steps),
        "temperature": float(world.p.temperature),
        "integrator": "well-mixed" if world.well_mixed else world.p.integrator,
        "diffusion": {t.name: t.diffusion for t in types},
        "radius": {t.name: t.radius for t in types},
        "mass": {t.name: t.mass for t in types},
        "sites": {t.name: [{"name": s.name, "dist": s.dist, "theta": s.theta, "phi": s.phi,
                            "states": list(s.states)} for s in t.sites] for t in types},
        "ideal_angles": {t.name: [[round(v, 6) for v in row] for row in t.ideal] for t in types},
        "tolerances": {"distance": world.t_dist, "angle_degrees": world.t_ang,
                       "orientation": world.orientation},
        "refractory_time": float(world.p.refractory_time),
        "reaction_criterion":
            "a bimolecular rule fires between two elementary molecules only if their distance is "
            "within distance_tolerance of the ideal bond length d_ij + d_kl and, for every bond "
            "either of them already carries, the angle between the direction to the partner and "
            "that bond is within angle_tolerance of the ideal angle; an unbound molecule has no "
            "orientation and passes (orientation='sampled': with probability (1 - cos t_ang)/2). "
            "A compatible pair then reacts with probability 1 - exp(-k dt).",
        "particles": [{"species": belongs[i], "type": world.kind[i].name,
                       "position": [round(float(v), 4) for v in world.pos[i]]}
                      for i in range(world.n)],
    }


def _law(p) -> dict:
    return {
        "name": "SRSim coarse-grained molecular dynamics (LAMMPS)",
        "state": "each elementary molecule is a sphere with a position, a velocity, a mass, a "
                 "radius and components at fixed polar coordinates",
        "potentials": {
            "bond": "E_d = K_d (d - d_ij - d_kl)^2",
            "angle": "E_a = K_a (alpha - alpha_ijk)^2 for every pair of bonds of one molecule",
            "repulsion": "E_s = A [1 + cos(pi r / r_c)] for r < r_c = r_i + r_j (soft sphere)",
        },
        "motion": "F_i = F_i^S - gamma_0 v_i + sqrt(2 k_B T gamma_0) xi_i(t) with "
                  "gamma_0 = k_B T / D (Langevin), or its overdamped limit "
                  "dx = D F dt / k_B T + sqrt(2 D dt) xi (Brownian); k_B = 1",
        "constants": {"K_d": float(p.k_bond), "K_a": float(p.k_angle), "A": float(p.k_repulsion),
                      "temperature": float(p.temperature)},
        "time_step_rule": "the paper's rule of thumb: sqrt(6 D dt) should be about a tenth of a "
                          "particle diameter",
    }


def _kinetics(world: World) -> dict:
    return {
        "k1_micro": "k1_mic = k1_mac for monomolecular rules (additional file 2)",
        "k2_micro": "k2_mic = k2_mac / (V_reactor P) with P = V_react / V_reactor; both partners "
                    "unoriented gives V_react = (2 pi / 3)(3 d0^2 t_dist + t_dist^3)"
                    "(1 - cos t_ang)^2, orientations fixed by other bonds give the shell "
                    "(8 pi / 3)(3 d0^2 t_dist + t_dist^3)",
        "reaction_rate": "the k of every observed reaction is the macroscopic constant of the rule "
                         "that fired it: k_mic * V_react for a bimolecular rule (volume per time, "
                         "so that k [A][B] is the well-stirred rate) and k_mic otherwise",
        "reactor_volume": float(world.volume),
        "rules": {r.label: {"k_micro": r.k, "order": r.order, "kind": r.kind,
                            "ideal_bond_length": r.d0, "reactive_volume": r.v_react,
                            "k_macro": _k_macro(r)} for r in world.rules},
    }


def _k_macro(rule: Rule) -> float:
    return rule.k * (rule.v_react if rule.order == 2 and rule.kind != "unbind" else 1.0)


def evolve(p, rng):
    """The spatial run: a frame per molecular-dynamics step (frame 0 is the initial state)."""
    if p.model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}, got {p.model!r}")
    if p.integrator not in INTEGRATORS:
        raise ValueError(f"integrator must be one of {INTEGRATORS}, got {p.integrator!r}")
    if p.orientation not in ORIENTATIONS:
        raise ValueError(f"orientation must be one of {ORIENTATIONS}, got {p.orientation!r}")
    spec = build_model(p)
    types = parse_types(spec["molecule_types"])
    rules = parse_rules(spec["rules"], types)
    for rule in rules:
        rule.k *= float(p.rate_scale)
    world = World(types, rules, spec["init"], p, rng)
    if world.max_reaction_distance > world.L / 2:
        raise ValueError(f"box = {p.box} is too small: the reaction distance is "
                         f"{world.max_reaction_distance:.2f} and the minimum image convention "
                         "needs a box of more than twice that")
    yield world.frame()
    for _ in world.run():
        yield world.frame()

    species: dict[str, str] = {}
    reactions, labels = [], []
    for (before, after), (left, right, count, rule_labels) in world.events.items():
        k = sum(_k_macro(next(r for r in world.rules if r.label == label)) for label in rule_labels)
        reactions.append(Reaction(dict(left), dict(right),
                                  {"law": "mass-action", "k": float(k),
                                   "rules": ",".join(rule_labels)}, int(count)))
        labels.append(list(rule_labels))
        for sid in (*before, *after):
            species[sid] = world.species[sid]
    space = _space(world, types)
    for sid in (*world.initial, *world.counts):
        species[sid] = world.species[sid]

    conservation = []
    for t in types:
        vector = {sid: text.count(f"{t.name}(") for sid, text in species.items()}
        if any(vector.values()):
            conservation.append({"name": f"elementary molecules of type {t.name}",
                                 "vector": vector})
    return Network(
        species=[Species(sid, structure=text) for sid, text in species.items()],
        reactions=reactions,
        status="observed",
        initial_state={sid: float(n) for sid, n in world.initial.items()},
        extras={
            "model": p.model,
            "molecule_types": spec["molecule_types"],
            "rules": list(spec["rules"]),
            "init": dict(spec["init"]),
            "bngl": bngl(spec),
            "reaction_rules": labels,
            "space": space,
            "interaction_law": _law(p),
            "kinetics": _kinetics(world),
            "conservation": conservation,
            "species_encoding":
                "id = structure = canonical site-graph text (colour-refined breadth-first order, "
                "sites grouped by name, bonds numbered by appearance), so isomorphic complexes are "
                "one species; a complex whose text exceeds 120 characters gets formula#hash",
            "analysis": _analysis(world),
        },
    )
