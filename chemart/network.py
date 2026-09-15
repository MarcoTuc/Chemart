"""The network record every chemistry returns.

A `Network` is a materialised chemical reaction network: species, reactions
with separate reactant and product stoichiometry, optional rate laws, and the
provenance needed to reproduce it (chemistry id, parameters, seed, status).
Everything in it is plain JSON data, so `Network.from_dict(net.to_dict())`
round-trips exactly.

Status values:

- ``complete``  the whole network the chemistry defines (or its full closure)
- ``truncated`` a closure that was cut off by a size budget
- ``observed``  the reactions that actually fired in a simulation; each
  reaction then carries a firing ``count``
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable

from chemart import kinetics

STATUSES = ("complete", "truncated", "observed")

#: `outflow` value for the book's standard flow reactor: a non-selective
#: dilution flux Phi(t) that keeps the total concentration constant
#: (x_i' = production_i - x_i * Phi / sum_k x_k).
CONSTANT_TOTAL = "constant-total"

#: `extras` keys with a fixed meaning; any other key is free-form.
RESERVED_EXTRAS = (
    "space", "compartments", "energies", "conservation", "analysis",
    "interaction_law",
)

_EXTRAS_PROVIDE = {
    "space": "space",
    "compartments": "compartments",
    "energies": "energies",
    "conservation": "mass-conservation",
}


@dataclass
class Species:
    id: str
    structure: str | None = None


@dataclass
class Reaction:
    reactants: dict[str, int]
    products: dict[str, int]
    rate: dict[str, Any] | None = None
    count: int | None = None

    @classmethod
    def of(cls, reactants: Iterable[str], products: Iterable[str],
           rate: dict | None = None, count: int | None = None) -> Reaction:
        """Build a reaction from species-id sequences: ``of(["A", "A"], ["B"])``."""
        return cls(dict(Counter(reactants)), dict(Counter(products)), rate, count)

    @property
    def catalysts(self) -> dict[str, int]:
        """Species present on both sides, with the multiplicity that survives."""
        return {
            s: min(n, self.products[s])
            for s, n in self.reactants.items() if s in self.products
        }

    def to_text(self) -> str:
        text = f"{_side(self.reactants)} -> {_side(self.products)}"
        if self.rate:
            args = " ".join(f"{k}={v}" for k, v in self.rate.items() if k != "law")
            text += f"  [{self.rate['law']}{' ' + args if args else ''}]"
        if self.count is not None:
            text += f"  (x{self.count})"
        return text


def _side(stoich: dict[str, int]) -> str:
    if not stoich:
        return "∅"
    return " + ".join(s if n == 1 else f"{n} {s}" for s, n in stoich.items())


@dataclass
class Network:
    species: list[Species]
    reactions: list[Reaction]
    status: str = "complete"
    initial_state: dict[str, float] | None = None
    inflow: dict[str, float] | None = None
    outflow: dict[str, float] | float | str | None = None   # or CONSTANT_TOTAL
    extras: dict[str, Any] = field(default_factory=dict)
    # Provenance; filled in by chemart.generate_network.
    chemistry: str = ""
    params: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None

    def __post_init__(self) -> None:
        problems = self.problems()
        if problems:
            more = f" (+{len(problems) - 5} more)" if len(problems) > 5 else ""
            raise ValueError("invalid network: " + "; ".join(problems[:5]) + more)

    # ------------------------------------------------------------------
    def problems(self) -> list[str]:
        out: list[str] = []
        if self.status not in STATUSES:
            out.append(f"status {self.status!r} not in {STATUSES}")
        ids = [s.id for s in self.species]
        known = set(ids)
        if len(known) != len(ids):
            dup = next(i for i, n in Counter(ids).items() if n > 1)
            out.append(f"duplicate species id {dup!r}")
        for s in self.species:
            if not isinstance(s.id, str) or not s.id:
                out.append(f"species id must be a non-empty string, got {s.id!r}")
        for i, r in enumerate(self.reactions):
            if not r.reactants and not r.products:
                out.append(f"reaction {i} is empty on both sides")
            for side in (r.reactants, r.products):
                for s, n in side.items():
                    if s not in known:
                        out.append(f"reaction {i} uses unknown species {s!r}")
                    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
                        out.append(f"reaction {i}: stoichiometry of {s!r} must be a positive int, got {n!r}")
            if r.rate is not None:
                problem = kinetics.rate_problem(r.rate)
                if problem:
                    out.append(f"reaction {i}: {problem}")
            if r.count is not None and (not isinstance(r.count, int) or r.count < 0):
                out.append(f"reaction {i}: count must be a non-negative int")
        for name in ("initial_state", "inflow"):
            for s in getattr(self, name) or {}:
                if s not in known:
                    out.append(f"{name} uses unknown species {s!r}")
        if isinstance(self.outflow, dict):
            for s in self.outflow:
                if s not in known:
                    out.append(f"outflow uses unknown species {s!r}")
        elif isinstance(self.outflow, str) and self.outflow != CONSTANT_TOTAL:
            out.append(f"outflow must be a dict, a number or {CONSTANT_TOTAL!r}, got {self.outflow!r}")
        return out

    # ------------------------------------------------------------------
    @property
    def provides(self) -> list[str]:
        """Capability tags derived from what this network actually contains."""
        tags = {"topology"}
        if self.reactions:
            tags.add("stoichiometry")
        if any(r.catalysts for r in self.reactions):
            tags.add("catalysts")
        if any(r.rate for r in self.reactions):
            tags.add("rate-constants")
        if any(r.rate and r.rate["law"] != "mass-action" for r in self.reactions):
            tags.add("rate-law")
        if self.inflow or self.outflow:
            tags.add("flow")
        if self.initial_state:
            tags.add("initial-state")
        tags.update(tag for key, tag in _EXTRAS_PROVIDE.items() if self.extras.get(key))
        return sorted(tags)

    def matrices(self):
        """Return ``(species_ids, R, P)`` as sparse int matrices (species x reactions).

        The net stoichiometric matrix is ``S = P - R``.
        """
        import numpy as np
        from scipy.sparse import coo_array

        ids = [s.id for s in self.species]
        index = {s: i for i, s in enumerate(ids)}
        shape = (len(ids), len(self.reactions))

        def build(side: str):
            rows, cols, vals = [], [], []
            for j, r in enumerate(self.reactions):
                for s, n in getattr(r, side).items():
                    rows.append(index[s])
                    cols.append(j)
                    vals.append(n)
            return coo_array(
                (np.array(vals, dtype=np.int64), (np.array(rows, dtype=np.int64), np.array(cols, dtype=np.int64))),
                shape=shape,
            ).tocsc()

        return ids, build("reactants"), build("products")

    def summary(self) -> str:
        head = self.chemistry or "network"
        lines = [
            f"{head}: {len(self.species)} species, {len(self.reactions)} reactions, "
            f"status={self.status}",
            f"provides: {', '.join(self.provides)}",
        ]
        if self.seed is not None:
            lines.append(f"seed: {self.seed}")
        if not self.reactions:
            law = "; see extras['interaction_law']" if "interaction_law" in self.extras else ""
            lines.append(f"no reactions{law}")
        extra = sorted(self.extras)
        if extra:
            lines.append(f"extras: {', '.join(extra)}")
        return "\n".join(lines)

    def to_text(self) -> str:
        return "\n".join(r.to_text() for r in self.reactions)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Network:
        data = dict(data)
        data["species"] = [Species(**s) for s in data["species"]]
        data["reactions"] = [Reaction(**r) for r in data["reactions"]]
        return cls(**data)
