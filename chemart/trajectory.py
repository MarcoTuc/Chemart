"""The record of a run: a network and its frames over time.

Every way of running a chemistry returns a `Trajectory`: integrating a network
(`chemart.simulate.ode`), sampling it stochastically (`chemart.simulate.ssa`),
or evolving a Turing gas (`chemart.evolve`). A trajectory is a list of frames;
each frame holds the time, the amount of every present species, the reactions
that fired since the previous frame, and any quantities the chemistry itself
reports. Like `Network`, it is plain JSON and round-trips exactly.

Frame fields:

- ``t``            the time, in the trajectory's ``clock`` unit
- ``state``        {species: amount}; zeros are left out. Amounts are per unit
                   volume, so they are counts when volume and Avogadro's number
                   are 1 (the stochastic simulator and the gases)
- ``fired``        [[reactant ids], [product ids], count] for each reaction that
                   fired since the previous frame (empty for ODE runs)
- ``observables``  chemistry-specific scalars, e.g. BFF's high-order entropy
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Any

from chemart.network import Network, Reaction, Species

METHODS = ("ode", "ssa", "evolve")


@dataclass
class Frame:
    t: float
    state: dict[str, float]
    fired: list[list] = field(default_factory=list)
    observables: dict[str, Any] = field(default_factory=dict)


@dataclass
class Trajectory:
    network: Network
    frames: list[Frame]
    method: str
    clock: str = "time"
    #: Provenance: seed, parameters, solver options, the rate and initial-state specs.
    settings: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}, got {self.method!r}")

    # ------------------------------------------------------------------
    def times(self) -> list[float]:
        return [f.t for f in self.frames]

    def series(self, name: str) -> list[Any]:
        """An observable over time; None in frames that do not report it."""
        if not any(name in f.observables for f in self.frames):
            known = sorted({k for f in self.frames for k in f.observables})
            raise KeyError(f"no observable {name!r}; this trajectory has {known}")
        return [f.observables.get(name) for f in self.frames]

    def array(self, species: list[str] | None = None):
        """Dense (species_ids, times, X) with X[time, species]; absent species are 0.

        By default every species that is ever present, in order of appearance.
        """
        import numpy as np

        if species is None:
            species = list(dict.fromkeys(s for f in self.frames for s in f.state))
        X = np.array([[f.state.get(s, 0.0) for s in species] for f in self.frames], dtype=float)
        return species, np.array(self.times(), dtype=float), X.reshape(len(self.frames), len(species))

    def window(self, i: int, width: int = 1) -> Network:
        """The observed network of the reactions fired in frames i-width+1 .. i.

        `width=None` takes everything up to frame i (the cumulative network).
        """
        if i < 0:
            i += len(self.frames)
        start = 0 if width is None else max(0, i - width + 1)
        fired: dict[tuple, list] = {}
        for f in self.frames[start:i + 1]:
            for lhs, rhs, n in f.fired:
                key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
                fired.setdefault(key, [lhs, rhs, 0])[2] += n
        names = dict.fromkeys(s for lhs, rhs, _ in fired.values() for s in (*lhs, *rhs))
        known = {s.id: s for s in self.network.species}
        return Network(
            species=[known.get(s, Species(s)) for s in names],
            reactions=[Reaction.of(lhs, rhs, count=n) for lhs, rhs, n in fired.values()],
            status="observed",
            chemistry=self.network.chemistry,
            params=self.network.params,
            seed=self.network.seed,
        )

    def turnover(self) -> list[float]:
        """A clock shared by every gas: cumulative reactions fired per molecule.

        Collisions, epochs and iterations cannot be compared across chemistries;
        how many times the population has turned over can.
        """
        out, total = [], 0
        for f in self.frames:
            total += sum(n for _, _, n in f.fired)
            size = sum(f.state.values())
            out.append(total / size if size else 0.0)
        return out

    # ------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return {
            "network": self.network.to_dict(),
            "frames": [asdict(f) for f in self.frames],
            "method": self.method,
            "clock": self.clock,
            "settings": self.settings,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Trajectory:
        if not isinstance(data, dict):
            raise ValueError(f"a trajectory must be a JSON object, got {type(data).__name__}")
        try:
            return cls(
                network=Network.from_dict(data["network"]),
                frames=[Frame(**f) for f in data["frames"]],
                method=data["method"],
                clock=data.get("clock", "time"),
                settings=data.get("settings", {}),
            )
        except (TypeError, KeyError) as err:
            raise ValueError(f"malformed trajectory: {err}") from None
