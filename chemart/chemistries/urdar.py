"""Urdar (Gerlee & Lundh 2010): cross-feeding digital organisms. Catalog id: urdar.

Organisms are elementary cellular automaton rules (256 species, Wolfram
numbering). Metabolites are binary strings of length L in a well-stirred pool
of N_R strings. One update (Evolution 64:2716, "Implementation"):

1. every organism picks a random string of the pool and applies one CA step
   of its rule to it (periodic boundary); with transform "on-gain" (default)
   the product replaces the string only if its energy is lower, with
   "always" it replaces it unconditionally;
2. the energy drop dE = E(old) - E(new), with E = 1 - s and s the block
   entropy difference S_3 - S_2 (appendix B, m = 2), lets the organism
   reproduce with probability P(dE) = (1 - exp(-dE/beta)) / (1 - exp(-beta))
   for dE > 0 and 0 otherwise (eq. 1); the offspring replaces a random
   organism, so the population size is constant;
3. with probability mu the offspring is another rule, uniformly;
4. after all organisms have acted, every string is replaced with probability
   gamma by a fresh low-entropy string (bits 1 with probability p0, or bits 0
   with probability p0, half and half).

The run is returned as the observed network of events: metabolism
``R<a> + m<old> -> R<a> + m<new>``, metabolism with a birth
``R<a> + m<old> + R<v> -> R<a> + R<child> + m<new>`` (the victim R<v> is
replaced), and the flow events ``m -> ∅`` / ``∅ -> m``. The energy of every
string is in extras.energies.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from chemart.helpers.params import apportion
from chemart.network import Network, Reaction, Species

N_RULES = 256
TOL = 1e-9   # energy drops below this are numerical noise, not metabolism


# ----------------------------------------------------------------------------
# the chemistry: CA step, entropy, reproduction probability, inflow strings

def rule_table(rule: int) -> str:
    """Outputs for neighbourhoods 111, 110, ..., 000 (fig. A1): rule 30 -> '00011110'."""
    return format(int(rule), "08b")


def ca_step(strings: np.ndarray, rules) -> np.ndarray:
    """One synchronous elementary-CA step of each row with its rule, periodic boundary."""
    x = np.asarray(strings, dtype=np.uint8)
    single = x.ndim == 1
    x = np.atleast_2d(x)
    rules = np.broadcast_to(np.asarray(rules, dtype=np.uint16), (x.shape[0],))
    code = (np.roll(x, 1, axis=1) << 2) | (x << 1) | np.roll(x, -1, axis=1)
    out = ((rules[:, None] >> code) & 1).astype(np.uint8)
    return out[0] if single else out


def _block_entropy(counts: np.ndarray, length: int) -> np.ndarray:
    p = counts / length
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, -p * np.log2(np.where(p > 0, p, 1.0)), 0.0)
    return terms.sum(axis=1)


def entropy(strings: np.ndarray) -> np.ndarray:
    """s = S_3 - S_2 (eq. A1-A2 with m = 2), substrings counted cyclically."""
    x = np.atleast_2d(np.asarray(strings, dtype=np.uint8))
    n, length = x.shape
    code = (np.roll(x, 1, axis=1) << 2) | (x << 1) | np.roll(x, -1, axis=1)
    flat = (code.astype(np.int64) + 8 * np.arange(n)[:, None]).ravel()
    c3 = np.bincount(flat, minlength=8 * n).reshape(n, 8)
    c2 = c3[:, 0::2] + c3[:, 1::2]          # the pair (left, centre) of each triple
    s = np.clip(_block_entropy(c3, length) - _block_entropy(c2, length), 0.0, 1.0)
    return s if np.asarray(strings).ndim > 1 else s[:1]


def energy(strings: np.ndarray) -> np.ndarray:
    """E = 1 - s: ordered strings are energy rich."""
    return 1.0 - entropy(strings)


def reproduction_probability(dE, beta: float) -> np.ndarray:
    """Eq. 1. Values above 1 (small beta) mean certain reproduction."""
    dE = np.asarray(dE, dtype=float)
    p = -np.expm1(-np.maximum(dE, 0.0) / beta) / -np.expm1(-beta)
    return np.where(dE > TOL, p, 0.0)


def fresh_strings(rng, n: int, length: int, p0: float) -> np.ndarray:
    """Inflow strings: each bit is the majority symbol with probability p0; majority 1 or 0 equally often."""
    majority_one = rng.random(n) < 0.5
    bits = rng.random((n, length)) < p0
    return np.where(majority_one[:, None], bits, ~bits).astype(np.uint8)


def shannon_index(population: np.ndarray) -> float:
    """H = sum_i f_i ln(1/f_i) over species (rules)."""
    f = np.bincount(population, minlength=N_RULES) / len(population)
    f = f[f > 0]
    return float(-(f * np.log(f)).sum())


# ----------------------------------------------------------------------------
# the simulation

class _Recorder:
    def __init__(self):
        self.metabolites: dict[bytes, str] = {}
        self.structures: dict[str, str] = {}
        self.energies: dict[str, float] = {}
        self.reactions: dict[tuple, list] = {}
        self.rules: set[int] = set()

    def keys(self, rows: np.ndarray, energies: np.ndarray) -> list[bytes]:
        packed = np.packbits(rows, axis=1)
        out = []
        for i, row in enumerate(packed):
            key = row.tobytes()
            if key not in self.metabolites:
                sid = "m" + key.hex()
                self.metabolites[key] = sid
                self.structures[sid] = "".join("1" if b else "0" for b in rows[i])
                self.energies[sid] = float(energies[i])
            out.append(key)
        return out

    def fire(self, reactants: list[str], products: list[str]) -> None:
        # deduplicate as multisets: e.g. a birth whose victim is of the child's species
        key = (frozenset(Counter(reactants).items()), frozenset(Counter(products).items()))
        entry = self.reactions.get(key)
        if entry is None:
            self.reactions[key] = entry = [reactants, products, 0]
        entry[2] += 1


def simulate(rng, population, n_strings: int, string_len: int, flow_rate: float,
             mutation_rate: float, beta: float, p0: float, updates: int,
             record: bool = True, gain_only: bool = False) -> dict:
    """Run Urdar. `population` is the initial list of rule numbers (one per organism)."""
    pop = np.asarray(population, dtype=np.int64).copy()
    n_org = len(pop)
    pool = fresh_strings(rng, n_strings, string_len, p0)
    pool_energy = energy(pool)
    depth = np.zeros(n_strings, dtype=np.int64)   # successful metabolic steps per string
    rec = _Recorder() if record else None
    if record:
        pool_keys = rec.keys(pool, pool_energy)
        rec.rules.update(int(r) for r in pop)
        initial = Counter(f"R{r}" for r in pop)
        initial.update(rec.metabolites[k] for k in pool_keys)
    history = {k: [] for k in ("shannon", "births", "energy_uptake", "efficiency",
                               "removed_depth", "pool_energy")}
    history["shannon"].append(shannon_index(pop))

    for _ in range(updates):
        order = rng.permutation(n_org)
        picks = rng.integers(n_strings, size=n_org)
        draws = rng.random(n_org)
        mutate = rng.random(n_org) < mutation_rate
        shifts = rng.integers(1, N_RULES, size=n_org)
        victims = rng.integers(n_org, size=n_org)

        # organisms act in `order`; several may pick the same string, in that order
        by_string = np.argsort(picks, kind="stable")
        sorted_picks = picks[by_string]
        first = np.searchsorted(sorted_picks, sorted_picks, side="left")
        rank = np.empty(n_org, dtype=np.int64)
        rank[by_string] = np.arange(n_org) - first

        acts_rule, acts_pos, acts_dE, acts_old, acts_new = [], [], [], [], []
        for r in range(int(rank.max()) + 1 if n_org else 0):
            pos = np.flatnonzero(rank == r)
            slots = picks[pos]
            rules = pop[order[pos]]
            old = pool[slots]
            new = ca_step(old, rules)
            e_new = energy(new)
            dE = pool_energy[slots] - e_new
            gain = dE > TOL
            if gain_only:
                # an organism that cannot extract energy leaves the string unchanged
                new = np.where(gain[:, None], new, old)
                e_new = np.where(gain, e_new, pool_energy[slots])
            if record:
                old_keys = [pool_keys[s] for s in slots]
                new_keys = list(old_keys)
                changed = np.flatnonzero(gain) if gain_only else np.arange(len(slots))
                for i, k in zip(changed, rec.keys(new[changed], e_new[changed])):
                    new_keys[i] = k
                for s, k in zip(slots, new_keys):
                    pool_keys[s] = k
                acts_old += old_keys
                acts_new += new_keys
            pool[slots] = new
            pool_energy[slots] = e_new
            depth[slots] += gain
            acts_rule.append(rules)
            acts_pos.append(pos)
            acts_dE.append(dE)

        rules = np.concatenate(acts_rule) if acts_rule else np.zeros(0, dtype=np.int64)
        pos = np.concatenate(acts_pos) if acts_pos else np.zeros(0, dtype=np.int64)
        dE = np.concatenate(acts_dE) if acts_dE else np.zeros(0)
        chrono = np.argsort(pos, kind="stable")
        born = reproduction_probability(dE, beta) > draws[pos]

        new_pop = pop.copy()
        n_births = 0
        for i in chrono:
            p, a = pos[i], int(rules[i])
            if born[i]:
                child = (a + int(shifts[p])) % N_RULES if mutate[p] else a
                v = int(victims[p])
                victim = int(new_pop[v])
                new_pop[v] = child
                n_births += 1
                if record:
                    rec.rules.update((child, victim))
                    old, new = rec.metabolites[acts_old[i]], rec.metabolites[acts_new[i]]
                    rec.fire([f"R{a}", old, f"R{victim}"], [f"R{a}", f"R{child}", new])
            elif record and acts_old[i] != acts_new[i]:
                old, new = rec.metabolites[acts_old[i]], rec.metabolites[acts_new[i]]
                rec.fire([f"R{a}", old], [f"R{a}", new])
        pop = new_pop

        # flow: each string is replaced with probability gamma
        slots = np.flatnonzero(rng.random(n_strings) < flow_rate)
        new = fresh_strings(rng, len(slots), string_len, p0)
        e_in = energy(new) if len(slots) else np.zeros(0)
        e_out = pool_energy[slots]
        if record and len(slots):
            new_keys = rec.keys(new, e_in)
            for s, k in zip(slots, new_keys):
                old, fresh = rec.metabolites[pool_keys[s]], rec.metabolites[k]
                rec.fire([old], [])
                rec.fire([], [fresh])
                pool_keys[s] = k
        history["removed_depth"].append(float(depth[slots].mean()) if len(slots) else None)
        history["energy_uptake"].append(float((e_in - e_out).sum()))
        history["efficiency"].append(float((e_in - e_out).sum() / e_in.sum()) if len(slots) else None)
        pool[slots], pool_energy[slots], depth[slots] = new, e_in, 0
        history["births"].append(n_births)
        history["shannon"].append(shannon_index(pop))
        history["pool_energy"].append(float(pool_energy.mean()))

    out = {"population": pop, "history": history}
    if record:
        out.update(recorder=rec, initial=initial)
    return out


# ----------------------------------------------------------------------------

def initial_population(p) -> list[int]:
    rules = list(range(N_RULES)) if not p.rules else p.rules
    bad = [r for r in rules if not isinstance(r, int) or isinstance(r, bool) or not 0 <= r < N_RULES]
    if bad:
        raise ValueError(f"rules must be elementary CA rule numbers 0..255, got {bad!r}")
    if len(set(rules)) != len(rules):
        raise ValueError(f"rules must be distinct, got {rules!r}")
    weights = p.abundances or [1.0] * len(rules)
    if len(weights) != len(rules):
        raise ValueError(f"abundances needs one weight per rule ({len(rules)}), got {len(weights)}")
    if any(not isinstance(w, (int, float)) or isinstance(w, bool) or w < 0 for w in weights):
        raise ValueError(f"abundances must be non-negative numbers, got {weights!r}")
    counts = apportion(p.n_organisms, [float(w) for w in weights])
    return [r for r, c in zip(rules, counts) for _ in range(c)]


def generate(p, rng):
    if p.beta <= 0:
        raise ValueError(f"beta must be > 0, got {p.beta}")
    population = initial_population(p)
    run = simulate(rng, population, p.n_strings, p.string_len, p.flow_rate,
                   p.mutation_rate, p.beta, p.p0, p.updates,
                   gain_only=p.transform == "on-gain")
    rec = run["recorder"]
    rules = sorted(rec.rules)
    species = [Species(f"R{r}", structure=rule_table(r)) for r in rules]
    species += [Species(sid, structure=rec.structures[sid]) for sid in rec.metabolites.values()]
    reactions = [Reaction(dict(Counter(lhs)), dict(Counter(rhs)), count=n)
                 for lhs, rhs, n in rec.reactions.values()]
    final = Counter(f"R{r}" for r in run["population"])
    h = run["history"]
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={k: float(v) for k, v in run["initial"].items()},
        extras={
            "energies": rec.energies,
            "urdar": {
                "organisms": [f"R{r}" for r in rules],
                "energy": "E = 1 - (S_3 - S_2), block entropies of cyclic substrings (bits)",
                "reproduction": "P(dE) = (1 - exp(-dE/beta)) / (1 - exp(-beta)) if dE > 0 else 0",
                "flow": {"rate": p.flow_rate, "p0": p.p0, "pool_size": p.n_strings},
                "reactor": "well-stirred: organisms pick strings uniformly at random from one pool",
            },
            "analysis": {
                "updates": p.updates,
                "final_population": {k: final[k] for k in sorted(final, key=lambda s: int(s[1:]))},
                "shannon_index": h["shannon"],
                "births_per_update": h["births"],
                "energy_uptake": h["energy_uptake"],
                "efficiency": h["efficiency"],
                "removed_string_depth": h["removed_depth"],
                "mean_pool_energy": h["pool_energy"],
            },
        },
    )
