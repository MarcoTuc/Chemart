"""The EVOLVE series of virtual ecosystems (book 8.2.3, refs [191, 195, 134, 135, 714]).

Michael Conrad's ecosystem models: the 1969/1970 one-dimensional world of
Conrad & Pattee, then EVOLVE II (Conrad & Strizich 1985), EVOLVE III (Rizki &
Conrad 1985) and EVOLVE IV (Brewster & Conrad 1998).  None of those papers is
openly available, so this is the book's paragraph implemented literally, with
every gap filled explicitly (see the catalog `decisions`).  The book says:

  "virtual organisms are also subject to energy and matter conservation.  They
   live in a 2D world where they respond to their environment, secrete
   chemicals, harvest energy in the form of light, and use this energy for
   various metabolic and genetic processes.  Enzymes and other proteins are
   represented as strings with critical sections that determine their shape.
   Their function is then fetched from a table by a matching mechanism that
   weights the critical section more heavily than the rest of the string...
   Starting with a population of only autotrophs, scavenger populations
   emerge, go extinct, and reappear, showing the emergence of an ecological
   niche aimed at decomposing the metabolic remains of autotrophs."

Matter is a closed inventory of units that cycle

    mineral --(fix: an autotroph spending light energy)--> organic
    organic --(respire: a scavenger releasing that energy)--> mineral

Bodies are organic: an organism holds `body` units while it lives, builds each
offspring out of `body` more, and returns them all to the organic pool as its
remains when it dies.  Energy enters only as light and is tracked to the last
unit (`extras["energies"]["ledger"]`).  An organism reaches the matter in its
own cell and the four cells it touches.

The returned network is the set of events that fired in one run, with counts:

    G + mineral -> G + organic            photosynthesis ("fix")
    G + organic -> G + mineral            scavenging ("respire")
    G + b organic -> G + O                reproduction (b = body)
    G + H + b organic -> G + H + O        conjugative reproduction, donor H
    G -> b organic                        death, recycling the body
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter

from chemart.helpers.params import apportion
from chemart.network import Network, Reaction, Species

#: The function table proteins are matched against.  "light" harvests sunlight,
#: "fix" builds organic matter out of mineral matter (autotrophy), "respire"
#: breaks organic remains back down to mineral matter (scavenging) and
#: "replicate" is the genetic machinery of reproduction.
FUNCTIONS = ("light", "fix", "respire", "replicate")

MINERAL = "mineral"
ORGANIC = "organic"
POOLS = (MINERAL, ORGANIC)
DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def reference(index: int, length: int, alphabet: int) -> tuple[int, ...]:
    """Reference string of table entry `index`: ref[k] = (index + k (index + 1)) mod a."""
    return tuple((index + k * (index + 1)) % alphabet for k in range(length))


def critical_slice(protein_length: int, critical_length: int) -> slice:
    """The critical section: a centred window of the protein string."""
    start = (protein_length - critical_length) // 2
    return slice(start, start + critical_length)


class Chemistry:
    """The string/function machinery: proteins, critical sections, the table."""

    def __init__(self, p) -> None:
        self.alphabet = p.alphabet
        self.protein_length = p.protein_length
        self.critical_length = p.critical_length
        self.weight = p.critical_weight
        self.threshold = p.match_threshold
        self.n_genes = p.n_genes
        self.gene_length = p.protein_length + 1          # sensitivity locus + protein
        self.genome_length = self.gene_length * p.n_genes
        self.critical = critical_slice(p.protein_length, p.critical_length)
        self.table = {
            name: reference(i, p.protein_length, p.alphabet)
            for i, name in enumerate(FUNCTIONS)
        }
        self.max_score = (self.weight * p.critical_length
                          + (p.protein_length - p.critical_length))
        self._phenotype: dict[str, dict] = {}

    # -- the matching mechanism ---------------------------------------------
    def score(self, protein, ref) -> float:
        """Weighted match: a symbol inside the critical section counts `weight` times."""
        total = 0.0
        lo, hi = self.critical.start, self.critical.stop
        for k, (a, b) in enumerate(zip(protein, ref)):
            if a == b:
                total += self.weight if lo <= k < hi else 1.0
        return total / self.max_score

    def function(self, protein) -> tuple[str | None, float]:
        """Fetch a protein's function from the table, with its match quality."""
        best_name, best = None, 0.0
        for name, ref in self.table.items():
            s = self.score(protein, ref)
            if s > best:
                best_name, best = name, s
        if best_name is None or best < self.threshold:
            return None, best
        return best_name, best

    def protein(self, gene: str) -> tuple[int, ...]:
        return tuple(ord(c) - ord("a") for c in gene[1:])

    # -- genomes -------------------------------------------------------------
    def genes(self, genome: str) -> list[str]:
        g = self.gene_length
        return [genome[i * g:(i + 1) * g] for i in range(self.n_genes)]

    def sensitivity(self, gene: str) -> int:
        """EVOLVE II: how much a mutation changes is itself coded in the gene."""
        return ord(gene[0]) - ord("a") + 1

    def phenotype(self, genome: str) -> dict:
        """Best enzyme efficiency per function, and the mean mutational sensitivity."""
        known = self._phenotype.get(genome)
        if known is not None:
            return known
        eff = dict.fromkeys(FUNCTIONS, 0.0)
        sens = []
        for gene in self.genes(genome):
            sens.append(self.sensitivity(gene))
            name, quality = self.function(self.protein(gene))
            if name is not None and quality > eff[name]:
                eff[name] = quality
        out = {
            "efficiency": eff,
            "sensitivity": sum(sens) / len(sens) if sens else 0.0,
            "autotroph": eff["fix"] > 0.0,
            "scavenger": eff["respire"] > 0.0,
        }
        self._phenotype[genome] = out
        return out

    # -- variation -----------------------------------------------------------
    def mutate(self, genome: str, rng, rate: float) -> str:
        symbols = list(genome)
        g = self.gene_length
        for i in range(self.n_genes):
            base = i * g
            if rng.random() < rate:                       # the sensitivity locus itself
                symbols[base] = self._symbol(rng)
            if rng.random() < rate:
                steps = min(ord(symbols[base]) - ord("a") + 1, self.protein_length)
                for k in rng.choice(self.protein_length, size=steps, replace=False):
                    symbols[base + 1 + int(k)] = self._symbol(rng)
        return "".join(symbols)

    def recombine(self, a: str, b: str, rng) -> str:
        """Conjugation (Conrad & Pattee 1970): one crossover at a gene boundary."""
        cut = int(rng.integers(1, self.n_genes)) * self.gene_length
        return a[:cut] + b[cut:]

    def _symbol(self, rng) -> str:
        return chr(ord("a") + int(rng.integers(self.alphabet)))

    def spell(self, name: str) -> str:
        return "".join(chr(ord("a") + s) for s in self.table[name])

    def founder(self, rng, noise: int, sensitivity: int) -> str:
        """A pure autotroph: one gene per function except "respire", then filler.

        Scavenging is the one table entry the founders lack, so the scavenger
        guild of the book's runs can only arise by mutation.
        """
        locus = chr(ord("a") + sensitivity - 1)
        genes = [locus + self.spell(name) for name in FUNCTIONS if name != "respire"]
        while len(genes) < self.n_genes:
            while True:
                filler = "".join(self._symbol(rng) for _ in range(self.protein_length))
                if self.function(tuple(ord(c) - ord("a") for c in filler))[0] is None:
                    break
            genes.append(locus + filler)
        symbols = list("".join(genes))
        for _ in range(noise):                 # mutate protein symbols, not the loci
            i = int(rng.integers(self.n_genes)) * self.gene_length
            symbols[i + 1 + int(rng.integers(self.protein_length))] = self._symbol(rng)
        return "".join(symbols)


def genotype_id(genome: str) -> str:
    return "g" + hashlib.blake2s(genome.encode(), digest_size=4).hexdigest()


class Organism:
    __slots__ = ("cell", "genome", "energy", "age")

    def __init__(self, cell, genome: str, energy: int):
        self.cell, self.genome, self.energy, self.age = cell, genome, energy, 0


class _World:
    def __init__(self, p, rng, chem: Chemistry):
        self.p, self.rng, self.chem = p, rng, chem
        self.width, self.height = p.width, p.height
        self.pool = {
            MINERAL: {(x, y): p.initial_mineral
                      for x in range(p.width) for y in range(p.height)},
            ORGANIC: {(x, y): p.initial_organic
                      for x in range(p.width) for y in range(p.height)},
        }
        self.organic_energy = 0            # energy banked in the organic pool
        self.by_cell: dict[tuple[int, int], Organism] = {}
        self.fired: dict[tuple, list] = {}
        self.ledger = Counter()
        self.events = Counter()
        self.genomes: dict[str, str] = {}   # genotype id -> genome
        self.history: dict[str, list] = {
            "step": [], "light": [], "organisms": [], "autotrophs": [], "scavengers": [],
            "mineral": [], "organic": [], "genotypes": [],
            "mean_light_efficiency": [], "mean_fix_efficiency": [],
            "mean_respire_efficiency": [], "mean_sensitivity": [],
        }

    # -- bookkeeping ---------------------------------------------------------
    def name(self, genome: str) -> str:
        gid = genotype_id(genome)
        self.genomes.setdefault(gid, genome)
        return gid

    def event(self, reactants: dict[str, int], products: dict[str, int], kind: str) -> None:
        key = (frozenset(reactants.items()), frozenset(products.items()))
        entry = self.fired.get(key)
        if entry is None:
            self.fired[key] = entry = [dict(reactants), dict(products), 0, kind]
        entry[2] += 1
        self.events[kind] += 1

    def neighbours(self, cell):
        x, y = cell
        return [((x + dx) % self.width, (y + dy) % self.height) for dx, dy in DIRECTIONS]

    def reach(self, cell):
        """The matter an organism can touch: its own cell, then the four it borders."""
        return [cell, *self.neighbours(cell)]

    def available(self, cell, pool: str) -> int:
        return sum(self.pool[pool][c] for c in self.reach(cell))

    def take(self, cell, pool: str, amount: int) -> bool:
        if self.available(cell, pool) < amount:
            return False
        for c in self.reach(cell):
            if amount == 0:
                break
            taken = min(amount, self.pool[pool][c])
            self.pool[pool][c] -= taken
            amount -= taken
        return True

    def light_at(self, t: int) -> int:
        if self.p.light_period <= 0:
            return self.p.light_input
        swing = self.p.light_swing * math.sin(2.0 * math.pi * t / self.p.light_period)
        return max(0, int(round(self.p.light_input * (1.0 + swing))))

    def total(self, pool: str) -> int:
        return sum(self.pool[pool].values())

    # -- one organism's step --------------------------------------------------
    def live(self, org: Organism, light: int) -> None:
        p, chem = self.p, self.chem
        pheno = chem.phenotype(org.genome)
        eff = pheno["efficiency"]
        gid = self.name(org.genome)
        cell = org.cell

        # 1. harvest light ("harvest energy in the form of light")
        gain = int(eff["light"] * light)
        org.energy += gain
        self.ledger["light_incident"] += light
        self.ledger["light_harvested"] += gain
        self.ledger["light_lost"] += light - gain

        # 2. scavenge: respire one unit of organic remains back to mineral matter
        if eff["respire"] > 0.0 and self.available(cell, ORGANIC) > 0:
            units = self.total(ORGANIC)
            share = self.organic_energy // units if units else 0
            self.take(cell, ORGANIC, 1)
            self.pool[MINERAL][cell] += 1
            self.organic_energy -= share
            got = int(eff["respire"] * share)
            org.energy += got
            self.ledger["respired"] += got
            self.ledger["dissipated"] += share - got
            self.event({gid: 1, ORGANIC: 1}, {gid: 1, MINERAL: 1}, "respire")

        # 3. fix: spend light energy to turn mineral matter into the organic
        #    matter bodies are made of, and secrete it into the organism's cell.
        #    An organism builds up one offspring's worth and then stops.
        if (eff["fix"] > 0.0 and org.energy >= p.fix_cost + p.maintenance
                and self.available(cell, MINERAL) > 0
                and self.available(cell, ORGANIC) < p.body):
            org.energy -= p.fix_cost
            self.take(cell, MINERAL, 1)
            self.pool[ORGANIC][cell] += 1
            self.organic_energy += p.fix_cost
            self.ledger["fix_spent"] += p.fix_cost
            self.event({gid: 1, MINERAL: 1}, {gid: 1, ORGANIC: 1}, "fix")

        # 4. maintenance; an organism that cannot pay it starves
        paid = min(org.energy, p.maintenance)
        org.energy -= paid
        self.ledger["maintenance_spent"] += paid
        if paid < p.maintenance:
            self.die(org, gid)
            return
        org.age += 1
        if p.max_age > 0 and org.age > p.max_age:
            self.die(org, gid)
            return

        # 5. reproduction
        self.reproduce(org, gid, pheno)

    def reproduce(self, org: Organism, gid: str, pheno: dict) -> None:
        p, chem, rng = self.p, self.chem, self.rng
        if org.energy < p.repro_cost + p.offspring_energy:
            return
        free = [c for c in self.neighbours(org.cell) if c not in self.by_cell]
        if not free or self.available(org.cell, ORGANIC) < p.body:
            return
        if rng.random() >= pheno["efficiency"]["replicate"]:
            self.events["failed_replication"] += 1
            return

        donor = None
        if rng.random() < p.recombination_rate:
            mates = [self.by_cell[c] for c in self.neighbours(org.cell)
                     if c in self.by_cell]
            if mates:
                donor = mates[int(rng.integers(len(mates)))]

        source = chem.recombine(org.genome, donor.genome, rng) if donor else org.genome
        child_genome = chem.mutate(source, rng, p.mutation_rate)
        self.take(org.cell, ORGANIC, p.body)
        org.energy -= p.repro_cost + p.offspring_energy
        self.ledger["repro_dissipated"] += p.repro_cost
        cell = free[int(rng.integers(len(free)))]
        self.by_cell[cell] = Organism(cell, child_genome, p.offspring_energy)
        cid = self.name(child_genome)

        # built up key by key: a clonal offspring has cid == gid, and a dict
        # literal would collapse the two into one instead of counting both
        reactants = Counter({gid: 1, ORGANIC: p.body})
        products = Counter({gid: 1})
        products[cid] += 1
        if donor is not None:
            did = self.name(donor.genome)
            reactants[did] += 1
            products[did] += 1
        self.event(dict(reactants), dict(products),
                   "conjugation" if donor is not None else "reproduction")

    def die(self, org: Organism, gid: str) -> None:
        """Mass-conserving recycling: the body becomes organic remains."""
        self.pool[ORGANIC][org.cell] += self.p.body
        self.organic_energy += org.energy
        self.ledger["death_deposited"] += org.energy
        org.energy = 0
        del self.by_cell[org.cell]
        self.event({gid: 1}, {ORGANIC: self.p.body}, "death")

    # -- the world's step ------------------------------------------------------
    def diffuse(self) -> None:
        if self.p.diffusion <= 0.0:
            return
        for pool in POOLS:
            for cell, amount in list(self.pool[pool].items()):
                moving = int(self.p.diffusion * amount)
                if moving <= 0:
                    continue
                self.pool[pool][cell] -= moving
                for target, share in zip(self.neighbours(cell),
                                         apportion(moving, [1.0] * 4)):
                    self.pool[pool][target] += share

    def sample(self, t: int, light: int) -> None:
        chem = self.chem
        pop = list(self.by_cell.values())
        phenos = [chem.phenotype(o.genome) for o in pop]
        n = len(pop) or 1
        h = self.history
        h["step"].append(t)
        h["light"].append(light)
        h["organisms"].append(len(pop))
        h["autotrophs"].append(sum(1 for f in phenos if f["autotroph"]))
        h["scavengers"].append(sum(1 for f in phenos if f["scavenger"]))
        h["mineral"].append(self.total(MINERAL))
        h["organic"].append(self.total(ORGANIC))
        h["genotypes"].append(len({o.genome for o in pop}))
        for name in ("light", "fix", "respire"):
            h[f"mean_{name}_efficiency"].append(
                round(sum(f["efficiency"][name] for f in phenos) / n, 6))
        h["mean_sensitivity"].append(round(sum(f["sensitivity"] for f in phenos) / n, 6))

    def run(self) -> None:
        self.ledger["founders_endowment"] = sum(o.energy for o in self.by_cell.values())
        self.initial = Counter(self.name(o.genome) for o in self.by_cell.values())
        self.initial[MINERAL] = self.total(MINERAL)
        self.initial[ORGANIC] = self.total(ORGANIC)
        self.sample(0, self.light_at(0))
        for t in range(1, self.p.steps + 1):
            light = self.light_at(t)
            order = list(self.by_cell.values())
            for k in self.rng.permutation(len(order)):
                org = order[int(k)]
                if self.by_cell.get(org.cell) is org:
                    self.live(org, light)
            self.diffuse()
            self.sample(t, light)
        self.final = Counter(self.name(o.genome) for o in self.by_cell.values())


# ---------------------------------------------------------------------------
def _check(p) -> None:
    if p.critical_length > p.protein_length:
        raise ValueError(
            f"critical_length ({p.critical_length}) cannot exceed protein_length "
            f"({p.protein_length}): the critical section is part of the protein string"
        )
    if p.critical_weight < 1.0:
        raise ValueError(
            "critical_weight must be >= 1: the book weights the critical section "
            f"more heavily than the rest of the string, got {p.critical_weight}"
        )
    if p.n_genes < len(FUNCTIONS):
        raise ValueError(
            f"n_genes ({p.n_genes}) must be at least {len(FUNCTIONS)}, one gene per "
            f"entry of the function table {list(FUNCTIONS)}"
        )
    if p.initial_sensitivity > p.alphabet:
        raise ValueError(
            f"initial_sensitivity ({p.initial_sensitivity}) cannot exceed alphabet "
            f"({p.alphabet}): the sensitivity locus is a single symbol"
        )
    if p.n_founders > p.width * p.height:
        raise ValueError(
            f"n_founders ({p.n_founders}) exceeds the {p.width}x{p.height} grid "
            f"({p.width * p.height} cells); one organism lives per cell"
        )
    if not isinstance(p.founder_genomes, list) or not all(
            isinstance(g, str) for g in p.founder_genomes):
        raise ValueError(
            f"founder_genomes must be a list of genome strings, got {p.founder_genomes!r}")


def generate(p, rng):
    _check(p)
    chem = Chemistry(p)
    for genome in p.founder_genomes:
        if len(genome) != chem.genome_length:
            raise ValueError(
                f"founder genome {genome[:12]!r} has {len(genome)} symbols; "
                f"n_genes x (protein_length + 1) = {chem.genome_length} are required"
            )
        bad = sorted({c for c in genome if not 0 <= ord(c) - ord("a") < p.alphabet})
        if bad:
            raise ValueError(
                f"founder genome uses symbols {bad} outside the {p.alphabet}-letter "
                f"alphabet a..{chr(ord('a') + p.alphabet - 1)}"
            )

    world = _World(p, rng, chem)
    cells = [(x, y) for x in range(p.width) for y in range(p.height)]
    chosen = rng.choice(len(cells), size=p.n_founders, replace=False)
    for k, index in enumerate(chosen):
        genome = (p.founder_genomes[k % len(p.founder_genomes)] if p.founder_genomes
                  else chem.founder(rng, p.founder_noise, p.initial_sensitivity))
        cell = cells[int(index)]
        world.by_cell[cell] = Organism(cell, genome, p.offspring_energy)
    world.run()
    return _network(world, p, chem)


def _episodes(counts: list[int]) -> list[dict]:
    """The contiguous stretches during which a guild is present."""
    out, start = [], None
    for t, n in enumerate(counts):
        if n > 0 and start is None:
            start = t
        elif n == 0 and start is not None:
            out.append({"from": start, "to": t - 1, "peak": max(counts[start:t])})
            start = None
    if start is not None:
        out.append({"from": start, "to": len(counts) - 1, "peak": max(counts[start:])})
    return out


def _network(world: _World, p, chem: Chemistry) -> Network:
    names = sorted(world.genomes)
    species = [Species(n, structure=world.genomes[n]) for n in names]
    species += [
        Species(MINERAL, structure="one unit of free inorganic matter"),
        Species(ORGANIC, structure="one unit of organic matter: secretions and remains"),
    ]
    reactions = [Reaction(lhs, rhs, None, count) for lhs, rhs, count, _ in world.fired.values()]
    events = [{"kind": kind, "count": count} for _, _, count, kind in world.fired.values()]

    matter = {n: p.body for n in names}
    matter[MINERAL] = matter[ORGANIC] = 1

    ledger = {k: int(v) for k, v in sorted(world.ledger.items())}
    ledger["organisms_energy"] = int(sum(o.energy for o in world.by_cell.values()))
    ledger["organic_energy"] = int(world.organic_energy)
    ledger["organisms_balanced"] = (
        ledger.get("founders_endowment", 0) + ledger.get("light_harvested", 0)
        + ledger.get("respired", 0)
        == ledger.get("fix_spent", 0) + ledger.get("maintenance_spent", 0)
        + ledger.get("repro_dissipated", 0) + ledger.get("death_deposited", 0)
        + ledger["organisms_energy"])
    ledger["organic_balanced"] = (
        ledger.get("fix_spent", 0) + ledger.get("death_deposited", 0)
        == ledger.get("respired", 0) + ledger.get("dissipated", 0)
        + ledger["organic_energy"])

    phenotypes = {}
    for gid in names:
        f = chem.phenotype(world.genomes[gid])
        phenotypes[gid] = {
            "efficiency": {k: round(v, 6) for k, v in f["efficiency"].items()},
            "sensitivity": round(f["sensitivity"], 6),
            "guild": ("mixotroph" if f["autotroph"] and f["scavenger"]
                      else "autotroph" if f["autotroph"]
                      else "scavenger" if f["scavenger"] else "none"),
        }
    final_guilds = Counter(phenotypes[g]["guild"] for g in world.final.elements())

    h = world.history
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={k: float(v) for k, v in sorted(world.initial.items())},
        extras={
            "conservation": [{
                "name": "matter",
                "vector": matter,
                "notes": f"a living organism holds {p.body} units of organic matter; they "
                         "return to the organic pool as its remains when it dies",
            }],
            "space": {
                "dimensions": 2, "lattice": "torus", "width": p.width, "height": p.height,
                "neighbourhood": "von-neumann",
                "cell": "one organism plus a local mineral and organic pool",
                "reach": "an organism draws matter from its own cell and the four it "
                         "borders, and secretes into its own cell",
                "diffusion": f"each step a fraction {p.diffusion} of each pool is shared "
                             "equally among the four neighbouring cells",
            },
            "energies": {
                "units": "integer energy units",
                "light": {
                    "input": p.light_input, "period": p.light_period, "swing": p.light_swing,
                    "law": "light_input per cell per step when light_period is 0, else "
                           "light_input (1 + swing sin(2 pi t / light_period))",
                    "harvested": "floor(light-enzyme efficiency x incident light); the rest "
                                 "is lost",
                },
                "costs": {"fix": p.fix_cost, "maintenance": p.maintenance,
                          "reproduction": p.repro_cost, "offspring": p.offspring_energy},
                "storage": "fixing one organic unit banks fix_cost energy in the organic "
                           "pool and death banks the dead organism's remaining energy; "
                           "respiring one unit withdraws the pool's mean share and returns "
                           "the respire enzyme's fraction of it, dissipating the rest",
                "ledger": ledger,
            },
            "chemistry": {
                "functions": list(FUNCTIONS),
                "table": {name: chem.spell(name) for name in FUNCTIONS},
                "matching": "function = the table entry of highest weighted match score "
                            "(critical_weight per matching symbol inside the critical "
                            "section, 1 outside, normalised to [0, 1]); a best score below "
                            "match_threshold leaves the protein without a function, and the "
                            "score is the enzyme's efficiency",
                "critical_section": [chem.critical.start, chem.critical.stop],
                "genome": f"{p.n_genes} genes of {chem.gene_length} symbols over the "
                          f"alphabet a..{chr(ord('a') + p.alphabet - 1)}; the first symbol "
                          "of a gene is its mutational sensitivity (EVOLVE II: the magnitude "
                          "of phenotypic change per mutation is itself a property of the "
                          "gene), the rest is the protein string",
                "trophic": {
                    "fix": "G + mineral -> G + organic (autotrophy, paid for with light)",
                    "respire": "G + organic -> G + mineral (scavenging, releases energy)",
                    "reproduce": f"G + {p.body} organic -> G + O (the offspring's body)",
                    "death": f"G -> {p.body} organic (the remains scavengers decompose)",
                },
            },
            "phenotypes": phenotypes,
            "events": events,
            "final_state": {k: int(v) for k, v in sorted(world.final.items())},
            "analysis": {
                "steps": p.steps,
                "event_counts": {k: int(v) for k, v in sorted(world.events.items())},
                "genotypes_seen": len(world.genomes),
                "survivors": int(sum(world.final.values())),
                "final_guilds": {k: int(v) for k, v in sorted(final_guilds.items())},
                "history": h,
                "scavenger_episodes": _episodes(h["scavengers"]),
                "founder_light_efficiency": h["mean_light_efficiency"][0],
                "founder_fix_efficiency": h["mean_fix_efficiency"][0],
                "founder_sensitivity": h["mean_sensitivity"][0],
                "final_light_efficiency": h["mean_light_efficiency"][-1],
                "final_fix_efficiency": h["mean_fix_efficiency"][-1],
                "final_sensitivity": h["mean_sensitivity"][-1],
                "matter_total": int(sum(world.total(pool) for pool in POOLS)
                                    + p.body * len(world.by_cell)),
            },
        },
    )
