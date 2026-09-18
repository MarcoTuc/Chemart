# String Metabolic Network (SMN)

`smn` · *Ono, Fujiwara & Yuta, 2005*

*Also known as:* *artificial metabolic system*

A metabolism made of strings. Metabolites are strings over a small alphabet, joined and cut, or recombined by swapping tails between two of them - and every reaction has its own enzyme, so an organism's genome is simply the set of enzymes it owns. Letter counts are conserved by all three operations, which means the total amount of substance is fixed and makes the trade-off explicit: hold the amount constant and more letters means longer, fewer compounds.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.3.1 |
| **refs** | [638], https://doi.org/10.1007/11553090_72 |
| **provides** | `topology`, `stoichiometry`, `mass-conservation`, `flow`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (implicit): metabolites = strings over an alphabet of up to 8 letters {a..h}, arbitrary length (bounded here by max_length)

**R — reactions** (explicit, arity 2): Ligation/cleavage A + B <-> AB (eq. 18.17) and recombination AB + CD <-> AD + CB (eq. 18.18), e.g. adbg + ef <-> adbgef and facbaha + eefg <-> facbfg + eeaha. Each reaction is catalysed by its own enzyme, and the genome of an organism is its list of enzymes; the metabolism is the network of those reactions (both directions) over the organism's compound set.

**A — reactor**: ode
 · *dilution:* regulated inflow/outflow keeps the total amount of substance in the cell constant (constant-total)

## What you get

```python
net = chemart.generate_network("smn", seed=1)
```

```
smn: 21 species, 24 reactions, status=complete
provides: flow, initial-state, mass-conservation, stoichiometry, topology
seed: 1
extras: conservation, fitness, flow_law, genome, initial_metabolites, kinetics
```

First reactions:

```
cg + bcg -> cgbcg
cgbcg -> cg + bcg
ea + eg -> eaeg
eaeg -> ea + eg
2 cg -> cgcg
cgcg -> 2 cg
cg + d -> cgd
cgd -> cg + d
… and 16 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `alphabet_size` | `int` | `8` | structural | number of letters (a, b, ...); the book's SMN uses 8 <br>`2` … `8` |
| `n_metabolites` | `int` | `6` | population | number of distinct random initial metabolites (ignored when initial_metabolites or genome is given) <br>`1` … `1000` |
| `max_metabolite_length` | `int` | `3` | population | random initial metabolites have a length drawn uniformly from 1 to this value <br>`1` … `50` |
| `initial_metabolites` | `list` | `` | population | explicit initial metabolites (resources present at the start of the organism's life); empty: the educts of an explicit genome, or random strings |
| `n_enzymes` | `int` | `12` | structural | size of the random initial genome; each enzyme acts on compounds already in the compound set (ignored when genome is given) <br>`0` … `10000` |
| `ligation_fraction` | `float` | `0.5` | stochastic | probability that a random enzyme is a ligation rather than a recombination <br>`0` … `1` |
| `genome` | `list` | `` | structural | explicit enzymes: 'A + B' for the ligation A + B <-> AB, 'A\|B + C\|D' for the recombination AB + CD <-> AD + CB (e.g. 'facb\|aha + ee\|fg'); replaces the random genome |
| `mutations` | `int` | `` | stochastic | applications of the book's duplicate-and-mutate operator to the genome, without selection <br>`0` … `100000` · *range:* the book's control: networks grown by random mutation without selection are not modular |
| `max_length` | `int` | `12` | structural | longest compound allowed; random enzymes and mutations whose products are longer are redrawn (the book allows arbitrary length) <br>`2` … `1000` |
| `initial_amount` | `float` | `1.0` | population | initial concentration of each initial metabolite (the book gives no amounts) <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- letter counts are conserved by ligation, cleavage and recombination
- with the amount of substance held constant, a larger letter mass means longer compounds, the fitness criterion
- evolved networks show a power-law community-size distribution close to the average of 100 organisms with large metabolic networks (book figure 18.6, 30 SMN runs)
- networks grown by random mutation without selection pressure are not modular

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The original paper (Ono, Fujiwara & Yuta, Artificial metabolic system: an evolutionary model for community organization in metabolic networks, ECAL 2005, LNAI 3630, pp. 716-724, https://doi.org/10.1007/11553090_72) is closed access: no open copy was found via Springer, Unpaywall, Semantic Scholar or the Internet Archive. Everything here follows book section 18.3.1.
- The v1 origin 'Ono, Fujimoto & Ikegami' is wrong; bibliography entry [638] gives N. Ono, Y. Fujiwara and K. Yuta.
- The book equates an enzyme with the reaction it catalyses, so enzymes are not species. Each enzyme gives two directed reactions (forward and reverse); enzymes catalysing the same reaction share them. extras.genome lists every enzyme with its notation, its reversible reaction, the indices of its two reactions and, for mutants, its parent enzyme.
- No kinetics: the book says only that synthesising longer molecules is harder and gives no coefficients, so every rate is None and the v1 parameter rate_bias is dropped. The fact is recorded in extras.kinetics.
- The regulated inflow/outflow keeps the amount of substance constant (the book: a larger letter mass therefore means longer compounds), which is outflow = constant-total. The composition of the inflow is not given, so there is no explicit inflow.
- The book does not say how the initial genome or the initial metabolites are made. Here the initial metabolites are n_metabolites distinct random strings (or given explicitly), and the genome grows one random enzyme at a time, each acting on compounds already in the compound set (initial metabolites plus earlier products), a ligation with probability ligation_fraction and otherwise a recombination with uniform cut points. Draws that are elastic or exceed max_length are redrawn. Initial amounts are all initial_amount, a Chemart choice.
- Mutation follows the book's operator: a random enzyme is duplicated, one of its two educts (chosen at random) is replaced by a different compound from the current compound set (of length >= 2 for a recombination), a new recombination point is drawn uniformly in the replaced educt only (the book says 'a new point of recombination', singular), and the mutant is added to the genome. Mutants that are elastic or exceed max_length are redrawn.
- Evolution (fitness = letter mass after a fixed number of timesteps, selection across generations) needs the paper's kinetics and acts across networks, so it is not part of one network. The v1 parameters mutation (callable text) and fitness (callable) are dropped; mutations applies the operator without selection, and the fitness definition is recorded in extras.fitness.
- Species are strings without internal structure beyond their letters; letter counts are conserved by both reaction classes (extras.conservation, one vector per letter).

## Notes

One network is the metabolism of one organism. Community-size distributions (figure 18.6) are statistical results of the evolutionary algorithm with the paper's kinetics and are not reproduced here. See `bnc-cell` for the related bond-number chemistry of Hintze & Adami.

---

*Specification: `catalog/chemistries/smn.yaml` · generator: `chemart/chemistries/smn.py` · tests: `tests/chemistries/test_smn.py`*
