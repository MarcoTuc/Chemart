# Nuclear reaction networks

`nuclear-reaction-networks`

Real nuclear physics in the same record as the artificial chemistries, which makes it an unusually good validation target. Nuclides and elementary particles are the species, fusion and decay the reactions, and measured cross-section libraries exist to check against. The CNO cycle is also a textbook catalytic organisation: carbon, nitrogen and oxygen are consumed and regenerated while hydrogen is converted to helium, so the catalysts appear on both sides and would vanish from a net stoichiometric matrix.

| | |
|---|---|
| **family** | non-chemical |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 20.2 |
| **refs** | [187] |
| **provides** | `topology`, `stoichiometry`, `mass-conservation` |

## Molecules, reactions, reactor

**S — molecules** (explicit): nuclides and elementary particles (n, 1H, 2H, ..., e+, e-, nu_e, gamma)

**R — reactions** (explicit, arity [1, 2]): e.g. 1H + 1H -> 2H + e+ + nu_e. The CNO cycle is a genuine CATALYTIC cycle: 12C + 1H -> 13N + gamma ; 13N -> 13C + e+ + nu_e ; 13C + 1H -> 14N + gamma ; 14N + 1H -> 15O + gamma ; 15O -> 15N + e+ + nu_e ; 15N + 1H -> 12C + 4He.

**A — reactor**: ode
 · *dilution:* open: 1H flows in, 4He, e+, nu_e flow out

## What you get

```python
net = chemart.generate_network("nuclear-reaction-networks", seed=1)
```

```
nuclear-reaction-networks: 11 species, 6 reactions, status=complete
provides: mass-conservation, stoichiometry, topology
seed: 1
extras: boundary, conservation
```

First reactions:

```
12C + 1H -> 13N + gamma
13N -> 13C + e+ + nu_e
13C + 1H -> 14N + gamma
14N + 1H -> 15O + gamma
15O -> 15N + e+ + nu_e
15N + 1H -> 12C + 4He
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `network` | `enum` | `cno-cycle` | structural | which reaction network to generate <br>one of `cno-cycle`, `pp-chain`, `big-bang-nucleosynthesis` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- the CNO cycle is a closed, self-maintaining organisation with C, N, O as catalysts
- observed light-element abundances from big-bang nucleosynthesis

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Topology only. Rates depend on temperature and density through cross sections (e.g. REACLIB), which the book does not give, so cross_sections, temperature and density are not parameters.
- cno-cycle is book eq. 20.4, with the book's open boundary (1H in; 4He, e+, nu_e out) recorded in extras.boundary rather than as numeric flows. pp-chain is eq. 20.3 completed with the standard pp-I branch. big-bang-nucleosynthesis is the standard 12-reaction network among n, 1H, 2H, 3H, 3He, 4He, 7Li, 7Be, standing in for figure 20.4, which only draws the network.
- extras.conservation lists baryon number, charge (fully ionised nuclei) and electron lepton number, all exact.

## Notes

An excellent external validation target: real cross-section libraries exist (REACLIB), and the CNO cycle is a textbook chemical organisation.

---

*Specification: `catalog/chemistries/nuclear-reaction-networks.yaml` · generator: `chemart/chemistries/nuclear_reaction_networks.py` · tests: `tests/chemistries/test_nuclear_reaction_networks.py`*
