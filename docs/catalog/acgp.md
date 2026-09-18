# Algorithmic Chemistry GP (ACGP)

`acgp` · *Banzhaf & Lasarczyk, 2004-2008*

*Also known as:* *sequenceless register machine programs*

Run this one backwards from the others: instead of building a chemistry, it reads one out of a program. The registers of a virtual machine are the species and their contents are concentrations, so an instruction like r3 <- r1 op r2 simply *is* the reaction r1 + r2 -> r3. Genetic programming then evolves the instruction multiset, and the evolved classifier can be inspected as a reaction network. Because an instruction is a reaction rather than a step, order matters less and deleting instructions degrades performance gracefully.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 16.6 |
| **refs** | [74], [75], [497], [498] |
| **provides** | `topology`, `stoichiometry`, `catalysts` |

## Molecules, reactions, reactor

**S — molecules** (explicit): registers r0..r(n-1) of a virtual machine (species) and their contents (concentrations)

**R — reactions** (explicit, arity [1, 2]): an instruction r3 <- r1 op r2 is read as a reaction r1 + r2 -> r3

**A — reactor**: well-stirred-multiset
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("acgp", seed=1)
```

```
acgp: 12 species, 100 reactions, status=complete
provides: catalysts, stoichiometry, topology
seed: 1
extras: input_registers, output_registers, program
```

First reactions:

```
r5 + r6 -> r10
r0 + r1 -> r10
r2 + r3 -> r10
r3 + r9 -> r6
r7 + r6 -> r4
r10 + r9 -> r10
r9 + r3 -> r7
r1 + r3 -> r4
… and 92 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n_registers` | `int` | `12` | structural | number of registers <br>`2` … `128` |
| `n_instructions` | `int` | `100` | structural | multiset size; the multiplicity of an instruction is its importance, which gives fault tolerance <br>`1` … `10000` |
| `input_registers` | `list` | `[0, 1, 2, 3]` | structural | read-only registers where input is fed; never written |
| `output_registers` | `list` | `[11]` | structural | registers the answer is read from |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- evolves classifiers (e.g. UCI thyroid) as bipartite register/instruction graphs
- robust to instruction deletion
- trivially parallelisable - no coordination between instruction executions

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The generator samples one random program: n_instructions instructions with uniformly random operand registers, a random writable destination (not an input register) and an operator from {+, -, *, /}. Each instruction is the reaction of eq. 16.9, repeated instructions stay repeated (multiset), and extras.program lists the instructions in reaction order.
- Evolution of programs (population, generations, executions_per_eval) acts on populations of such networks and is not part of one network, so those parameters are dropped.

## Notes

The inverse case: here a CRN is READ OUT of a program rather than generated. The bipartite graph (registers = species, instructions = reactions) is exactly the SBML/BioNetGen bipartite view.

---

*Specification: `catalog/chemistries/acgp.yaml` · generator: `chemart/chemistries/acgp.py` · tests: `tests/chemistries/test_acgp.py`*
