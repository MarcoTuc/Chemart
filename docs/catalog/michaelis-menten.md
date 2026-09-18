# Michaelis-Menten enzyme kinetics

`michaelis-menten`

The canonical enzyme mechanism, and a lesson about output formats. Written elementarily it is three mass-action reactions among four species: enzyme and substrate bind reversibly, then the complex releases product. Eliminate the complex and you get the familiar saturating rate law over just two species - fewer species and fewer reactions, but no longer mass action. The same chemistry thus has two faithful representations, and any schema that assumes mass action can only hold one of them.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.2.1; also 4.4 (as an abridgement method) |
| **refs** | [42], [81], [503] |
| **provides** | `topology`, `stoichiometry`, `rate-constants`, `rate-law`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): E, S, ES, P (elementary) or S, P (abridged)

**R — reactions** (explicit, arity [1, 2]): E + S <-ka/ka'-> ES --kb--> E + P (irreversible product formation)

**A — reactor**: ode, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("michaelis-menten", seed=1)
```

```
michaelis-menten: 4 species, 3 reactions, status=complete
provides: initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
E + S -> ES  [mass-action k=1.0]
ES -> E + S  [mass-action k=1.0]
ES -> E + P  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `form` | `enum` | `elementary` | structural | elementary mass-action steps, or the single abridged reaction S -> P with a Michaelis-Menten rate law <br>one of `elementary`, `abridged` |
| `ka` | `float` | `1.0` | kinetic | binding E + S -> ES <br>≥ `0` |
| `ka_rev` | `float` | `1.0` | kinetic | unbinding ES -> E + S <br>≥ `0` |
| `kb` | `float` | `1.0` | kinetic | catalytic turnover ES -> E + P <br>≥ `0` |
| `E0` | `float` | `1.0` | population | total enzyme, free + bound <br>≥ `0` |
| `S0` | `float` | `10.0` | population | initial substrate <br>≥ `0` |

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives no numerical rates; all default to 1, with E0 = 1 and S0 = 10.
- k_m is derived, not a parameter. In the extracted text of eq. 18.4 the prime on ka' is lost ('k_m = (k_a + k_b)/k_a'); the quasi-steady-state derivation gives k_m = (ka_rev + kb) / ka, which is used.

## Notes

The abridged form has FEWER species and reactions than the elementary form but a NON-MASS-ACTION rate law; a schema that assumes mass action silently changes the model.

---

*Specification: `catalog/chemistries/michaelis-menten.yaml` · generator: `chemart/chemistries/michaelis_menten.py` · tests: `tests/chemistries/test_michaelis_menten.py`*
