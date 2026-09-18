# RAF sets (reflexively autocatalytic, F-generated)

`raf` · *Hordijk & Steel, 2004-2015*

A test rather than a chemistry: given a catalytic reaction system and a food set, is there a self-sustaining subset? A set is reflexively autocatalytic if every reaction in it has a catalyst the set itself can make, and food-generated if every reactant can be built from the food; a RAF is both. The striking result is how cheap this is - the maximal RAF can be found in polynomial time, and the catalysis needed for one to appear grows only linearly with system size rather than exponentially.

| | |
|---|---|
| **family** | origin-of-life |
| **kind** | analysis |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 6.3.1 |
| **refs** | [400], [402], [599] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `mass-conservation` |

## Molecules, reactions, reactor

**S — molecules** (explicit): any catalytic reaction system Q = (X, R, C) plus a designated food set F; by default X = all strings up to length n over an alphabet of size B (Kauffman's binary polymer model) and F = all strings up to length t

**R — reactions** (explicit, arity variable): a set R' is RA if every reaction in R' has a catalyst in cl_R'(F); F-generated if every reactant of R' is in cl_R'(F); RAF = both. Default reactions: the reversible cleavage/ligation pairs A + B <-> AB of the binary polymer model, each catalysed by molecule E (A + B + E -> AB + E) with probability p = f / |R|

**A — reactor**: n/a - static graph analysis
 · *dilution:* n/a

## What you get

```python
net = chemart.generate_network("raf", seed=1)
```

```
raf: 254 species, 746 reactions, status=complete
provides: catalysts, mass-conservation, stoichiometry, topology
seed: 1
extras: analysis, buffered, conservation, food
```

First reactions:

```
2 0 + 001 -> 00 + 001
00 + 001 -> 2 0 + 001
2 0 + 010 -> 00 + 010
00 + 010 -> 2 0 + 010
01 + 0 + 1111101 -> 010 + 1111101
010 + 1111101 -> 01 + 0 + 1111101
01 + 1 + 00111 -> 011 + 00111
011 + 00111 -> 01 + 1 + 00111
… and 738 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n` | `int` | `7` | structural | maximum polymer length; the system has B^1 + ... + B^n molecules and (n-2) 2^(n+1) + 4 cleavage/ligation reactions for B = 2 <br>`1` … `14` · *range:* the papers go to n = 20 on a cluster; n = 8, 10, 12 carry published results |
| `t` | `int` | `2` | population | the food set F is every polymer of length at most t <br>`1` … `6` · *range:* t = 2 throughout the RAF papers, giving \|F\| = 6 for B = 2 |
| `B` | `int` | `2` | structural | alphabet size of the polymers <br>`2` … `4` · *range:* all published RAF results use the binary model B = 2 |
| `f` | `float` | `1.5` | structural | level of catalysis: the average number of reactions catalysed by one molecule, f = p \|R\|, from which the per-pair catalysis probability p is derived <br>`0` … `1000` · *range:* RAFs appear around f = 1-2 for n up to 20; f(n) = 1.0970 + 0.0189 n gives probability 0.5 |
| `irreducible_rafs` | `int` | `1` | stochastic | how many irreducible RAFs to sample inside the maximal RAF (each run of the deletion search can return a different one); 0 skips the search <br>`0` … `100` |
| `reactions` | `enum` | `catalysed` | structural | which reactions of the CRS the returned network contains: only the catalysed ones (the only ones an RAF can contain) or every cleavage/ligation pair as well <br>one of `catalysed`, `all` |
| `system` | `dict` | `` | structural | analyse a user-supplied catalytic reaction system instead of the polymer model: {'food': [ids], 'reactions': [{'id', 'reactants', 'products', 'catalysts', 'reversible'}]}; empty means use the polymer model |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- polynomial-time detection of the unique maximal RAF (O(|R|^2 log|R|) worst case, sub-quadratic in practice)
- only LINEAR growth of the catalysis level with system size is needed for RAFs to appear: f(n) = 1.0970 + 0.0189 n at probability 0.5, i.e. 1-2 reactions catalysed per molecule up to n = 20, against Kauffman's exponential expectation
- sharp transition in the probability of a RAF with f: at n = 10 no RAFs at all below f = 1.20 and RAFs almost surely a little above it
- maxRAF size grows linearly with f while the size of an irreducible RAF stays roughly constant; a maxRAF can contain exponentially many irrRAFs

## Sources

- W. Hordijk and M. Steel. Detecting autocatalytic, self-sustaining sets in chemical reaction systems. Journal of Theoretical Biology 227(4):451-461, 2004. https://doi.org/10.1016/j.jtbi.2003.11.020 (paywalled; its definitions, algorithm and results are restated in the open sources below)
- W. Hordijk, S. A. Kauffman and M. Steel. Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. Int. J. Mol. Sci. 12(5):3085-3101, 2011. https://doi.org/10.3390/ijms12053085 - modified RAF definition (sec. 3.1), the RAF algorithm used here (sec. 3.2), Table 1 and Figures 2-3
- W. Hordijk, J. Hein and M. Steel. Autocatalytic sets and the origin of life. Entropy 12(7):1733-1742, 2010. https://doi.org/10.3390/e12071733 - the 1-2 reactions per molecule result
- M. Steel, W. Hordijk and J. Smith. Minimal autocatalytic networks. J. Theor. Biol. 332:96-107, 2013. arXiv:1212.4450 https://arxiv.org/abs/1212.4450 - maxRAF/subRAF/irrRAF definitions, the algorithm as a fixed point, and the n = 10 size data of Figure 4
- W. Hordijk, J. I. Smith and M. Steel. Algorithms for detecting and analysing autocatalytic sets. Algorithms for Molecular Biology 10:15, 2015. https://doi.org/10.1186/s13015-015-0042-8 - pseudo-RAFs, the irrRAF sampling, and the p values that give P_n = 0.5 for n = 8, 10, 12
- W. Hordijk and M. Steel. Autocatalytic sets extended: dynamics, inhibition, and a generalization. J. Syst. Chem. 3:5, 2012. arXiv:1206.1017 https://arxiv.org/abs/1206.1017 - the worked 4-reaction example and the u-RAF/inhibition results

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The 2011 (modified) definition is implemented: the RA condition asks for a catalyst in cl_R'(F) rather than in supp(R'), so a reaction catalysed only by a food molecule can be part of an RAF. Hordijk, Kauffman & Steel prove (Lemma 3.1) every RAF of the 2004 definition is one here.
- Cleavage and ligation form a single reversible reaction, as in the model papers, so the closure applies both directions and the F-generated test asks for every reactant AND product of a reaction to be in the closure.
- The level of catalysis f = p |R| is the parameter rather than p, since the published results are stated in f; p = f / |R| is recorded in extras.analysis. Catalysis events are drawn exactly: their number is Binomial(|X||R|, p) and the pairs a uniform subset of that size.
- Only catalysed reactions can belong to an RAF, so by default the network contains just the catalysed reactions (catalyst on both sides, one reaction per catalyst - the 'expanded CRS' of the papers); reactions='all' adds the uncatalysed cleavage/ligation pairs.
- The irrRAF search is the published deletion search (try to remove the reactions of the maxRAF in a random order, keeping any removal whose remaining maxRAF is non-empty); the 2004 paper that first states it is paywalled, so it was reconstructed from its description in the 2013 and 2015 papers.
- Inhibition (u-RAFs, NP-hard in general, fixed-parameter tractable in the number of inhibitors) is not implemented; the v1 parameter allow_inhibition is dropped rather than left as a flag that does nothing.
- No rates and no dilution: the framework is explicitly static. The food set is recorded in extras.food and, since it is assumed freely available, also in extras.buffered.

## Notes

An ANALYSIS to expose alongside chemical organisation theory: it runs on any catalytic reaction system (parameter `system`), including the network of `kauffman-autocatalytic-sets`. Note the scope limits the authors themselves state: no dynamics, no compartments, no heredity/selection.

---

*Specification: `catalog/chemistries/raf.yaml` · generator: `chemart/chemistries/raf.py` · tests: `tests/chemistries/test_raf.py`*
