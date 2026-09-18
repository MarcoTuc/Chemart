# HBCB / PSD degradation-and-reuse chemistry

`hbcb-psd` · *Oohashi, Ueno, Maekawa, Kawai, Nishina & Honda, 2009; PSD model proposed by Oohashi et al. 1987*

*Also known as:* *hierarchical biomolecular covalent bond model*, *programmed self-decomposition model*

A question about decomposition, posed as a chemistry. Biomolecules are arranged in a hierarchy - biopolymers at the top, monomers below them, then progressively smaller classes - and hydrolysis walks down while synthesis walks up, each step with an energy cost. The model asks how deep an organism should break down its own dead material before rebuilding. Stopping at the monomer leaves exactly one synthesis step to get back, which is why programmed decomposition to that level pays off in a closed, finite world.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.3.2 |
| **refs** | [643], doi:10.1162/artl.2009.15.1.15103 |
| **provides** | `topology`, `stoichiometry`, `energies`, `mass-conservation`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (implicit): One class per level of the HBCB hierarchy, per family of biomolecules: the biopolymer BP at the top, the biomonomer BM directly below it, then the lower classes C(n-3) ... C0 down to the base class. A class-l molecule is `units_per_class` molecules of the class below, so a species is fixed by (family, level); its id is '<class>_<family>' (BP_0, BM_0, C1_0, C0_0) and its structure records the composition and the number of base units it contains.

**R — reactions** (implicit, arity [1, 4]): Hydrolysis (cleavage) X_l -> units_per_class X_(l-1) and its reverse, synthesis units_per_class X_(l-1) -> X_l, between every pair of neighbouring classes; plus one programmed self-decomposition per family, BP -> units_per_class^d X_t, which cleaves a biopolymer in one genetically regulated step down to the class named by psd_target (with psd_target = BM this is the single hydrolysis step the model calls the main PSD process).

**A — reactor**: lattice-2d, well-stirred-multiset
 · *dilution:* Closed ecosystem: matter is finite and recycled, so there is no inflow or outflow and the base units of every family are conserved exactly.

## What you get

```python
net = chemart.generate_network("hbcb-psd", seed=1)
```

```
hbcb-psd: 16 species, 24 reactions, status=complete
provides: energies, initial-state, mass-conservation, stoichiometry, topology
seed: 1
extras: analysis, conservation, energies, hbcb, psd
```

First reactions:

```
BP_0 -> 4 BM_0
4 BM_0 -> BP_0
BM_0 -> 4 C1_0
4 C1_0 -> BM_0
C1_0 -> 4 C0_0
4 C0_0 -> C1_0
BP_1 -> 4 BM_1
4 BM_1 -> BP_1
… and 16 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n_levels` | `int` | `4` | structural | number of classes in the HBCB hierarchy, from the base class up to the biopolymer; the top class is BP, the one below it BM, and the rest are the lower classes C(n-3) ... C0 <br>`2` … `10` |
| `units_per_class` | `int` | `4` | structural | molecules of class l-1 in one molecule of class l; a BP therefore holds units_per_class^(n_levels-1) base units <br>`2` … `8` |
| `n_families` | `int` | `4` | structural | independent families of biomolecules, each with its own hierarchy (proteins, nucleic acids, ...); families do not interconvert <br>`1` … `8` · *range:* the SIVA ecosystem of the same group (Sayama 1998, fig. 3 and table 1) holds four finite substances; the number used in the 2009 paper is not known |
| `psd_target` | `enum` | `BM` | selection | depth of the programmed self-decomposition: 'BM' cleaves a BP to biomonomers, which can be re-used directly (the species the paper finds superior); 'sub-BM' one class below BM and 'base' all the way to the base class, both of which have to be rebuilt through every intermediate class; 'none' is a species that does not self-decompose <br>one of `none`, `BM`, `sub-BM`, `base` |
| `bond_energies` | `list` | `[4.0, 2.0, 1.0]` | thermodynamic | energy of one covalent bond at each level of the hierarchy, entry i being the bonds that link class-i units into a class-(i+1) molecule; needs n_levels - 1 entries and must decrease strictly, the HBCB ordering (deeper bonds are stronger, the BM-BM links of a BP weakest). The ordering is the model; the numbers are a Chemart default |
| `initial_bp` | `int` | `8` | population | biopolymers of each family present at the start, the matter of the closed ecosystem (initial_state); the paper gives no amounts <br>`0` … `1000000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- organisms that decompose their biopolymers to reusable biomonomers have an evolutionary advantage over organisms that decompose them into classes below BM
- decomposing to BM leaves exactly one synthesis to return to a biopolymer, while any deeper class needs strictly more synthesis steps and strictly more covalent bond energy
- the hierarchy is ordered by covalent bond energy: the BM-BM links of a biopolymer are the weakest bonds and hydrolysis cleaves them first
- matter is finite and recycled: the base units of every family are conserved by hydrolysis, synthesis and programmed self-decomposition alike

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, MIT Press, section 18.3.2, p. 388: the paragraph on Oohashi et al.'s HBCB and PSD models and their virtual ecosystem.
- Oohashi, T., Ueno, O., Maekawa, T., Kawai, N., Nishina, E. & Honda, M. (2009). An effective hierarchical model for the biomolecular covalent bond: an approach integrating artificial chemistry and an actual terrestrial life system. Artificial Life 15(1):29-58. https://doi.org/10.1162/artl.2009.15.1.15103 - NOT OBTAINED (see decisions); only the published abstract (via Crossref and Europe PMC, PMID 18855570) was available.
- Sayama, H. (1998). Development of a highly-parallel simulator of individuals of virtual automata for verification of the 'programmed self-decomposition model'. IPSJ Transactions 39(6):1782-1789 (in Japanese), open access: https://ipsj.ixsq.nii.ac.jp/record/13055/files/IPSJ-JNL3906025.pdf - the SIVA-2r ecosystem: Oohashi et al.'s self-reproduction~self-decomposition (SRSD) system as von Neumann's automaton with a decomposer module FZ (figs. 1-2), the 256x256 lattice with a temperature gradient and four finite substances, the chromosome genes (table 1), and the result that self-decomposing populations keep adapting while non-decomposing ones stagnate once space and matter are used up (figs. 10-11, 2 x 13 x 10 = 260 runs over 13 mutation rates).

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The 2009 paper (book ref [643]) could not be obtained. It is bronze open access on MIT Press, but direct.mit.edu answers both curl and WebFetch with a Cloudflare challenge / HTTP 403; the Internet Archive (Wayback and its CDX API) was globally offline during this session; Springer, ingenta and every other indexed location are closed access; Unpaywall, OpenAlex, Semantic Scholar, Europe PMC, CORE, OpenAIRE, fatcat and CiNii list no other full text. The implementation therefore follows book 18.3.2 and the paper's published abstract, which is why fidelity is book+decisions and not reconstructed. Everything below that the abstract and the book do not state is a Chemart reconstruction.
- From the abstract and the book, the model's content is: biomolecules are arranged hierarchically by the energy of their covalent bonds (HBCB); self-decomposition is a genetically regulated, endergonic process whose main step is hydrolysis, which decomposes a BP into BMs; and virtual species whose self-decomposition cleaves a BP to BMs are evolutionarily superior to species that cleave it into classes lower than BM. The depth of decomposition is therefore the parameter the model compares (psd_target), and the advantage is a re-use cost, which is what extras.analysis reports.
- The number of classes, the branching factor units_per_class, and which concrete biomolecules sit at each class are not known from the abstract; they are parameters (n_levels, units_per_class), with only BP and BM named as the model names them. The lower classes are written C(n-3) ... C0 rather than given chemical identities.
- Bond energies: the HBCB ordering (bonds deeper in the hierarchy are stronger, the BM-BM links of a biopolymer weakest, so hydrolysis cleaves them first) is the model and is enforced - bond_energies must decrease strictly upwards or generation raises ValueError. The numeric default [4.0, 2.0, 1.0] is a Chemart choice and is not a published value. A molecule's total bond energy (extras.energies) is accumulated down the hierarchy, a class-l molecule being units_per_class class-(l-1) molecules joined by units_per_class - 1 bonds; that arithmetic is also a Chemart convention.
- Each class-l molecule is treated as exactly units_per_class molecules of the class below, and cleavage as complete (a molecule falls apart into all its subunits at once, matching 'hydrolysis decomposes a BP into BMs'). Partial oligomers are not modelled: the paper's cleavage pattern is not known, and inventing one would invent chemistry.
- Programmed self-decomposition is one reaction per family, BP -> units_per_class^d X_t, the organism's single genetically regulated process; the stepwise cleavages of the hierarchy are present as the chemistry in their own right. With psd_target = 'BM' the two coincide, and the PSD reaction is that same hydrolysis reaction rather than a duplicate.
- No kinetics are published, so every rate is None and the entry provides no rate-constants; the fact that PSD is endergonic and genetically regulated is recorded in the notes rather than turned into a rate law or an energy barrier.
- The ecosystem is closed (no inflow or outflow): the book calls it a virtual ecosystem with recycling and the SIVA lineage gives it finite matter, so the base units of each family are conserved exactly (extras.conservation, one law per family, the vector being units_per_class^level). initial_bp is a Chemart amount; the paper gives none.
- Selection between decomposition strategies acts across competing organisms and cannot live inside one reaction network, so it is not generated. The cost of returning to a biopolymer from each class (number of syntheses and covalent bond energy per BP) is reported in extras.analysis.rebuild_cost_by_class instead, and that is the quantity in which decomposing to BM is optimal. The v1 parameters bond_hierarchy (type matrix) and decomposition_strategy become bond_energies and psd_target.
- v1 said constructive: true. The species set is finite, closed and known from the parameters (|S| = n_families x n_levels: the closure of the cleavage rule terminates at the base class), so it is corrected to false.
- n_families defaults to 4 because the SIVA ecosystem of the same group holds four finite substances (Sayama 1998); the 2009 paper's number is unknown, and the families do not interact, so this only replicates the hierarchy.

## Notes

A thermodynamic-hierarchy chemistry rather than a constructive one: the interest is in which level of the hierarchy an organism decomposes its own matter down to, and in what that costs to undo. The book's paragraph is four sentences long and the original paper could not be obtained, so the hierarchy's size and its bond energies are parameters and the published result is reproduced as a re-use cost rather than as an evolutionary run. See `smn` and `bnc-cell` for the neighbouring string-based metabolic chemistries of book 18.3.

---

*Specification: `catalog/chemistries/hbcb-psd.yaml` · generator: `chemart/chemistries/hbcb_psd.py` · tests: `tests/chemistries/test_hbcb_psd.py`*
