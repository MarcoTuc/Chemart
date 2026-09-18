# Brane calculi

`brane-calculi` · *Cardelli, 2004*

*Also known as:* *Phago/Exo/Pino calculus*, *Mate/Bud/Drip calculus*, *Basic Brane Calculus*

Membranes as the active element rather than the container. A system is a nest of membranes, and each membrane carries a program of actions: phagocytosis engulfs a neighbour, exocytosis expels contents and fuses, pinocytosis buds an empty vesicle inward, mating fuses two membranes. A reaction is one such rearrangement of the nesting. The payoff is that viral entry and reproduction - engulf, fuse with the endosome, release, re-envelope - become short derivations rather than special cases.

| | |
|---|---|
| **family** | rewriting |
| **kind** | formalism |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.7 |
| **refs** | [157], [269], doi:10.1007/978-3-540-25974-9_24 |
| **provides** | `topology`, `stoichiometry`, `compartments`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): a system is a nest of membranes s[P] whose branes s are compositions of actions (s|t, !s, a.s), plus free molecules; ASCII notation with ',' for the composition of systems (paper ∘), [ ] for membranes, co-actions prefixed 'co' (paper ⊥) and p1(p2)=>q1(q2) for bind&release. A species is a whole configuration: the canonical text of the system normalised under structural congruence; its structure is the indented membrane tree

**R — reactions** (implicit, arity 1): one reduction step of the whole configuration: phago_n.s|s0[P], cophago_n(r).t|t0[Q] -> t|t0[r[s|s0[P]], Q]; coexo_n.t|t0[exo_n.s|s0[P], Q] -> P, s|s0|t|t0[Q]; pino(r).s|s0[P] -> s|s0[r[], P]; mate_n.s|s0[P], comate_n.t|t0[Q] -> s|s0|t|t0[P, Q]; cobud_n(r).t|t0[bud_n.s|s0[P], Q] -> r[s|s0[P]], t|t0[Q]; drip(r).s|s0[P] -> r[], s|s0[P]; p1, (p1(p2)=>q1(q2)).s|s0[p2, P] -> q1, s|s0[q2, P]

**A — reactor**: compartments
 · *dilution:* none; configurations can grow without bound (replicated actions), so the closure is truncated by max_species

## What you get

```python
net = chemart.generate_network("brane-calculi", seed=1)
```

```
brane-calculi: 4 species, 3 reactions, status=complete
provides: compartments, initial-state, stoichiometry, topology
seed: 1
extras: analysis, compartments, definitions, main, program, reaction_rules, system
```

First reactions:

```
!coexo|!cophago(mate)[!coexo|!comate[],Z],phago.exo[X|!bud[vRNA]] -> !coexo|!cophago(mate)[!coexo|!comate[],Z,mate[exo[X|!bud[vRNA]]]]
!coexo|!cophago(mate)[!coexo|!comate[],Z,mate[exo[X|!bud[vRNA]]]] -> !coexo|!cophago(mate)[!coexo|!comate[exo[X|!bud[vRNA]]],Z]
!coexo|!cophago(mate)[!coexo|!comate[exo[X|!bud[vRNA]]],Z] -> !coexo|!cophago(mate)[!coexo|!comate[],X|!bud[vRNA],Z]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `system` | `enum` | `viral-infection` | structural | published system of Cardelli 2004: viral-infection (sec. 3.3, fig. 7: virus, cell), viral-reproduction (sec. 3.3, fig. 8: membrane[nucap, envelope-vesicle, Z']), viral-replication (sec. 4.6, fig. 11: nucap, cytosol with vRNA replication, capsomer translation and ER), eat-me and seek-and-store (sec. 4.5), plant-vacuole (sec. 4.4, fig. 10, with ATP, Cl-, Na+ outside), custom (program) <br>one of `viral-infection`, `viral-reproduction`, `viral-replication`, `eat-me`, `seek-and-store`, `plant-vacuole`, `custom` |
| `program` | `str` | `` | structural | custom only: definitions 'name := term' (one per line, continuation lines allowed, # comments) that must define main, the initial system; syntax in the module docstring <br>*range:* e.g. 'A := =>n.phago[P] B := n=>.cophago(rho)[Q] main := A, B' |
| `max_species` | `int` | `1000` | structural | budget of configurations in the reduction closure; beyond it the network is truncated <br>`1` … `100000` · *range:* viral-infection closes at 4 configurations, viral-reproduction 4, eat-me 4, seek-and-store 5, plant-vacuole 5; viral-replication is infinite (vRNA replication) and truncated |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- viral infection (fig. 7): virus, cell -> Phago -> Mate with the endosome -> Exo, leaving membrane[nucap, cytosol]: the nucleocapsid reaches the cytosol
- viral reproduction (fig. 8): an envelope vesicle fuses by Exo with the cell membrane, which then buds the nucleocapsid out as a new virus
- nucleocapsid replication (sec. 4.6): disassembly releases vRNA, which is replicated, translated into capsomers that assemble new nucaps (B&R, Drip, B&R) and into envelope vesicles by the ER (B&R, Drip)
- Mate, Bud and Drip are derivable from Phago, Exo and Pino with the three-step encodings of sec. 3.2
- bitonality: Phago, Exo, Pino, Mate, Bud and Drip preserve the nesting parity of every subsystem
- molecular pumps and channels (sec. 4.4) and molecularly triggered membrane interactions (Eat Me, Seek and Store, sec. 4.5)

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 9.7 (one paragraph: brane calculi as process calculi on nested membranes; endocytosis, exocytosis, budding, mating; viral infections, molecular pumps, ion channels; Fellermann et al. [269]).
- Cardelli, L. (2005). Brane calculi: interactions of biological membranes. In V. Danos & V. Schachter (eds.), Computational Methods in Systems Biology (CMSB 2004), LNCS 3082, pp. 257-278, Springer (book [157]). Author's copy dated 2004-06-01: syntax and structural congruence (sec. 2.1), Phago/Exo/Pino (sec. 3.1, fig. 3), Mate/Bud/Drip and their encodings with derivations (sec. 3.2), viral infection part 1 with figs. 7-8 (sec. 3.3), bind&release (sec. 4.1), chemical reactions (4.2), plant vacuole (4.4), Eat Me and Seek and Store (4.5), viral infection part 2 (4.6). http://lucacardelli.name/Papers/Brane%20Calculi.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Source copy: the author's 2004-06-01 PDF was used; pdftotext loses every calculus glyph (∘ ◇ ⦇⦈ ⟲ ⊚ ⇉ ⊥), so all rules, definitions and derivations were read from rendered page images. The two exocytosis-like glyphs were told apart by the derivations themselves (virus = phago.exo[nucap]: Phago fires first in fig. 7, Exo last).
- [269] (Fellermann, Hadorn, Bönzli & Rasmussen, Living Machines 2012, LNCS 7375 pp. 343-344) is a two-page poster abstract with no open copy; its follow-up is a different formalism (the chemtainer calculus, Fellermann & Cardelli, J. R. Soc. Interface 11:20130987, 2014) and is not implemented here.
- Calculus offered: Phago/Exo/Pino, Mate/Bud/Drip as primitives (the paper says they should be primitives in practice and derives them from PEP; the derivations are reproduced in the tests), and bind&release of molecule multisets (sec. 4.1). Not offered: complexes m1:m2 (4.7), communication, choice, In/Out/Wrap (sec. 5 extensions), restriction.
- Species granularity: a species is a whole configuration and every step is a unimolecular reaction configuration -> configuration. Sub-membranes cannot be species because Phago, Mate and bind&release are contextual (they need a sibling or molecules outside the membrane) and every membrane interaction rebuilds the enclosing term, so reactions on parts do not compose into a multiset CRN.
- Structural congruence is decided by a normal form: compositions flattened and sorted by canonical text, units and 0[0] dropped, ! distributed over composition and !! collapsed, and a replicated element absorbs its plain copies (P,!P ≡ !P); replicated elements form a set. Reduction under ! happens only on an unfolded copy, which leaves the replicated term in place.
- Names: omitted subscripts match each other ('we omit them when there is no ambiguity'); pino and drip may carry a name, which has no effect. Metavariables of the paper (X, Z, Z', σ0, ρ, P, Q, Nucleus) are inert placeholders: identifiers on a brane that are not actions, or molecules in a system.
- Errata in the derivations of sec. 3.2 (page 6): Mate's second step prints ⟲⊥n''.τ|τ0⦇⟲n''|σ|σ0⦇◇⦈⦈ P ∘ Q⦈, whose extra bracket must be removed (coexo_n''.t|t0[exo_n''|s|s0[], P, Q], what the Exo rule gives); the Drip definition prints an unbalanced ⊚(⊚(ρ).⟲n)).⟲⊥n.σ, read as pino(pino(r).exo_n).coexo_n.s as its derivation requires.
- viral-reproduction starts from membrane[nucap, envelope-vesicle, Z'] with the section 3.3 definitions (nucap = !bud|X[vRNA]); after budding, the new virus can be phagocytosed again by the same membrane, a branch the paper does not discuss (it is stuck without an endosome in Z').
- viral-replication uses the section 4.6 definitions; vRNA-repl, written in the paper as the chemical reaction vRNA -> vRNA ∘ vRNA, is the section 4.2 encoding !vRNA=>vRNA,vRNA[] (an empty catalyst membrane). Its closure is infinite and truncated by max_species.
- plant-vacuole: the paper defines the vacuole membrane (proton pump, ion channel, proton antiporter) but no environment; Chemart puts one ATP, one Cl- and one Na+ outside an empty vacuole. Ion names are ASCII (Cl- for Cl–).
- Rates: none; the calculus is qualitative. max_species replaces any simulation size. extras: system, program, definitions, main, reaction_rules (the rule of each reaction: phago, exo, pino, mate, bud, drip, bind&release), compartments (membrane tree of the initial configuration) and analysis (the published configurations with whether they are reached, terminal configurations of complete closures, and for viral-replication the largest numbers of nucap and envelope-vesicle copies).

## Notes

Contrast with P systems: the membrane is the active element, not the container. The module exports parse_system, parse_brane, parse_program, steps (one-step reductions), reachable, layout, tree, membranes and molecule_parities for direct use.

---

*Specification: `catalog/chemistries/brane-calculi.yaml` · generator: `chemart/chemistries/brane_calculi.py` · tests: `tests/chemistries/test_brane_calculi.py`*
