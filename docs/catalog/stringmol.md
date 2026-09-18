# Stringmol

`stringmol` · *Hickinbotham, Clark, Stepney et al. (University of York), 2009-2011*

*Also known as:* *molecular microprograms*, *Stringmol automata chemistry*

An automata chemistry built to be run for a long time. Strings over a small alphabet bind by complementary template matching, and the bound pair then executes: one acts as enzyme, copying, cleaving or modifying the other, before dissociating. A seed replicase copies itself, and single-point mutations in a search template produce genuinely new behaviours rather than broken molecules - which is why it has been used to study open-ended evolution and parasite dynamics.

| | |
|---|---|
| **family** | bio-inspired |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 11.1.2 |
| **refs** | [379], [380], [381], [783], doi:10.1007/978-3-642-21283-3_37 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): Strings over 33 symbols: 26 template codes A..Z and 7 function codes $ > ^ ? = % }. Species id = the sequence with the function codes written as lowercase letters ($ s, > m, ^ t, ? i, = c, % x, } e); structure = the sequence. While bound, the active molecule also holds four pointers (instruction, flow, read, write), each with a toggle saying whether it points into its own string or the partner's.

**R — reactions** (implicit, arity 2): enzyme + substrate -> enzyme' + substrate' + cleaved products (one bind-execute-dissociate event)

**A — reactor**: well-stirred-multiset
 · *dilution:* constant per-step decay of every molecule; population held by the balance of energy influx and decay (about 350 molecules in ALife XII)

## What you get

```python
net = chemart.generate_network("stringmol", seed=1)
```

```
stringmol: 1 species, 1 reactions, status=observed
provides: catalysts, flow, initial-state, stoichiometry, topology
seed: 1
extras: aborted, active_counts, analysis, decayed, epochs, extinct, final_state, in_progress, time_steps
```

First reactions:

```
2 OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROOREsBLUBOtBmCscimssBLUBOxeOYHOB -> 3 OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROOREsBLUBOtBmCscimssBLUBOxeOYHOB  (x101)
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `container` | structural | container: the authors' time-stepped container with energy and decay, observed bind-execute-dissociate events with counts; soup: chemart.soup.soup with instantaneous pair reactions (bind test, program run to the end) and constant population; closure: all reactions reachable from the distinct seed molecules, with exact copying (Chemart additions) <br>one of `container`, `soup`, `closure` |
| `molecules` | `dict` | `{'OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>…` | population | initial multiset {sequence: count} (closure: the distinct sequences form the seed set) <br>*range:* seed replicase of spec v0.2 app. B.1 (default); upstream configs use WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB x 150; ALife XII species 9 is OBEQBXUUUDYGRHBBOSEOLHHHRLUEUOBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB |
| `steps` | `int` | `3000` | population | container: time steps; soup: collisions (elastic ones included); ignored by closure <br>`0` … `100000000` · *range:* ALife XII: until extinction, modal 750000 and up to about 15e6 time steps |
| `energy_per_step` | `int` | `25` | kinetic | energy units added to the container after every time step (container only); binding and each instruction cost one unit <br>≥ `0` · *range:* ALife XII paper: 25; upstream default ESTEP 20; later spatial configs 2500 |
| `cell_radius` | `float` | `2500.0` | spatial | container radius; bind propensity 1-(1-(agent_radius/cell_radius)^2)^n with n unbound molecules not yet visited in the step (container only) <br>≥ `0.001` · *range:* upstream CELLRAD 2500; spatial configs 1 |
| `agent_radius` | `float` | `10.0` | spatial | molecule radius in the bind propensity (container only); must not exceed cell_radius <br>≥ `0.0` · *range:* upstream AGRAD 10 |
| `decay` | `float` | `0.00023668639` | kinetic | probability per time step that a visited molecule is deleted; a decaying bound molecule takes its partner with it (container only) <br>`0.0` … `1.0` · *range:* spec v0.2 and ALife XII: 1/65^2; upstream configs 0.0005-0.0015 |
| `substitution_rate` | `float` | `1e-05` | stochastic | per-copy probability that '=' writes a neighbour of the read symbol in the symbol loop instead of the symbol <br>`0.0` … `1.0` · *range:* ALife XII: 1e-5; upstream configs (MUTATE) 1e-4-2e-3, with indel_rate equal to it |
| `indel_rate` | `float` | `3.06125e-08` | stochastic | per-copy probability of an insertion (the copy plus a random symbol) or a deletion (the symbol is skipped), 50:50 <br>`0.0` … `1.0` · *range:* upstream ALife XII default 3.06125e-8 (paper: p_s/(10 n) = 3.03e-8) |
| `max_length` | `int` | `2000` | structural | maximum string length: a copy with R or W at this index ends the reaction <br>`2` … `100000` · *range:* ALife XII and upstream MAXLEN: 2000; ECAL 2009: 512 |
| `traceback` | `enum` | `matrix` | structural | Smith-Waterman traceback: matrix is the proper trace matrix (spec v0.2, SmithWatermanV2); highest-neighbour walks back through the largest neighbouring cell (the v0.1 implementation, spec app. B.1) <br>one of `matrix`, `highest-neighbour` |
| `max_exec_steps` | `int` | `20000` | structural | soup and closure only: instructions after which an isolated reaction that has not ended counts as non-terminating (no reaction); the seed replicase needs a few hundred <br>≥ `1` |
| `max_species` | `int` | `100` | structural | closure only: species budget (status truncated when exceeded) <br>≥ `1` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- the seed replicase copies itself: R + R -> 3 R, with the active molecule copying the passive one (spec app. B.1; ALife XII: 240 time steps per copy)
- a single-point mutation in a search template ($BLUBO -> $BLUBP) makes a molecule that pastes a copy of its partner over its own end, giving a nearly double-length molecule (ALife XII fig. 7, species 29 + 9 -> 30; spec app. B.1)
- the double-length molecule makes truncated copies because the bind site shifts: 30 + 9 -> 30 + 9 + 31, where 31 lacks the first 14 symbols (ALife XII fig. 7)
- binding is complementary, so a replicase binds copies of itself only through two different sites; mutants with better-matching sites bind more strongly (ECAL 2009 invasion when rare: 88 of 100 trials)
- extinction: every one of 1000 ALife XII trials ends when a parasite that is copied but cannot copy spreads and the replicases decay (modal extinction time 750000 steps)
- characteristic sweeps (under 50000 steps), slow sweeps, neutral drift, sub-populations of non-lethal parasites, rapid sweep sequences
- emergent, spontaneous and multispecies hypercycles when a dominant species loses self-replication and depends on a partner (ALife XII: hypercycles in 30 of 1000 trials)

## Sources

- Hickinbotham, S., Clark, E., Stepney, S., Clarke, T., Nellis, A., Pay, M., Young, P. (2010). Specification of the Stringmol chemical programming language version 0.2. Technical Report YCS-2010-458, University of York: alphabet and complements (table 2), substitution matrix (eqs. 7-9), Smith-Waterman (eq. 10, app. B), bind probability (eq. 13), bind handshake (eq. 18), pointers and opcodes (secs. 7-10), mutation (eq. 29), decay (sec. 4, app. C), seed replicase (app. B.1). Shipped as smtr0.2.pdf in https://github.com/franticspider/stringmol
- Hickinbotham, S. et al. (2010). Diversity from a monoculture: effects of mutation-on-copy in a string-based artificial chemistry. ALife XII, MIT Press, 24-31: 25 energy units per step, about 350 molecules, p_s = 1e-5, p_i = p_s/(10n), seed replicase 65 symbols and 240 time steps per copy, species 9/29/30/31 and fig. 7 (origin of species 31). https://www-users.york.ac.uk/~ss44/bib/ss/nonstd/alife12_plazzmid.pdf
- Hickinbotham, S. et al. (2011). Molecular microprograms. ECAL 2009, LNCS 5777, 291-298: replicase, bind probability (s/l)^l of v0.1, length-dependent decay 1/L^2, invasion when rare. https://www-users.york.ac.uk/~ss44/bib/ss/nonstd/ecal09c.pdf
- Hickinbotham, S. stringmol, C++ source (GPL-3), commit 52b5625 (2021): stringPM.cpp (make_next, testbind, get_bprob, set_exec, exec_step, hcopy, cleave, rewind_bad_ptrs, testdecay), instructions.cpp (HSearch, IfLabel), alignment.cpp (SmithWatermanV2, SmithWaterman, default_table, sym_from_adj, align_event), stringmanip.cpp (AlphaComp), agents_base.cpp (eqn_prop), stringmol.cpp (SmPm_AlifeXII main loop). https://github.com/franticspider/stringmol
- Smith, T. F. & Waterman, M. S. (1981). Identification of common molecular subsequences. J. Mol. Biol. 147:195-197 (book [783]).

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The machine is a line-by-line port of the upstream C++ source, including its single-precision float arithmetic and the order in which random numbers are drawn. Compiled upstream with rand0to1 replaced by a shared LCG, it and the port agree on the substitution table, on 1503 alignments (score, start and end of both tracebacks, bind probability) and step by step (every molecule, pointer, toggle and the energy after every time step) on 7 container runs: ALife XII parameters (150 seed replicases, 3000 steps), dense binding with high mutation rates, random molecules with many function codes (3 seeds), max_length 130, and species 9 with mutant 29. Upstream crashes for max_length < 128 (HSearch clears 128 bytes of a max_length buffer); the port has no such limit.
- Where the spec (and the book figure copied from it) and the source disagree, the source wins, since it produced the published runs: (1) symbol order of the substitution matrix and mutation loop is ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ (default_table, ALXII.mtx), not the spec's ABC$DEF>GHIJ^KLM=NOP?QRS}TUV%WXYZ; complements are the same in both; (2) mismatch scores are -d*33/272 printed to 3 decimals (eq. 8 as written); spec table 3 (-0.12, -0.25, ...) comes from an integer division in table_from_string, which the source does not use; (3) P(bind) is 0 for aligned length l <= 2 (spec: l <= 3) and uses m = 1.124 in eq. 13; (4) the bind propensity uses radii, (agent_radius/cell_radius)^2 = 1.6e-5 with CELLRAD 2500 and AGRAD 10, where spec eq. 1 uses the areas 10/2500, and n counts the unbound molecules not yet visited in the current step.
- Spec erratum: eq. 18 makes M_A active when its bind site starts nearer the string start, contradicting the prose (the executing string is the one whose bind site is furthest from its start); the source (set_exec) and the ALife XII paper follow the prose, and on ties the molecule that initiated the bind is active.
- Opcode details from the source that the spec states differently or not at all: '$' succeeds with probability min(score, l-1.124)/(l-1.124) with no cutoff for short alignments (the cutoff line is overwritten), puts F just after the match on the searched string, and on failure puts F on the last template symbol of the label on I's string (spec: the next instruction); a '$' with no template sets F to I; '>' moves I to F + 1 (spec: to F, or +1 if already there); '^' increments I after toggling, so '^A' increments the pointer on the other string; '?' with no or a one-symbol template skips one instruction if R is past the end, otherwise it skips with probability (score/len)^len of the template against the symbols at R (align_event), not eq. 13.
- Copy: '=' overwrites the symbol at W (the ECAL paper says it inserts; W is normally at the end of the string, so the string grows). A read past the end writes nothing but I still advances. R or W at max_length ends the reaction. Mutation (hcopy): with probability indel_rate an indel, then 50:50 an insertion that writes the copy and then a random symbol (two symbols) or a deletion that skips the read symbol and also advances I once more; otherwise with probability substitution_rate a random neighbour of the read symbol in the symbol loop. indel_rate defaults to the source's ALife XII value 3.06125e-8, the paper's formula p_s/(10*33) gives 3.03e-8.
- Cleave splits the string that F is on (spec: 'the' string); if F is inside it, the tail becomes a new molecule, pointers past the new end move to the end, and if either string is left empty that molecule is destroyed and the reaction ends (cleave-zero). Decay: every visited molecule decays with constant probability per step whether bound or not (spec v0.2 DC+ALL); only the visited member of a pair is tested, and its partner is deleted with it.
- Energy: nothing but decay happens while energy <= 0; each successful bind and each instruction (including '}') costs one unit; the container starts with 0 energy and energy_per_step is added after each step (SmPm_AlifeXII). The ALife XII paper says 25 units per step; the source's default ESTEP is 20; 25 is the default.
- Seed replicase: the ALife XII figure's 65-symbol seed is only an image. The default is the 64-symbol seed that spec app. B.1 gives for the matrix traceback (OOGEOL...). ALife XII species 9, printed in the paper's fig. 7, is kept as a constant for the tests. The v1 catalog parameters complement_offset (fixed at 13 in AlphaComp), bind_probability (it follows from the alignment) and capacity (the population follows from energy and decay) are dropped; seed_replicase becomes molecules.
- Traceback option: spec app. B.1 says v0.1 walked back through the highest neighbouring cell and v0.2 uses the trace matrix; both are in the source and both are offered. The ALife XII fig. 7 reactions (29 + 9 -> 30 + 9 and 30 + 9 -> 30 + 9 + 31) come out the same with both.
- Observed network (container): one reaction per complete event, reactants = the two sequences at binding, products = the two sequences at dissociation plus every molecule cleaved off during the event; reactions are deduplicated as multisets, with firing counts; extras.active_counts gives, per reaction, how often each reactant species was the active one. Events cut short by decay are not reactions: extras.aborted lists them (reactants, molecules already released, count); extras.decayed counts deletions of unbound molecules; extras.in_progress lists complexes still bound at the end; extras.final_state counts the unbound molecules. Then initial_state + sum(count * (products - reactants)) + sum over aborted and in_progress of (released - reactants) - decayed = final_state exactly. outflow is the per-step decay probability (a first-order removal of every species). extras.analysis samples population, distinct species and energy (about 500 samples); extras.epochs lists when the most abundant species changes (sweeps).
- Rates: none. Binding probability depends on the alignment and reaction duration on the program and energy, so no mass-action constant is published or implied.
- Soup and closure are Chemart additions for sampling the reaction space without the container: an isolated reaction runs the pair's program to the end with unlimited energy and no decay; the first molecule initiates the bind (soup draws pairs in random order). Soup tests the bind with its probability and mutates on copy; closure takes every bind with P > 0 and copies exactly, but still draws '$' and '?' outcomes at random once per ordered pair. A reaction still running after max_exec_steps instructions counts as non-terminating (no reaction; counted in extras.non_terminating).
- Undefined behaviour in the source, fixed choices: an alignment against an empty string (the C result is uninitialised) scores 0; reads beyond the string buffer return the terminator.
- The v1 reactor lattice-2d is removed: the book calls the topology dispensable, [379]-[381] use the aspatial container, and the later spatial Stringmol (smspatial.cpp, GRIDX/GRIDY) is not implemented.

## Notes

Paper-scale runs (about 350 molecules for 10^5-10^7 time steps) are far beyond pure Python; the defaults run a few thousand steps with 100 seed replicases, where mutations are rare. Raise substitution_rate (the upstream configs use 1e-4 to 2e-3) to see new species within seconds.

---

*Specification: `catalog/chemistries/stringmol.yaml` · generator: `chemart/chemistries/stringmol.py` · tests: `tests/chemistries/test_stringmol.py`*
