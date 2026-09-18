# Machine-tape chemistry

`ikegami-hashimoto` · *Ikegami & Hashimoto, 1995*

*Also known as:* *machine-tape interaction*, *coevolution of machines and tapes*

Von Neumann's separation of constructor and description, made into a chemistry. Two populations coexist: machines, which act, and tapes, which are read. A machine reading a tape produces a new machine and a new tape, so reproduction requires both. A machine whose description is missing washes out of the reactor, which puts pressure on closed loops of machines that between them produce each other's descriptions - including the minimal case of a machine and tape that reproduce as a pair.

| | |
|---|---|
| **family** | automata |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 10.5.4 |
| **refs** | [419], [420], doi:10.1007/3-540-59496-5_302 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): Two separate populations. Tapes: circular 7-bit strings, species T<2 hex digits> = the bits read clockwise from the tape's source (the papers name a tape by its smallest rotation, extras.paper_names). Machines: 16-bit words, species M<4 hex digits> = T' column, M' column, head (4 bits), tail (4 bits) of a transition table (tape bit, machine state) -> (tape bit', state').

**R — reactions** (implicit, arity 2): M + T -> M + T + M' + T'

**A — reactor**: ode, well-stirred-multiset
 · *dilution:* each population has capacity N; every generation a fraction c of each is replaced by reaction products (c = d = 0.6)

## What you get

```python
net = chemart.generate_network("ikegami-hashimoto", seed=1)
```

```
ikegami-hashimoto: 164 species, 1482 reactions, status=observed
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: analysis, final_state, noise_induced, paper_names
```

First reactions:

```
M7922 + T22 -> M7922 + 2 T22 + M0aa1  (x2)
M7922 + T22 -> M7922 + T22 + M2ee5 + T72  (x2)
M7922 + T22 -> M7922 + 2 T22 + M2885  (x2)
M7922 + T22 -> M7922 + T22 + M0cc1 + T14  (x2)
M7922 + T22 -> M7922 + T22 + M0881 + T02  (x1)
M7922 + T22 -> M7922 + T22 + M2aa5 + T49  [mass-action k=0.3 frame_length=3]  (x5)
M7922 + T22 -> M7922 + T22 + M0ee1 + T52  [mass-action k=0.3 frame_length=3]  (x5)
M24e7 + T6f -> M24e7 + T6f + M2ee5 + T72  (x2)
… and 1474 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `dynamics` | structural | dynamics: the papers' population dynamics, observed reactions with firing counts; closure: the noise-free reaction network reachable from the seed machines and tapes (Chemart addition) <br>one of `dynamics`, `closure` |
| `machines` | `list` | `` | structural | seed machines as 16-bit hex strings (T' column, M' column, head, tail); empty draws n_machines random machines <br>*range:* e.g. the minimal self-replicating loop [1002] with tapes [01] |
| `tapes` | `list` | `` | structural | seed tapes as 7-bit hex strings read from the source; empty draws n_tapes random tapes <br>*range:* e.g. [01]; the papers' T1 read from another source is 04, 10, ... |
| `n_machines` | `int` | `10` | population | number of random seed machines when machines is empty <br>`0` … `65536` · *range:* papers: about 10 |
| `n_tapes` | `int` | `3` | population | number of random seed tapes when tapes is empty <br>`0` … `128` · *range:* papers: 2 or 3 |
| `N` | `int` | `1000` | population | dynamics only: capacity of each population (machines, tapes); seeds share it equally <br>`1` … `10000000` · *range:* not stated in the papers |
| `c` | `float` | `0.6` | kinetic | fraction of each population replaced by reaction products per generation (c = d_m = d_t); also the rate scale of eq. 4 <br>`0.0` … `1.0` · *range:* papers: c = d = 0.6 |
| `noise` | `float` | `0.05` | stochastic | dynamics only: external noise mu_P, bit-flip probability per bit of the reading frame <br>`0.0` … `1.0` · *range:* papers: 0.04 (minimal loop), 0.055 (oscillation), 0.07-0.08 (core networks), scans 0-0.1 |
| `generations` | `int` | `150` | population | dynamics only: number of generations <br>`0` … `100000` · *range:* papers: 1600-3000 |
| `noise_off` | `int` | `-1` | stochastic | dynamics only: generation at which external noise is turned off; -1 keeps it on <br>`-1` … `100000` · *range:* papers: 2000 |
| `source` | `enum` | `random` | stochastic | dynamics only: source of a tape whose circular pattern is new to the population: random site (papers) or the site it was written from <br>one of `random`, `inherit` |
| `max_species` | `int` | `400` | structural | closure only: species budget (the noise-free closure has at most 128 tapes and 128 + seed machines) <br>`2` … `100000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- a machine without its description tape is washed out; a closed loop of machines producing each other's tapes is needed
- low noise: minimal self-replicating loop M1002 + T1 (Eigen-Schuster type, zero active mutation), metastable: it survives when noise is turned off
- higher noise (0.04-0.055): parasitic machines (M3006 with T5, M1222 with T3) invade and populations oscillate, with bursts of active mutation
- mid noise (0.05-0.1), after about 600-1100 generations: a core network of about 40-60 machines with high active mutation that keeps its diversity after noise is turned off (fixed point or oscillating); double autocatalytic loops (RNA-like editing) besides Eigen-Schuster loops (DNA-like replication)
- high noise: back to a minimal loop or extinction
- cell assemblies of core networks differentiate ([420], not implemented)

## Sources

- Ikegami, T. & Hashimoto, T. (1995). Active mutation in self-reproducing networks of machines and tapes. Artificial Life 2(3):305-318. Author's manuscript dated May 31, 1996: model (sec. 2, eqs. 2-9, fig. 1), results (figs. 2-9, sec. 4 list of self-replicating pairs). https://drive.google.com/file/d/1WwFXhB10rx8FiTPZZjUVpQIHl6BtN73D/view (linked from https://www.sacral.c.u-tokyo.ac.jp/publications)
- Ikegami, T. & Hashimoto, T. (1995). Coevolution of machines and tapes. In Moran et al. (eds.), Advances in Artificial Life (ECAL 95), LNAI 929, 234-245. Eq. 4 with c_ij = 0 for non-reading pairs, source rule, figs. 1-6. https://web.archive.org/web/20170809083155/http://sacral.c.u-tokyo.ac.jp/pdf/ikegami_alife_1995.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Bit layout from paper fig. 1, not from the prose. The prose says head and tail come first and translation starts at the reading frame; the figure translates the rewritten tape [1110110] from its source into M ebbd with the table columns in bits 0-7 (T', M' alternating) and head/tail in bits 8-15 (alternating). That layout reproduces fig. 1 exactly, gives the machines of the published pairs (M1002 <- T1, M3006 <- T5, M1222 <- T3), makes the five self-replicating pairs of sec. 4 self-replicate, and makes M1002 the machine produced by the largest number of machines (sec. 5).
- The rewrite runs over sites h .. tail-1 (fig. 1: head at site 3, tail at site 2, six steps), so L is 1-6 (the <L> axis of figs. 2-3 ends at 6). The tail is searched from h + 1 and may overlap the head.
- Tapes are strings read from their source, so the 7 rotations of a circular tape are different species (one tape encodes different machines). In the dynamics, rotations share one source ('every translational invariant tape has the same source'): a product whose circular pattern is present is merged into that rotation; a new pattern gets a random source (source: random) as in the papers. Translation of the product machine is from the source of the tape that was read, as in fig. 1. The closure has no population, so it keeps rotations apart.
- Normalisation: the Artificial Life text makes c N new machines and tapes per generation, the book says the population is kept constant, while ECAL eq. 4 divides by the sum over all pairs (so non-reading pairs waste production). Chemart follows the book and the Artificial Life text: the c N products are shared among reading pairs, d = c, and populations stay at N (up to integer parts). Each network reaction gets mass-action k = c/2 per initial state producing it (c when both states give the same products), and outflow constant-total: the rate equations are the continuous-time version of eqs. 4-6 without noise and integer parts, up to a rescaling of time by the sum of f_M f_T over reading pairs (checked in tests).
- Noise: the reaction terms are multiplied by 1 - eps with eps = 1 - (1 - mu_P)^L (eq. 7), and eps times the flux of each reaction and state becomes mutants, rounded stochastically to an integer ('at most eps c N mutant populations', Monte Carlo). A mutant flips each frame bit of the rewritten tape with probability mu_P conditioned on at least one flip, and its machine is translated from the mutated tape. The papers do not say whether the error hits the read or the written bit; Chemart flips the written tape.
- Integer parts: populations are floored after each generation (so objects below one copy, f < 1/N, are removed); mutants are whole objects. Seeds share N equally; the papers do not give initial populations or N (default 1000).
- Observed network: every reaction that fired, with count = number of generations in which it fired. Error-free reactions carry the mass-action rate above and frame_length L (the noise factor (1 - mu_P)^L is not folded in, since noise_off changes it); reactions produced only by noise have no rate and are listed by index in extras.noise_induced. extras.analysis records per generation the distinct machines and tapes (figs. 4-6), the active mutation rate <mu_A> of eq. 9 (w/L averaged over reading pairs with weights m_i t_j and over both initial states) and the average frame length <L>.
- Dropped v1 params: tape_len and machine_len (the encoding is only defined for 7-bit tapes and 16-bit machines); compartments (the cell model of [420] is not implemented: the ALife V paper is not openly available).
- Paper errata: sec. 4 lists the self-replicating pair 'Mbdd1 with T37', but bdd1's tail 0001 cannot bind T37 (two zeros); Mbdd7 is the self-replicator of T37. It also writes 'T9dd3' for M9dd3. Fig. 1 calls the rewritten tape [1110110] Te6; its value is 76 hex (T37 as a circular tape). ECAL fig. 1 names M1222's tape T41, a rotation of T3.
- Sec. 4 says there are 5 possible self-replicating pairs (M1002/T1, M2004/T1, Mdffb/T3f, M9dd3/T1d, Mbdd7/T37). With this encoding all five self-replicate from initial state 1, but so do 10 more (e.g. M0000/T0, Mffff/T7f, M3446/Td); the papers do not say how the five were selected.

## Notes

A von Neumann style replicator chemistry with separate description (tapes) and constructor (machines). Machines are only 16 bits and tapes 7 bits, so every reaction is a table lookup; the noise-free closure of any seed has at most 128 tapes. The papers' core-network phenomena need long noisy runs (thousands of generations) with a capacity N they do not state.

---

*Specification: `catalog/chemistries/ikegami-hashimoto.yaml` · generator: `chemart/chemistries/ikegami_hashimoto.py` · tests: `tests/chemistries/test_ikegami_hashimoto.py`*
