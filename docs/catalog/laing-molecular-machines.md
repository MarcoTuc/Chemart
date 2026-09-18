# Laing's artificial molecular machines

`laing-molecular-machines` · *Laing, 1972-1977*

*Also known as:* *Laing molecular tapeworms*, *hybrid cellular-kinematic automata*

The first molecular machine chemistry, from before anyone simulated such things. Strings come in two forms: passive tapes of bits, and active machines - folded chains of instructions where the fold is what brings distant instructions together. A machine attaches to a tape, runs as a read-head over it, and leaves the tape cut into pieces. Laing's goal was self-reproduction by *self-inspection*: a machine that reads its own structure rather than copying a stored description of it.

| | |
|---|---|
| **family** | automata |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 10.5.1 |
| **refs** | [484], [485], [486], [487], [488], doi:10.1016/0022-5193(77)90294-6 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): strings of constituents in two forms: passive binary tapes 't:0110' and active machines 'm:CT1.W1.H.TT1' (a folded chain of instructions; the fold joining CT<k> to TT<k> is the shared label). Species structure is 'passive' or 'active: ' followed by the instructions with the index each CT jumps to

**R — reactions** (implicit, arity 2): m + t -> m + t_1 + ... + t_n: the machine attaches at its first instruction to one unit of the tape, runs until it halts, and the tape is replaced by its pieces; machine + machine and tape + tape do not react

**A — reactor**: sequential-vm, well-stirred-multiset
 · *dilution:* closure: none; soup: a random molecule is removed after each productive event, keeping the population size constant (Chemart addition)

## What you get

```python
net = chemart.generate_network("laing-molecular-machines", seed=1)
```

```
laing-molecular-machines: 17 species, 15 reactions, status=truncated
provides: catalysts, stoichiometry, topology
seed: 1
extras: seed
```

First reactions:

```
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:0 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:1
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:1 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:01
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:01 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:11
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:11 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:001
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:001 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:101
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:101 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:011
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:011 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:111
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:111 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:0001
… and 7 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `closure` | structural | closure: every reaction reachable from the seed machines and tapes (status complete or truncated); soup: random pairs drawn from a population with chemart.soup, returning the reactions that fired with counts (status observed) <br>one of `closure`, `soup` |
| `machines` | `list` | `['CT1.W1.H.TT1.W0.R.CT1.W1']` | structural | the active machine strings, instructions separated by '.'; every CT<k> needs exactly one TT<k>. The default adds 1 to a little-endian binary number (Chemart example) <br>*range:* any program over W0, W1, L, R, H, NOP, D, CT<k>, TT<k>; turing_program() compiles a 2-symbol Turing machine, e.g. the busy beavers of the tests |
| `tapes` | `list` | `['0']` | structural | the passive seed tapes, strings over 0 and 1 |
| `binding` | `enum` | `leftmost` | stochastic | the tape unit the machine attaches to: leftmost, every unit (each outcome a separate reaction; closure only), or one drawn with the generator (per pair in the closure, per collision in the soup) <br>one of `leftmost`, `all`, `random` |
| `max_length` | `int` | `4` | structural | longest tape kept; a run that grows a tape beyond it gives no reaction and the status becomes truncated <br>`1` … `10000` |
| `max_steps` | `int` | `1000` | structural | instructions executed before a run is taken as non-halting; it then gives no reaction and the status becomes truncated <br>`1` … `10000000` |
| `max_species` | `int` | `200` | structural | closure only: species budget; the status becomes truncated when it cuts the closure off <br>`1` … `100000` |
| `copies` | `int` | `20` | population | soup only: copies of each seed machine and tape in the initial population, whose size the soup keeps constant <br>`1` … `100000` |
| `steps` | `int` | `500` | stochastic | soup only: number of pair draws <br>`0` … `10000000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- an active string plays a Turing machine read-head and a passive string its tape, so the system can carry out any Turing computation (Laing [485], quoted in Sipper 1998); here: the 2-state busy beaver writes 4 ones in 6 moves, the 3-state machines 6 ones in 14 moves and 5 ones in 21 moves (Rado 1962; Lin & Rado 1965)
- a tape grows when the machine slides off an end and a new 0 is recruited (Freitas & Merkle 4.8)
- detachment cuts a tape so that one reaction yields several product molecules
- self-reproduction by self-inspection and self-description (Laing [486-488]; not reproduced, see decisions)

## Sources

- Banzhaf & Yamamoto (2015), 10.5.1 and fig. 10.5 (reprinted from [486]): two forms of molecules, passive data and folded active machines; attachment at one position of each; the instruction at the contact runs and moves the contact on the tape (left, right), reads or writes, and the machine's contact advances to the next instruction; instruction set w(1), w(0), TT, CT, L, R, NOP.
- Freitas, R. A. Jr. & Merkle, R. C. (2004). Kinematic Self-Replicating Machines, section 4.8 'Laing Molecular Tapeworms (1974-1978)', quoting Laing's rigid-constituent version: tape of 0/1 molecules; ten instruction types W(0), W(1), L, R, H, NOP, CT (origination conditional transfer), TT (destination transfer target), detachment and synthesize; sliding off the tape recruits a new 0; the self-reproduction procedure. http://www.molecularassembler.com/KSRM/4.8.htm
- Sipper, M. (1998). Fifty years of research on self-replication: an overview. Artificial Life 4(3):237-257, section 4 and fig. 12 (Laing 1977 quoted: strings of primitive finite-state automata in sliding contact, active and passive primitives; 1977 primitives N, 0, 1, P0, P1, F, B, H, TN, T0, T1, A, AD, C). https://fab.cba.mit.edu/classes/865.18/replication/Sipper.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Laing's papers [484-488] were not accessible: Elsevier and J. Cybernetics are paywalled, the University of Michigan Deep Blue copy of [488] and Laing's 1977 Binghamton dissertation 'Automaton self-reference' both returned 403 to automated downloads, and the Wayback Machine was offline. The machine is built from the book and the two secondary sources above, keeping only what they state; every other rule is listed here.
- Instruction set: the book's fig. 10.5 set (W0, W1, L, R, NOP, CT, TT) plus H and D (detach) from Freitas & Merkle's list of Laing's 1975 constituents. Their 'synthesize' (turns a passive parts blank into a constituent) and Laing's 1977 activation, activate-and-detach and conversion primitives are not implemented, because no accessible source says which state or type they produce. So machines are never built: only tapes are created, and Laing's self-reproduction by self-description and self-inspection ([486-488]) is not reproduced.
- CT/TT: the fold that brings a CT into contact with its TT is written as a shared label (CT1 ... TT1); several CTs may share one TT, and each label has exactly one TT. A CT transfers when the contacted tape constituent is 1 and falls through on 0, as in Wang's B-machine conditional transfer (the book names Wang's machine as the model). Execution continues after the TT, which is a no-op.
- Attachment: the machine attaches at its first instruction (the book: the instruction at the binding position runs, then the contact advances to the next instruction). The tape unit is a parameter (binding), since the book only says the two molecules touch at one position each.
- Halting: H or running past the last instruction. Laing's machines need not halt; a run still going after max_steps instructions gives no reaction and makes the closure truncated.
- D (detach) severs the contacted constituent from its predecessors (Freitas & Merkle's quote); the part to the left becomes a separate tape, and the machine stays in contact with the same constituent, now leftmost of its part. D on a leftmost constituent does nothing.
- Stoichiometry: the machine is a catalyst and the tape is replaced by all its pieces at halt (m + t -> m + pieces). A run that leaves the tape unchanged is elastic and not listed. No rates are published (Laing gave no dynamics), so rates are None.
- Forms: a molecule is wholly passive (a 0/1 tape) or wholly active (instruction constituents). The book's caption says the machine is 'built from an interpretation of another binary string', but no binary code is given, so machines are written with their instruction names instead of a code; there is no transition between the forms.
- The v1 params instruction_set (a fixed string) and string_length (had no default) are dropped: the instruction set is fixed by the decisions above, and tape lengths are set by the seed tapes and max_length. The v1 'compartments' capability and reactor are dropped: Laing's artificial organisms as compartments are not modelled.
- Universal computation ([485]) is shown with turing_program, a Chemart compiler from 2-symbol Turing machines to Laing machines (a jump to a state is 'CTq1 W1 CTq0', restoring the 0 at TTq0), not with Laing's own construction.
- Soup (Chemart addition): chemart.soup with constant dilution; a machine and a tape drawn in either order react; outflow is constant-total.

## Notes

Historically the first molecular machine chemistry; Laing never simulated any dynamics. The closure of the default incrementer on tape 0 is a binary counter: the 16 tapes of up to 4 bits, cut off when 1111 would grow into 00001.

---

*Specification: `catalog/chemistries/laing-molecular-machines.yaml` · generator: `chemart/chemistries/laing_molecular_machines.py` · tests: `tests/chemistries/test_laing_molecular_machines.py`*
