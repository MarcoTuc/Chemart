# Fraglets

`fraglets` · *Tschudin, 2003*

Network protocols written as a chemistry. A fraglet is a string of symbols that is simultaneously a data packet and a fragment of program: its head symbol is an instruction, its tail is payload. Instructions transform a fraglet, split it, fork it, or match two fraglets together, and fraglets are located in network nodes and can be sent between them. A confirmed-delivery protocol is then a handful of reactions, and protocol correctness becomes a question about a reaction network.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 16.2.1, 17.3.1; appendix (Fraglets.py, NetFraglets.py) |
| **refs** | [862], [575], [576], [577], [578], [579], [580], [581], [594], [595], [864], [914], [935], [936], [942], [943], [753], http://www.fraglets.net/ |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `compartments`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): Fraglets: strings of symbols [s1 s2 ... sn] located in a network node; the head symbols are the instruction, the tail is data, so a fraglet is at once a data packet and a program fragment. Species id = node[s1,s2,...,sn]; structure = node[s1 s2 ... sn]. A node is a well-stirred vessel; nodes attached to the same segment can exchange fraglets with send.

**R — reactions** (implicit, arity [1, 2]): Transformations: dup [dup t a tail] -> [t a a tail]; exch [exch t a b tail] -> [t b a tail]; nop [nop tail] -> [tail]; nul [nul tail] -> []; pop [pop t a tail] -> [t tail]; fork [fork a b tail] -> [a tail], [b tail]; split [split seq1 * seq2] -> [seq1], [seq2]; send n[send m tail] -> m[tail] (pycellchem, book table 16.1) or n[send seg m tail] -> m[tail] (fraglets-2007). fraglets-2007 adds copy [copy tail] -> 2 [tail]; empty [empty y n] -> [y], [empty y n tail] -> [n tail]; eq/lt [lt y n a b tail] -> [y a b tail] if a < b else [n a b tail]; length [length t tail] -> [t |tail| tail]; sum/sub/mult/div/mod/pow [sum t a b tail] -> [t a+b tail]; abs [abs t a tail] -> [t |a| tail]; pop2 [pop2 h t a b tail] -> [h a], [t b tail]. Reactions: match [match t tail1] + [t tail2] -> [tail1 tail2]; matchp [matchp t tail1] + [t tail2] -> [matchp t tail1] + [tail1 tail2] (catalytic).

**A — reactor**: compartments, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("fraglets", seed=1)
```

```
fraglets: 7 species, 4 reactions, status=complete
provides: catalysts, compartments, initial-state, stoichiometry, topology
seed: 1
extras: compartments, scheduling, segments
```

First reactions:

```
a[matchp,cdp,send,b,split,send,a,ack,*] + a[cdp,data] -> a[matchp,cdp,send,b,split,send,a,ack,*] + a[send,b,split,send,a,ack,*,data]
a[send,b,split,send,a,ack,*,data] -> b[split,send,a,ack,*,data]
b[split,send,a,ack,*,data] -> b[send,a,ack] + b[data]
b[send,a,ack] -> a[ack]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `program` | `str` | `a a net a b net f a[matchp cdp send b split sen…` | population | the Fraglets program: initial fraglets per node (with multiplicities) and the network topology <br>*range:* any .fra text: 'f node[symbols]mult' fraglet lines, 'a node segment' attachment lines, '#' comments, 'e' end; symbols separated by spaces or ':'. Default: the active confirmed delivery protocol of AINS 2003 sec. IV.B / book fig. 16.10 |
| `dialect` | `enum` | `pycellchem` | structural | instruction set and edge cases: pycellchem = book table 16.1 and the book's reference interpreter (send dest tail); fraglets-2007 = upstream C interpreter fraglets0.32 and its 2007 instruction set (send seg dest tail, logic and arithmetic) <br>one of `pycellchem`, `fraglets-2007` |
| `method` | `enum` | `closure` | structural | closure: every reaction reachable from the program fraglets (complete/truncated); ssa: a PyCellChemistry-scheduled run from the program multiset, observed reactions with counts <br>one of `closure`, `ssa` |
| `max_species` | `int` | `200` | structural | closure only: species budget; status truncated when exceeded (elongating programs never close) <br>≥ `1` |
| `steps` | `int` | `1000` | population | ssa only: number of match reactions to schedule (the run stops earlier when the nodes are inert) <br>`0` … `10000000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- the active confirmed delivery protocol delivers [data] to node b and returns [ack] to node a, carrying its own code to b (AINS 2003 sec. IV.B, book fig. 16.10)
- the passive CDP with pre-installed code at b gives the same result (AINS 2003 sec. IV.B)
- matchp acts as a catalyst: persistent program storage
- a quine [fork nop b match b fork nop b] regenerates itself through fork, nop and match (PyCellChemistry quine)
- tail buildup: [split matchp a dup a * a a] makes ever longer [a a ... a] fraglets (PyCellChemistry codegrowth), so the closure never finishes
- lossy link emulation: two competing matchp rules drop a fraglet in proportion to their multiplicities (tutorial 2007: 44 and 74 of 100 delivered for 50% and 25% loss)
- Flow control with credits, self-healing multipath routing, resilient load balancing (the disperser), autocatalytic recovery from code deletion, and the evolvability problems (code bloat, tail buildup, broken pathways) are reported in [862, 575-581, 943] and are not reproduced here

## Sources

- Tschudin, C. (2003). Fraglets - a metabolistic execution model for communication protocols. Proc. 2nd Annual Symposium on Autonomous Intelligent Networks and Systems (AINS): instruction set (sec. III.B), header rewriting (fig. 3), passive and active confirmed delivery protocol with execution trace (sec. IV.B, fig. 4), flow control with credits (appendix). Shipped as doc/2003-ains.pdf in fraglets0.32.tgz, http://www.fraglets.net/pkg/fraglets0.32.tgz (Wayback Machine, 2011-07-23)
- Tschudin, C., Yamamoto, L. et al. Fraglets Instruction Set, 2007-09-24. http://www.fraglets.net/frag-instrset-20070924.txt (Wayback Machine, 2007-10-25)
- fraglets0.32 interpreter, C source (GPL), src/frag.c: f_nop, f_copy, f_fork, f_dup, f_empty, f_exch, f_lt, f_eq, f_length, f_op2, f_abs, f_nul, f_pop, f_pop2, f_split, f_send, rx_fraglet, f_execFromQueue, f_parsefrag, name_add; examples fra/examples/cdp.fra, acdp.fra, rescdp.fra
- Yamamoto, L. (2007). Fraglets Tutorial, BIONETS SP2, Brussels, 29 Jan 2007: header rewriting, append symbol, code mobility, lossy link emulation (b[msg]44 and 74), increment, prepend, append, del, count, getmin. http://www.fraglets.net/pres/tut2007-01-29.pdf (Wayback Machine, 2011-07-23)
- Yamamoto, L. PyCellChemistry (GPL-3), src/Fraglets.py (instruction set, r_* rules, inject, propensity, react, quine and codegrowth examples), src/NetFraglets.py (CDP), src/artchem/Cell.py (hierarchical Gillespie), src/artchem/KeyMultiset.py and Multiset.py (random picks weighted by multiplicity). https://github.com/laryamamoto/PyCellChemistry

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Two dialects instead of one: the book's table 16.1 and its reference code use 'n[send m tail] -> m[tail]' and the ten PyCellChemistry opcodes; the upstream interpreter and the 2007 instruction set use 'send seg dest tail' and add logic and arithmetic. The dialects also differ on short fraglets: PyCellChemistry turns [dup a] into [a], [exch a b] into [a b], [pop a] or [pop a b] into [a], and [fork a] into [a], and never stores a one-symbol instruction fraglet; frag.c drops all of these (and [pop a b] gives [a]). Each dialect follows its source exactly.
- Symbols: PyCellChemistry encodes each symbol as one character (its parser keeps only the first letter of a symbol, marked TMP in the source); whole symbols are used here, as in the book and frag.c. The PyCellChemistry CDP example 'matchp c send b split send a k *' is the book's 'matchp cdp send b split send a ack *' with abbreviated symbols.
- Topology: PyCellChemistry connects nodes explicitly (add_cnx); here both dialects read the upstream 'a node segment' lines, and in the pycellchem dialect a node can send to every other node attached to one of its segments (NetFraglets' CDP connects a->b and b->a, not a->a).
- fraglets-2007 numbers: a symbol starting with a digit (or '-' and a digit) is an integer, as in frag.c name_add (atol, so '007' and '7x' both read 7). Arithmetic uses Python integers (frag.c: C int, overflow not modelled); div and mod truncate toward zero like C; div/mod by 0, a negative exponent and 0^0 drop the fraglet. Matching numbers is by exact symbol: frag.c's number matching path is unfinished ('stopped here!').
- fraglets-2007 instructions that need time (delay, wait), node creation (newnode, inject, expel), printing (printsym; send stdout/stderr simply removes the fraglet), broadcast/anycast, splitat, newname or the host application ('_' trigger, out) are not ported; a program using them is rejected. The upstream examples cdp.fra and acdp.fra depend on the '_' trigger of sys-send-deliver and are therefore not reproduced; the AINS 2003 versions of both CDPs are.
- Scheduling: the upstream C reactor (maximally parallel 'min' algorithm) is not ported; method ssa uses the PyCellChemistry scheduler for both dialects. Transformations are instantaneous, so no finite rate can be written for them and all reaction rates are None; extras.scheduling states the propensity n_active * n_passive (k = 1) of matches. The mass-action rate constants of Meyer et al. [575, 580] are per-program choices not given in the book, so the v1 rate_constants parameter is dropped, as are the v1 matrix-typed topology/seed_code (now the program text) and scheduler (maximal-parallel and next-reaction are not implemented).
- ssa details: PyCellChemistry computes the propensities node by node, so a fraglet sent to an already-visited node waits one step; here all transformations in all nodes are settled before every draw and once more after the last one. Injection of new input during the run (NetFraglets adds a data fraglet every 4 iterations) is not modelled: all input is in the program. A transformation that reproduces itself (e.g. [dup dup dup]) would loop forever in both interpreters; it is left inert and listed in extras.elastic. extras.final_state, bimolecular_events, inert and stopped describe the run.
- Closure: chemart.expand.expand with arity 1 and 2 over ordered pairs, from the distinct program fraglets; a match needs the active and passive fraglets in the same node. Closure ignores multiplicities, so data-dependent recursions (counters) can show reactions no single run takes.

## Notes

The chemical disperser of section 17.3.1 is catalogued separately as disperser. Fraglet programs are write-only: the tutorial derives them bottom-up from the desired output, and the default program is kept small (7 species, 4 reactions) so the closure reads like the AINS trace.

---

*Specification: `catalog/chemistries/fraglets.yaml` · generator: `chemart/chemistries/fraglets.py` · tests: `tests/chemistries/test_fraglets.py`*
