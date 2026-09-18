# Reflexive artificial chemistry

`reflexive-ac` · *Salzberg, 2007*

*Also known as:* *graph-based reflexive artificial chemistry*, *graph-constructing graphs*, *reflexive composition of state machines*

Machines that compose with machines, including themselves. A molecule is a finite state machine drawn as a graph with a current state; a reaction feeds one machine's output into the other and keeps the reachable graph of paired states as the product. Composition multiplies state spaces, so products grow quickly and are self-similar. The reflexive case is the interesting one: a machine composed with itself is a machine reading its own structure, which is the sense in which it can act on itself.

| | |
|---|---|
| **family** | rewriting |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.8 |
| **refs** | [737], doi:10.1016/j.biosystems.2005.12.008, doi:10.7551/mitpress/1429.003.0084, arXiv:2505.07186 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): finite state machines (Mealy machines) drawn as labelled directed graphs: every node has at most one link per input symbol, labelled input/output, and one node is the current state (the pointer). Input data are graphs of the same kind (missing links allowed). Species id: canonical text with nodes numbered breadth-first from the pointer, node j listing destination/output per input, e.g. the adding machine M45 at S0 is 0/0,1/0;1/1,0/1

**R — reactions** (implicit, arity 2): s1 + s2 -> s1 + s2 + s1 o s2: the reachable composite state graph of s1 (sender) and s2 (receiver) under reflexive composition; elastic if it has more than max_states nodes

**A — reactor**: graph-rewrite, well-stirred-multiset
 · *dilution:* closure: none; soup: each reactive collision adds the product and removes a random molecule (population M constant)

## What you get

```python
net = chemart.generate_network("reflexive-ac", seed=1)
```

```
reflexive-ac: 6 species, 7 reactions, status=complete
provides: catalysts, stoichiometry, topology
seed: 1
extras: analysis, notation, seed
```

First reactions:

```
2 0/0,1/0;1/1,0/1 -> 2 0/0,1/0;1/1,0/1 + 0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1
0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1 + 0/0,1/0;1/1,0/1 -> 0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1 + 0/0,1/0;1/1,0/1 + 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0
2 0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1 -> 2 0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1 + 0/0,1/0;2/1,3/1;4/0,5/0;6/1,7/1;8/1,9/1;10/0,11/0;12/1,13/1;14/0,15/0;13/0,12/0;15/1,14/1;9/0,8/0;11/1,10/1;5/1,4/1;7/0,6/0;1/1,0/1;3/0,2/0
1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0 + 0/0,1/0;1/1,0/1 -> 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0 + 0/0,1/0;1/1,0/1 + 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;16/1,17/1;18/0,19/0;20/1,21/1;22/0,23/0;24/1,25/1;26/0,27/0;28/1,29/1;30/0,31/0;29/0,28/0;31/1,30/1;25/0,24/0;27/1,26/1;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0;21/0,20/0;23/1,22/1;17/0,16/0;19/1,18/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0
1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0 + 0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1 -> 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0 + 0/0,1/0;2/1,3/1;3/0,2/0;1/1,0/1 + 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;16/1,17/1;18/0,19/0;20/1,21/1;22/0,23/0;24/1,25/1;26/0,27/0;28/1,29/1;30/0,31/0;32/0,33/0;34/1,35/1;36/0,37/0;38/1,39/1;40/1,41/1;42/0,43/0;44/1,45/1;46/0,47/0;48/0,49/0;50/1,51/1;52/0,53/0;54/1,55/1;56/1,57/1;58/0,59/0;60/1,61/1;62/0,63/0;55/0,54/0;53/1,52/1;51/0,50/0;49/1,48/1;63/1,62/1;61/0,60/0;59/1,58/1;57/0,56/0;15/0,14/0;13/1,12/1;5/0,4/0;3/1,0/1;23/1,22/1;21/0,20/0;19/1,18/1;17/0,16/0;39/0,38/0;37/1,36/1;35/0,34/0;33/1,32/1;47/1,46/1;45/0,44/0;43/1,42/1;41/0,40/0;11/0,10/0;9/1,8/1;7/0,6/0;2/1,1/1;31/1,30/1;29/0,28/0;27/1,26/1;25/0,24/0
0/0,1/0;2/1,3/1;4/0,5/0;6/1,7/1;8/1,9/1;10/0,11/0;12/1,13/1;14/0,15/0;13/0,12/0;15/1,14/1;9/0,8/0;11/1,10/1;5/1,4/1;7/0,6/0;1/1,0/1;3/0,2/0 + 0/0,1/0;1/1,0/1 -> 0/0,1/0;2/1,3/1;4/0,5/0;6/1,7/1;8/1,9/1;10/0,11/0;12/1,13/1;14/0,15/0;13/0,12/0;15/1,14/1;9/0,8/0;11/1,10/1;5/1,4/1;7/0,6/0;1/1,0/1;3/0,2/0 + 0/0,1/0;1/1,0/1 + 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;16/1,17/1;18/0,19/0;20/1,21/1;22/0,23/0;24/1,25/1;26/0,27/0;28/1,29/1;30/0,31/0;32/0,33/0;34/1,35/1;36/0,37/0;38/1,39/1;40/1,41/1;42/0,43/0;44/1,45/1;46/0,47/0;48/0,49/0;50/1,51/1;52/0,53/0;54/1,55/1;56/1,57/1;58/0,59/0;60/1,61/1;62/0,63/0;55/0,54/0;53/1,52/1;51/0,50/0;49/1,48/1;63/1,62/1;61/0,60/0;59/1,58/1;57/0,56/0;15/0,14/0;13/1,12/1;5/0,4/0;3/1,0/1;23/1,22/1;21/0,20/0;19/1,18/1;17/0,16/0;39/0,38/0;37/1,36/1;35/0,34/0;33/1,32/1;47/1,46/1;45/0,44/0;43/1,42/1;41/0,40/0;11/0,10/0;9/1,8/1;7/0,6/0;2/1,1/1;31/1,30/1;29/0,28/0;27/1,26/1;25/0,24/0
1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;16/1,17/1;18/0,19/0;20/1,21/1;22/0,23/0;24/1,25/1;26/0,27/0;28/1,29/1;30/0,31/0;29/0,28/0;31/1,30/1;25/0,24/0;27/1,26/1;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0;21/0,20/0;23/1,22/1;17/0,16/0;19/1,18/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0 + 0/0,1/0;1/1,0/1 -> 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;16/1,17/1;18/0,19/0;20/1,21/1;22/0,23/0;24/1,25/1;26/0,27/0;28/1,29/1;30/0,31/0;29/0,28/0;31/1,30/1;25/0,24/0;27/1,26/1;9/1,8/1;11/0,10/0;2/1,1/1;7/0,6/0;21/0,20/0;23/1,22/1;17/0,16/0;19/1,18/1;13/1,12/1;15/0,14/0;3/1,0/1;5/0,4/0 + 0/0,1/0;1/1,0/1 + 1/0,2/0;0/0,3/0;4/1,5/1;6/1,7/1;8/0,9/0;10/1,11/1;12/0,13/0;14/1,15/1;16/1,17/1;18/0,19/0;20/1,21/1;22/0,23/0;24/1,25/1;26/0,27/0;28/1,29/1;30/0,31/0;32/0,33/0;34/1,35/1;36/0,37/0;38/1,39/1;40/1,41/1;42/0,43/0;44/1,45/1;46/0,47/0;48/0,49/0;50/1,51/1;52/0,53/0;54/1,55/1;56/1,57/1;58/0,59/0;60/1,61/1;62/0,63/0;55/0,54/0;53/1,52/1;51/0,50/0;49/1,48/1;63/1,62/1;61/0,60/0;59/1,58/1;57/0,56/0;15/0,14/0;13/1,12/1;5/0,4/0;3/1,0/1;23/1,22/1;21/0,20/0;19/1,18/1;17/0,16/0;39/0,38/0;37/1,36/1;35/0,34/0;33/1,32/1;47/1,46/1;45/0,44/0;43/1,42/1;41/0,40/0;11/0,10/0;9/1,8/1;7/0,6/0;2/1,1/1;31/1,30/1;29/0,28/0;27/1,26/1;25/0,24/0
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `closure` | structural | closure: every composition reachable from the seed machines (chemart.expand over ordered pairs), complete or cut off by max_species; soup: a stochastic flow reactor of M molecules, observed reactions with firing counts <br>one of `closure`, `soup` |
| `machines` | `list` | `['M45']` | structural | seed machines: elementary numbers M<n> (optionally @S0/@S1 for the pointer) or machine texts 'dest/out,...;...' with '-' for a missing link and optional '@j' pointer; all must share the alphabet. Empty: M random machines <br>*range:* e.g. ['M45', 'M61'], ['M54'], ['0/0,1/1;1/1,0/1@1']; the arXiv:2505.07186 machines of interest are M7, M44, M45, M54, M60, M61 |
| `M` | `int` | `50` | population | soup: reactor size (given machines are split equally); without machines, the number of random machines drawn (closure seeds or initial soup population) <br>`2` … `100000` |
| `states` | `int` | `2` | population | random machines: number of states before trimming to the part reachable from state 0; links are uniform over states and outputs <br>`1` … `1000` |
| `symbols` | `int` | `2` | structural | random machines: input/output alphabet size k <br>`1` … `16` · *range:* the sources use the binary alphabet {0, 1} |
| `max_states` | `int` | `64` | structural | largest product accepted: a composition with more reachable nodes is elastic <br>`1` … `100000` |
| `max_species` | `int` | `100` | structural | closure only: species budget; reactions that would add more are dropped and status becomes truncated <br>`1` … `100000` |
| `collisions` | `int` | `2000` | population | soup only: number of collisions (elastic ones included) <br>`0` … `10000000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- composition multiplies state spaces: repeated composition yields exponentially larger, self-similar state-space graphs (M45 gives 2, 4, 16, 32, 64 states)
- a graph composed with itself reads its own structure (reflexivity): s + s -> 2 s + s o s
- the n-fold composite of a machine is the lattice of reflexive composition (arXiv:2505.07186 eqs. 3-4)
- M45 (adding machine modulo 2) reproduces elementary CA rule 90 and the Sierpinski triangle on even steps (figs. 5-7)
- M54 and M60 run rule 90 backwards: each even step is a rule-90 preimage of the previous one (figs. 10-15); M45 and M54 are transposes of each other (fig. 17)
- M44 shows reversible billiard-ball dynamics, inverted by flipping the lattice (fig. 8)
- of the 256 two-state two-symbol machines, 76 are unique under complement and mirror (Table 2)
- Salzberg 2007 (abstract): continuous emergence of complex self-similar topologies, novel reaction pathways, seemingly open-ended diversity

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 9.8 (one sentence: finite state machines and their input data deconstructed into a graph that rewrites a graph, and in principle itself).
- Salzberg, C. (2007). A graph-based reflexive artificial chemistry. BioSystems 87(1):1-12 (book [737]); abstract only, the full text was not accessible: FSM machine graph and input string dissolved into one space of directed graphs, graphs interact with their own topology to generate a product, self-similar topologies, novel reaction pathways and open-ended diversity. https://pubmed.ncbi.nlm.nih.gov/16733079/
- Salzberg, C. & Sayama, H. (2025). Reflexive composition of elementary state machines, with an application to the reversal of cellular automata rule 90. arXiv:2505.07186: the formulation of the 'composite state machine' introduced in Salzberg 2006a, 2006b and 2007 (Mealy machine G = (Q, A, delta, phi); reflexive composition eqs. 1-4 with sender/receiver alternation; composition generating exponentially larger state-space graphs), the naming scheme of Table 1, Table 2, figs. 2-17. https://arxiv.org/abs/2505.07186
- Salzberg, C., Sayama, H. & Ikegami, T. (2004). A tangled hierarchy of graph-constructing graphs. Artificial Life IX, pp. 495-500, MIT Press (open access): Mealy machines with a shared input/output alphabet, tapes as chains of links, pointers to the current states, symmetric engagement of machine and tape graphs. https://doi.org/10.7551/mitpress/1429.003.0084

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The primary source [737] and its companions (Salzberg 2006a, ALife X; 2006b, Artificial Life 12:487-512) were not accessible (paywalled; no preprint, thesis or archived copy found). The reaction is reconstructed from the 2025 restatement of the formulation by the same author, which says the earlier studies composed distinct sender and receiver machines repeatedly into exponentially larger state-space graphs. The exact product construction, elastic conditions and reactor of the 2007 paper are therefore Chemart's reading, not a transcription.
- Molecule: a deterministic Mealy machine over the alphabet {0..k-1} with a pointer (the current state, the dark node of the 2025 figures), trimmed to the nodes reachable from the pointer. Missing links are allowed so that input data (e.g. the tape chains of the 2004 paper, with a symbol playing xi) are graphs of the same kind; a composite link exists only if both links it passes through exist.
- Product s1 o s2: the composite of eqs. 1-2 with s1 as sender, whose node (x, y) holds the sender-position and receiver-position states; after each input the machines swap positions (the 2025 'folding'). For distinct species the node also records which machine sends, so the product has at most 2|s1||s2| nodes; for s + s the machines are identical and it has at most |s|^2 (the uniform lattice of the 2025 paper). Its pointer is the pair of reactant pointers and only reachable nodes are kept.
- Reaction scheme: both reactants are catalysts, s1 + s2 -> s1 + s2 + s1 o s2, as in the book's constructive chemistries (AlChemy, eq. 9.4); the 2007 abstract says only that interacting graphs 'generate a product'. Ordered pairs: s1 o s2 and s2 o s1 differ.
- Not modelled: the 2004 graph-constructing-graph machinery (xi-flows chosen by a decision tree over N-loops, translation of the output string by codons into construction-head instructions that grow and fold the graph). The 2025 restatement does not use it, and the 2004 paper leaves conflict resolution open.
- Elastic boundary: products with more than max_states reachable nodes (Chemart choice, since state counts multiply at every composition). No rates are published, so reactions carry rate None.
- Canonical id: breadth-first numbering from the pointer with links taken in input order; for accessible deterministic graphs this is a complete isomorphism invariant, so no general graph canonisation is needed.
- Soup: chemart.soup with constant dilution draws an ordered pair (sender, receiver), adds the product and removes a random molecule; initial population: the given machines in equal copies, or M random machines (duplicates allowed). Random machines have uniform links over states and outputs, trimmed from state 0. These are Chemart choices; the defaults (closure from M45, M = 50, 2000 collisions, max_states 64, max_species 100) are small, not published values.
- The v1 entry had no parameters; all parameters are new. v1 said 'dilution: none', which holds for the closure; the soup is an added reactor.
- Published numbers checked by computation: Table 1 encoding (M61 = 00111101, M45 = 00101101), 76 machines unique under complement and mirror, the Table 2 equivalents, transposition M45 <-> M54 and M60 -> M60 (fig. 17, footnote 3), M45 equal to rule 90 on even steps with the fig. 7 boundary inputs, M54 even steps as rule-90 preimages with the eqs. 6-7 inputs, M44 inverted by flipping the lattice. The fig. 12 caption's input sequence (0 up to t = 5, then 1, 1) was not reproduced: with the constraint algorithm as described, the run from two centred 1 cells needs input 1 at t = 0; the figure's step indexing cannot be recovered from the text. Likewise fig. 13's statement that M60 with the same initial states and boundary inputs has the same even rows as M54 did not hold when M60 was fed M54's inputs, so only M54's reversal is tested.

## Notes

Chemart observation, not a published result: random soups of two-state machines usually collapse to one species within a few thousand collisions, typically a small machine that reproduces itself under self-composition (s o s = s). The module exports canonical, parse, compose, elementary and number for working with the machines directly.

---

*Specification: `catalog/chemistries/reflexive-ac.yaml` · generator: `chemart/chemistries/reflexive_ac.py` · tests: `tests/chemistries/test_reflexive_ac.py`*
