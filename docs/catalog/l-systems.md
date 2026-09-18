# L-systems

`l-systems` · *Lindenmayer, 1968; Prusinkiewicz & Lindenmayer, 1990*

*Also known as:* *Lindenmayer systems*, *D0L-systems*, *0L-systems*, *IL-systems*, *bracketed L-systems*, *stochastic L-systems*

Growth rather than chemistry. A word over an alphabet is rewritten by applying productions to *every* letter simultaneously - which is what distinguishes it from a grammar that rewrites one symbol at a time, and what makes it a model of cells dividing in parallel. Interpret the resulting letters as turtle commands and branching plants appear; the algae rule a -> ab, b -> a generates Fibonacci word lengths. There is no conservation and no kinetics, which is honest: it models development, not reactions.

| | |
|---|---|
| **family** | rewriting |
| **kind** | formalism |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.8 |
| **refs** | [512], [677], [376], [443], http://algorithmicbotany.org/papers/abop/abop.pdf |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): words over an alphabet whose symbols are one character plus an optional _subscript (a_r, F_b); brackets delimit branches and F f + - are turtle commands. reading words: species w<i> in discovery order, structure the word; reading symbols: species are the letters

**R — reactions** (explicit, arity 1): productions applied in PARALLEL to every letter of the word (unlike Chomsky grammars, which rewrite sequentially): a -> ab (0L), b < a -> b and F_a > F_b -> F_b (1L), 0 < 1 > 0 -> 1[+F1F1] (2L), F -> F[+F]F : 0.33 (stochastic). Letters without an applicable production are replaced by themselves. reading words: one reaction w -> w' per derivation step (per distinct successor word for a stochastic system); reading symbols: one reaction a -> multiset(chi) per production

**A — reactor**: maximally-parallel
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("l-systems", seed=1)
```

```
l-systems: 6 species, 5 reactions, status=complete
provides: initial-state, stoichiometry, topology
seed: 1
extras: analysis, axiom, ignore, probabilities, productions, reading, source, stationary, stochastic, system
```

First reactions:

```
w0 -> w1
w1 -> w2
w2 -> w3
w3 -> w4
w4 -> w5
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `system` | `enum` | `algae` | structural | published system from ABOP chapter 1: algae (eq. 1.2), anabaena (eq. 1.1), koch-snowflake (fig. 1.1 construction), koch-island (fig. 1.6), quadratic-koch-island and quadratic-snowflake (fig. 1.7), islands-and-lakes (fig. 1.8), koch-variant-a..f (fig. 1.9), dragon-curve and sierpinski-gasket (fig. 1.10), plant-a..f (fig. 1.24), stochastic-plant (sec. 1.7), signal-propagation (sec. 1.8), acropetal-signal and basipetal-signal (fig. 1.30), square-root-growth (eq. 1.5), hogeweg-hesper-a..e (fig. 1.31); custom uses axiom, productions, ignore and angle <br>one of `algae`, `anabaena`, `koch-snowflake`, `koch-island`, `quadratic-koch-island`, `quadratic-snowflake`, `islands-and-lakes`, `koch-variant-a` … and 24 more (32 total) |
| `reading` | `enum` | `words` | structural | words: species are the derivation words and each step is a reaction w -> w'; symbols: species are letters and each production a -> chi is a reaction (ABOP's growth matrix Q; context-free systems only) <br>one of `words`, `symbols` |
| `axiom` | `str` | `` | population | custom only: the nonempty starting word <br>*range:* e.g. 'F-F-F-F' or 'a_r' |
| `productions` | `list` | `` | structural | custom only: productions 'l < a > r -> successor : probability' (contexts, '*' and the probability optional; 'ε' or nothing for the empty successor); alternatives for the same predecessor and context need probabilities summing to 1 <br>*range:* e.g. ['a -> ab', 'b -> a'], ['b < a -> b', 'b -> a'], ['F -> F[+F]F : 0.5', 'F -> F[-F]F : 0.5'] |
| `ignore` | `str` | `` | structural | custom only: symbols skipped during context matching (#ignore) <br>*range:* ABOP 1.8 uses '+-' and '+-F' |
| `angle` | `float` | `` | spatial | custom only: turtle angle increment in degrees, recorded in extras.turtle (0 = none); named systems record their published angle <br>`0.0` … `360.0` |
| `iterations` | `int` | `5` | population | derivation length n: words up to n steps from the axiom (reading symbols: steps of extras.analysis.counts) <br>`0` … `200` · *range:* published derivation lengths: 2-10 for the curves, 4-7 for the plants, 24-30 for Hogeweg-Hesper (extras.turtle.published_iterations) |
| `max_species` | `int` | `2000` | structural | reading words: budget of words; a word with more combinations of stochastic choices than this is not expanded (status truncated) <br>`1` … `100000` · *range:* stochastic-plant: 301 words after 2 steps; 3 steps exceed any practical budget |
| `max_length` | `int` | `100000` | structural | reading words: successor words longer than this many symbols are not produced (status truncated) <br>`1` … `10000000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- parallel rewriting: the algae system a -> ab, b -> a gives a, ab, aba, abaab, abaababa with Fibonacci word lengths; from axiom b it gives b, a, ab, aba, abaab, abaababa, abaababaabaab (ABOP fig. 1.3)
- Anabaena catenula filament (ABOP eq. 1.1): a_r, a_l b_r, b_l a_r a_r, a_l a_l b_r a_l b_r, b_l a_r b_l a_r a_r b_l a_r a_r, cell counts following the Fibonacci numbers
- Koch constructions: the Koch curve has 4^n segments after n steps, the quadratic Koch island 4 * 8^n (fig. 1.6), the quadratic snowflake 5^n (fig. 1.7b)
- growth functions (ABOP 1.9): exponential for F -> FF (2^n), Fibonacci for eq. 1.2, polynomial with the Pascal triangle for a_i -> a_i a_{i+1}, floor(sqrt(n)) + 4 for the 2L-system of eq. 1.5
- context-sensitive signal propagation (ABOP 1.8): baaaaaaaa, abaaaaaaa, aabaaaaaa, ...; acropetal signals reach every segment above, basipetal signals only the path to the root
- stochastic L-systems give different specimens of the same plant species; the successor words of any word have probabilities summing to 1

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 9.8 (one paragraph: parallel versus sequential rewriting, recursion, turtle interpretation; no algorithm or example).
- Prusinkiewicz, P. & Lindenmayer, A. (1990). The Algorithmic Beauty of Plants. Springer (book [677]). Chapter 1: D0L definition and derivation (sec. 1.2, fig. 1.3, eq. 1.1 Anabaena), turtle interpretation and Koch curves (sec. 1.3, figs. 1.1, 1.6-1.9), edge rewriting (fig. 1.10), bracketed 0L (sec. 1.6.3, fig. 1.24), stochastic 0L (sec. 1.7), context-sensitive and bracketed context matching (sec. 1.8, figs. 1.29-1.31), growth functions (sec. 1.9, eqs. 1.2 and 1.5, Pascal table). http://algorithmicbotany.org/papers/abop/abop.pdf
- Prusinkiewicz, P. (1986). Graphical applications of L-systems. Proceedings of Graphics Interface '86, 247-253 (pL-systems: productions are tried in order, identity productions appended). http://algorithmicbotany.org/papers/graphical.gi86.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Book 9.8 only describes L-systems in prose, so the formalism and every example come from ABOP chapter 1. Productions printed inside figures (1.7-1.10, 1.24, 1.31) are lost by pdftotext and were read from rendered page images; figure 1.7a's successor is printed on two lines and is F+FF-FF-F-F+F+FF-F-F+F+FF+FF-F.
- Chemical reading: species are the words and each derivation step is a unimolecular reaction w -> w' (reading words, default); letters as species with one reaction per production is offered as reading symbols. The latter's stoichiometry is exactly ABOP's growth matrix Q (sec. 1.9), but it loses order and context, so context-sensitive systems raise an error there.
- Species ids are w0, w1, ... in derivation (breadth-first) order with the word as structure, because words grow exponentially and can be empty; words are deduplicated, so a derivation that cycles closes into a loop. The network is built with chemart.expand.expand (alternatives=True), cut at `iterations` steps from the axiom; status complete means every derivation of length <= iterations is present.
- A word that derives itself (a fixed point such as the end of the signal propagation) is an elastic step, dropped as in chemart.expand; its probability is kept in extras.stationary so that the alternatives of every expanded word still sum to 1.
- Stochastic 0L (ABOP 1.7): each letter occurrence chooses its production independently, so a word with k rewritable letters and 3 alternatives each has 3^k successor combinations. All distinct successor words become separate reactions; a word's probability is the sum over the combinations giving it. Probabilities are not rates (no kinetics is defined), so they are stored in extras.probabilities aligned with the reactions and the reactions carry rate None. max_species also caps the number of combinations enumerated for one word (status truncated), since stochastic-plant reaches 3^25 combinations at step 3.
- Nondeterministic 0L-systems without probabilities are not accepted: alternatives for the same predecessor and context must carry probabilities summing to 1 (ABOP's pi).
- Production notation is Chemart's text form of ABOP's: 'l < a > r -> chi : p' with spaces around < and >, '*' for no context (as in fig. 1.31) and the probability after ' : ' (ABOP prints it over the arrow). Subscripted letters are written with an underscore (a_r, F_b). The strict predecessor is one symbol; brackets replace themselves and cannot be rewritten; the axiom and successors must have balanced brackets.
- Context matching follows ABOP 1.8 and figure 1.29 (BC < S > G[H]M matches S in ABC[DE][SG[HI[JK]L]MNO]): left search skips ignored symbols, '[' and complete bracketed branches; right search skips ignored symbols and lateral branches, and ']' in the right context skips the rest of the branch. Among several matching context-sensitive productions for one letter, the first listed applies (GI'86 pL-system ordering); ABOP does not say.
- koch-snowflake (axiom F--F--F, F -> F+F--F+F, 60 degrees) is not printed as an L-system in ABOP chapter 1: it is the coding of the figure 1.1 construction (generator of 4 sides of length 1/3) that section 1.3 describes (initiator = axiom, generator = successor). Every other named system is transcribed from ABOP.
- algae: ABOP eq. (1.2) is a -> ab, b -> a from axiom a (the popular 'Lindenmayer algae' A -> AB, B -> A); figure 1.3 runs the same productions from axiom b. Lindenmayer (1968) [512] (J. Theor. Biol., paywalled) was not consulted; the examples and their numbers are verified against ABOP only.
- The eq. (1.5) growth function prints as 'fG(n) = sqrt n + 4' with the floor brackets lost; the derivation gives floor(sqrt(n)) + 4 symbols (X and the F's), which is what the test checks.
- Turtle interpretation is not executed: the published angle and derivation length of each named system are recorded in extras.turtle (angle, published_iterations); ABOP gives no step length d. The v1 'space' capability is removed because no positions are generated, and v1 'stochastic' is dropped (probabilities in the productions make a system stochastic).
- Parametric L-systems (ABOP 1.10), 3D turtle symbols and the growth of Anabaena with heterocysts (L-system 1.1) are not implemented; [376] (recurrence systems) and [443] (virtual cities) are cited by the book as context only.

## Notes

Growth/development rather than chemistry: no conservation and no kinetics, and a derivation step is unimolecular. The module also exports LSystem (successors, derive, derivation), parse_production, symbols and SYSTEMS for direct use.

---

*Specification: `catalog/chemistries/l-systems.yaml` · generator: `chemart/chemistries/l_systems.py` · tests: `tests/chemistries/test_l_systems.py`*
