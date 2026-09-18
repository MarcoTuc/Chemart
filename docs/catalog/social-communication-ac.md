# Artificial chemistry of social communication

`social-communication-ac` · *Dittrich, Kron & Banzhaf, 2003*

*Also known as:* *double contingency model*, *LuSi*, *luhmann3*

A social theory taken literally. Following Luhmann, the molecule is not a person but a *communication*: the population is the set of messages currently displayed, and a reaction is one agent observing another's message and choosing its own next one from its memory of what followed what. Agents are the context of a reaction, not its species. What emerges is social order - a shared, stable subset of the available messages - out of mutual expectation alone.

| | |
|---|---|
| **family** | non-chemical |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 20.4 |
| **refs** | [234], https://www.jasss.org/6/1/3.html, [235] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): the N possible activities (messages) 1..N, one species a_i per activity; a molecule is the message currently displayed on an agent's sign, so the population of displayed signs is a multiset of n_agents molecules. The species are communications, not persons: agents are the stateful context in which a reaction happens, and their memories are in extras.agents.

**R — reactions** (implicit, arity 2): a_i + a_k -> a_i + a_j : Ego observes Alter's displayed message i (a catalyst, Alter does not act), replaces its own displayed message k, and selects activity j with probability w_AP[j] computed from its two memory matrices; Ego, Alter and n_observers agents then memorise the event (i, j). The observed events of a run are returned with firing counts.

**A — reactor**: well-stirred-multiset
 · *dilution:* none: every agent always displays exactly one sign, so the population is constant at n_agents

## What you get

```python
net = chemart.generate_network("social-communication-ac", seed=1)
```

```
social-communication-ac: 10 species, 57 reactions, status=observed
provides: catalysts, initial-state, stoichiometry, topology
seed: 1
extras: agents, analysis, interaction_law
```

First reactions:

```
a10 + a3 -> a10 + a9  (x1)
a9 + a5 -> a9 + a6  (x1)
a10 + a6 -> a10 + a9  (x2)
a2 + a8 -> a2 + a8  (x1)
a4 + a9 -> a4 + a10  (x1)
a10 + a8 -> a10 + a9  (x1)
a9 + a1 -> a9 + a8  (x1)
a9 + a10 -> a9 + a6  (x2)
… and 49 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `N` | `int` | `10` | structural | number of possible activities (messages) the agents can show <br>`2` … `350` · *range:* paper: 2 (worked example), 10, 20, 40, 64 and up to 350 (fig. 3); software 2-100 |
| `n_agents` | `int` | `10` | population | number of agents M; also the number of displayed signs, i.e. the constant size of the molecule population <br>`2` … `100` · *range:* paper: 2 (dyadic) and 2-30 in figs. 7-9; software 2-100 |
| `interaction` | `enum` | `population` | structural | dyadic: the paper's two-agent world, Ego and Alter act alternately (Sec. 2.1); population: Ego and Alter are drawn at random each step (Sec. 6.1). dyadic requires n_agents = 2 <br>one of `dyadic`, `population` |
| `ee_memory` | `enum` | `alter` | structural | which memory the expectation-expectation reads: 'ego' is the basic model (Ego expects of itself what it did before, Sec. 2.5), 'alter' the variant of Sec. 6.6 (Ego expects of itself what others did when reacting to it). Only 'alter' can give scalable order <br>one of `ego`, `alter` |
| `n_observers` | `int` | `3` | structural | number of further agents that watch each interaction and store it in their alter-memory (Sec. 6.13); the information-proliferation mechanism that makes order scalable. Must be <= n_agents - 2 <br>`0` … `5` · *range:* paper: 0 (figs. 7-8) and 3 (fig. 9); software 0-5 |
| `alpha` | `float` | `` | selection | fraction of expectation-certainty in the activity value (eq. 2): 0 = expectation-expectation only, 1 = expectation-certainty only <br>`0.0` … `1.0` · *range:* paper: 0.0, 0.5, 1.0 |
| `gamma` | `float` | `2.0` | selection | exponent of the selection function (eq. 3): 0 is uniformly random choice, larger values make selection more deterministic <br>`0.0` … `40.0` · *range:* paper: 1 (proportional), 1.5, 2 (quadratic); gamma = infinity would be maximising selection and is not representable |
| `r_learn` | `float` | `0.5` | kinetic | learning rate: how much a memorised event raises its memory entry (eq. 8) <br>`0.0` … `2.0` · *range:* paper: 0.1 (worked example), 0.2 (figs. 2-6, 10, 11), 0.5 (figs. 7-9) |
| `r_forget` | `float` | `0.001` | kinetic | forgetting rate: added to every memory entry before the rows are renormalised (eq. 9), which pulls the memory back towards uniform <br>`0.0` … `1.0` · *range:* paper: 0.001 throughout, 0.01 in the worked example |
| `c_f` | `float` | `0.01` | selection | additive constant c_f of eq. 2; c_f/N keeps a small chance for every activity, so an agent that is as sure as possible still 'errs' with probability about c_f/(1+c_f) <br>`0.0` … `1.0` · *range:* paper's experiments: 0.01; the simulation software's default is 0.02, which is what the worked example of Sec. 2.13-2.27 used |
| `steps` | `int` | `1000` | stochastic | number of simulation steps; one step is one activity selection, i.e. one communication event <br>`0` … `100000` · *range:* paper: 1000 for the averaged figures, 4000 for the single runs behind figs. 10-11; software 50-10000 |
| `edge_threshold` | `float` | `0.1` | structural | edge threshold r of the activity graph (Sec. 7.4): only reactions of the average behaviour matrix with weight > r count when the activity systems are extracted <br>`0.0` … `1.0` · *range:* paper: 0.01 (fig. 10) and 0.1 (fig. 11); the software writes cutoffs 0.05, 0.10 and 0.15 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- dyadic order: with N = 20, gamma = 1, r_learn = 0.2 and alpha = 0 agents settle on about 15 of 20 activities per 50-step interval at average certainty O_AV about 0.58 (fig. 5, Sec. 5.14); gamma = 2 drives it down to about 2 activities
- order in the dyadic case does not need few messages: at gamma = 1, alpha = 0.5 the number of different activities saturates near 29-30 even for N = 300 (fig. 3, Sec. 5.11)
- order does not scale with the ego-memory: O_P falls from about 0.58 at M = 3 to 0.23 at M = 10 and 0.06 at M = 30 while single-agent certainty O_AV stays high (fig. 7), because every agent is sure but does something different (Sec. 7.10)
- with the alter-memory and alpha = 0, systems level order is maximal only in small populations: O_P is about 1 up to M = 6, 0.93 at M = 10, 0.50 at M = 20 and 0.21 at M = 30 (fig. 8)
- scalable social order: alter-memory, alpha = 0 and n_observers = 3 give O_P about 1 and about 3 different activities for every population size tested, M = 2 to 30 (fig. 9, Sec. 6.15)
- any expectation-certainty destroys scalability: with alpha = 0.5 or 1.0, O_P decays with M even with observers (fig. 9, Sec. 6.16)
- emergent activity systems: a small closed, self-maintaining set of activities that all agents share, e.g. {1, 5, 6} out of N = 10 with M = 10 (fig. 11); such a set is a chemical organisation in the sense of Fontana & Buss (Sec. 7.7)

## Sources

- Dittrich, P., Kron, T. & Banzhaf, W. (2003). On the scalability of social order: modeling the problem of double and multi contingency following Luhmann. Journal of Artificial Societies and Social Simulation 6(1):3. https://www.jasss.org/6/1/3.html (eqs. 1-3 and 8-15, the worked example of Sec. 2.13-2.27, figs. 3-11)
- Appendix to the above (Sec. 9.1 c_f, 9.4 memory models, 9.5 certainty measures, 9.6 the simulation software and its parameter file). https://www.jasss.org/6/1/3/appendix.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book (20.4) describes the model only in prose, so everything is reconstructed from [234]. The molecule is a communication, not a person: species are the N activities (the nodes of the book's fig. 20.8 activity graph), and an agent is the stateful context of a reaction. Since an activity is 'changing the number on his sign after observing the sign of the other agent', the event is the multiset reaction a_i + a_k -> a_i + a_j with Alter's displayed sign i as a catalyst (Alter does not act) and Ego's old sign k replaced by its selected activity j. This conserves the population (one sign per agent) exactly as the paper's algorithm does.
- Events in which Ego reselects the message it already displayed (j = k) are recorded as reactions with identical sides: they are genuine communication events, and keeping them makes the firing counts sum to `steps`.
- Eq. 2 adds c_f/N, but the printed worked example (Sec. 2.13-2.27) adds 0.01 at N = 2. The appendix's parameter file resolves it: the software's default is cf = 0.02, so c_f/N = 0.01 there, and eq. 2 is correct as printed. Sec. 2.13 lists every other setting of that example but not c_f. The default here is the 0.01 of the paper's experiments (all figure captions); the worked example is reproduced in the tests with c_f = 0.02.
- Eq. 10 as rendered in the HTML reads m_ij := m_ij + m_ij/sum_k m_ik, which is not a normalisation; the prose says 'we normalize every line' and m_ij := m_ij/sum_k m_ik reproduces the published memory matrices (0.545045, 0.454955) exactly, so the '+' is a typesetting error. Likewise eq. 12 prints 0.052803 for the activity value 0.502803.
- The figures were produced with memory model 05 (appendix 9.4), which is the simple neuronal matrix memory of Sec. 2.12 implemented here, and with the Shannon-entropy certainty measure (appendix 9.5). The other three memory models and three certainty measures of the appendix are not implemented.
- gamma = infinity (maximising selection) cannot be a JSON default; the software's own gamma range 0-40 is used instead. Sec. 6.13 does not say whether observers may be Ego or Alter themselves; they are drawn from the other agents, so an agent never stores its own activity as another's, and n_observers <= n_agents - 2 is enforced.
- Measurement follows the software's log files: the number of different activities is averaged over non-overlapping windows of intervalSize = 50 steps after startAverage = 500 (min(500, steps//2) for shorter runs), O_AV is the mean certainty of the acting agent's activity values over that phase (run.log column 3), and O_P (eq. 15) is evaluated on the final agent states, as the prediction game of Sec. 4 requires.
- An activity system (Sec. 7.6) is a set that is closed and self-maintaining in the activity graph thresholded at edge_threshold. The minimal such sets are exactly the bottom strongly connected components that feed themselves, which is what extras.analysis.activity_systems reports; enumerating all 2^N subsets is avoided.
- Reactions carry no rate: the paper gives selection probabilities that depend on agent state, not rate constants, so the selection rule is recorded in extras.interaction_law instead of being invented as a rate law.

## Notes

The book's fig. 20.8 (ten agents, ten activities) is this model's activity graph. It is the clearest case in the catalog of an artificial chemistry whose molecules are not chemical at all: Luhmann's elementary social unit is the communication, and social order appears as a chemical organisation over activities that reproduces itself independently of the agents carrying it.

---

*Specification: `catalog/chemistries/social-communication-ac.yaml` · generator: `chemart/chemistries/social_communication_ac.py` · tests: `tests/chemistries/test_social_communication_ac.py`*
