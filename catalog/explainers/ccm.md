## Introduction

The Chemical Casting Model is Yasusi Kanada's attempt, in the early 1990s, to
solve hard combinatorial problems the way a chemical system reaches order: by
many small, random, local events, with no part of the system ever looking at
the whole. Technically it is a *production system*, a set of if-then rewriting
rules like those of the expert systems of the time, but Kanada wanted rules
much more primitive than expert-system rules, so that behaviour emerges from
applying them over and over instead of being programmed into them. Each rule is
then less like a line of a program and more like the formula of a single
chemical reaction, which is why he called it a *chemical* casting model.

The running example is the **eight-queens puzzle**. Place eight chess queens
on an 8×8 board so that no two attack each other, which means no two may share
a row, a column or a diagonal. There are 92 ways to do it. A conventional
program searches the board systematically, or scores whole boards and improves
the score. CCM does neither. It treats the board as a test tube in which each
queen is a molecule, and lets a few queens at a time "react" by trading places,
but only when the trade does not make things worse *among the queens that took
part*. Nothing in the system knows how good the whole board is. Yet the board
reliably ends up solved.

The same machinery colours maps. Each region is a molecule, a reaction
recolours one region while its neighbours look on, and the run ends with a
valid colouring, for example of the 48 mainland US states in four colours.

Kanada used CCM to argue for *emergent computation*: a global solution arising
from local rules plus randomness, rather than from an algorithm that holds the
whole problem in view. Banzhaf and Yamamoto describe it in their chapter on
applications, under "Search and Optimization Algorithms Inspired by Chemistry"
(book §17.2.2), next to the molecular travelling salesman.

## How it works

### The working memory is one candidate solution

Most artificial chemistries hold a large population of molecules, and the
interesting thing is which kinds become common. CCM is different: its "test
tube" (Kanada calls it the *working memory*) holds exactly one atom per queen,
or one per map region. Together they make one complete candidate answer. A
reaction does not create or destroy atoms; it changes their *states*.

For N queens there is one queen per row, and a queen's state is the column it
sits in. The species `q3=4` means "the queen of row 3 is in column 4", with
rows and columns numbered from 0. A board of N queens therefore uses N of the
N² possible species at any moment. For map colouring, `v7=c2` means "vertex 7
has colour 2", and the map's borders are fixed links between vertices.

### Local order degree: how a reaction is judged

Every pair of atoms has a **local order degree** (LOD), a score that says
whether that pair is in a good relation:

- two queens score **1** if they cannot take each other, and **0** if they
  share a column or a diagonal;
- two vertices score **1** if they are not neighbours or have different
  colours, and **0** if they are neighbours of the same colour.

A reaction matches a few atoms. Adding up the LODs of the pairs among them
gives the *instance order degree*. The reaction fires only if that sum is not
lower afterwards than before. That is the whole acceptance test. The sum over
*all* pairs, the **global order degree**, is never computed by the dynamics.
For eight queens it runs from 0 (every pair attacks) to 28 (all 28 pairs safe,
a solution). Chemart reports it only so you can watch it.

### The rules, and what catalysts do

The n-queens rule swaps the columns of two queens. Because each queen keeps its
row and the two trade columns, there is always exactly one queen per row and
per column. Only diagonals can go wrong.

A swap can also involve **catalysts**: extra queens that do not move, but whose
pairs with the two swapped queens count in the test. Here is the first reaction
of the default run, with one catalyst:

```
q3=3 + q4=4 + q6=6 -> q3=4 + q4=3 + q6=6
```

The queens of rows 3 and 4 swap columns; the queen of row 6 is the catalyst.
Before the swap all three sit on the main diagonal, so every pair attacks and
the sum is 0. After it, q3 and q4 are still diagonal to each other (0), but
q3 at (3, 4) and q6 at (6, 6) are safe (1), and so are q4 at (4, 3) and q6
(1). The sum rises from 0 to 2, so the swap is accepted.

The number of catalysts sets how *local* the rule is. With none, no pair is
evaluated at all, so every swap is accepted and the board wanders at random for
ever. With more, each reaction looks at more of the board, fewer reactions are
needed, but finding a reaction that passes costs more attempts and the system
can more easily get stuck.

The *variable-catalyst* rule for queens is a different kind of move: it shifts
a single queen to a new column, with all the other queens as catalysts. That can
put two queens in one column, so here the column also counts in the LOD.

For map colouring, the rule recolours one vertex, with one, two or three of its
neighbours as catalysts, or, in the variable-catalyst rule, all of them.

### Getting unstuck: frustration

A purely non-decreasing test can trap the system in a state that no single
reaction improves. Kanada's fix, the **frustration accumulation method** (FAM),
works like a local, self-tuning temperature. Every atom carries a small
frustration value, initially `f0`. Each time a reaction that would have changed
the atom is tried and refused while some of the matched constraints are
violated, its frustration is multiplied by `c`. The frustration of the atoms a
reaction would change is subtracted from the "before" score in the test, so
a frustrated atom will eventually accept a move that makes things slightly
worse. When the atom finally reacts, its frustration drops back to `f0`. There
is no global temperature and no cooling schedule: annealing happens only where
the system is stuck.

### When it stops

Once every constraint among the matched atoms is satisfied, nothing changes any
more, so a solved board is stationary. Following Kanada's demo programs, a run
ends when a long streak of attempted reactions all fail, and `max_tests` is a
hard limit on the number of attempts.

## Using it

The default run is the classic experiment: eight queens on the diagonal, the
single-catalyst swap rule, and frustration on. It is solved after 188 accepted
reactions out of 1,637 attempts. The run's history is in `net.extras`:

```python
a = net.extras["analysis"]
a["solved"], a["reactions"], a["tests"]      # (True, 188, 1637)
net.extras["initial_assignment"]             # [0, 1, 2, 3, 4, 5, 6, 7]   the diagonal
net.extras["final_assignment"]               # [4, 6, 1, 5, 2, 0, 3, 7]   a solution
a["god_initial"], a["god_final"]             # (0, 28)
```

`final_assignment[i]` is the column of the queen in row `i`. The network itself
records each distinct accepted reaction once, with how often it happened, so a
count like 187 reactions in the summary is the number of *different* swaps,
while `a["reactions"]` counts every accepted one.

To watch the run rather than its outcome, `chemart.evolve` returns a
trajectory. It has a frame for the starting board, one after every accepted
reaction, and a last one at the final attempt when the run ended on a streak of
failures. Time is counted in attempted reactions (tests). Each frame holds the
board (`state`, one species per queen), the reaction just accepted (`fired`)
and the global order degree as the observable `god`, the curve Kanada plots:

```python
traj = chemart.evolve("ccm", seed=1)
god = traj.series("god")
len(traj.frames), traj.times()[-2:]          # (190, [573.0, 1637.0])
god[:10]                                     # [0, 12, 16, 21, 25, 25, 23, 24, 25, 26]
god[-3:]                                     # [24, 28, 28]
traj.frames[1].fired                         # [[['q3=3', 'q4=4', 'q6=6'], ['q3=4', 'q4=3', 'q6=6'], 1]]
```

The board was solved at the 573rd attempt. The run then went on failing,
because a solved board never changes, until the termination test stopped it at
attempt 1,637.

**Colouring the USA map.** Set `problem="graph-coloring"`. The default graph is
the 48 contiguous states with their 106 borders, starting all one colour:

```python
net = chemart.generate_network("ccm", seed=1, problem="graph-coloring",
                               rule="variable-catalyst")
net.extras["analysis"]["solved"], net.extras["analysis"]["reactions"]   # (True, 118)
```

**Seeing what catalysts do.** Compare rules on the same seeds. Over ten seeds
without frustration, eight queens needs on average about 146 reactions with one
catalyst, 130 with two and 53 with three, while the number of attempts falls
much less (about 1,570 to 1,200):

```python
import statistics
for rule in ("single-catalyst", "double-catalyst", "triple-catalyst"):
    runs = [chemart.generate_network("ccm", seed=s, rule=rule, frustration=False)
            .extras["analysis"] for s in range(10)]
    print(rule, statistics.mean(r["reactions"] for r in runs))
```

`rule="no-catalyst"` never stops on its own. It always runs to `max_tests`.

**Seeing what frustration does.** The variable-catalyst colouring rule without
frustration sometimes halts on a wrong colouring. Over seeds 0–9 with
`max_tests=30000`, it solves the map 7 times out of 10 with `frustration=False`
and 10 out of 10 with the default FAM.

**Other problems.** `N` sets the board size, `initial="random"` starts from a
random permutation instead of the diagonal, and `acceptance="increasing"`
switches to the book's stricter test (the score must strictly rise). For
colouring, `graph="random"` draws a graph with `V` vertices and edge
probability `edge_probability`, the model behind the DSJC benchmark graphs, and
`graph="custom"` takes your own `edges`.

## Results

**Eight queens, solved every time.** Kanada and Hirokawa's HICSS-27 paper
(1994) introduces the model with N queens. Starting from random layouts, their
runs "never fail to be solved", and they measured the average number of
matches and reactions up to N = 50. Their Table 1 puts the time for N queens at
roughly O(N^4.6). The same table lists other problems cast this way: travelling
salesperson (optimal in 97 of 100 trials for 10 cities), 0–1 knapsack (optimal
in 45 of 100 for 20 items), graph colouring and sorting, each with a single
rule and a single LOD. Chemart reproduces the n-queens claim: its tests solve
eight queens from both the diagonal and random layouts on every seed tried.

**Local improvement, global setbacks.** The paper's most interesting
observation is that N queens is a *conflicting system*. A swap that raises the
score among the three matched queens can lower the score of the whole board,
because the moved queens now clash with queens that were not looked at. The
global order degree therefore goes up and down rather than climbing steadily,
yet the run still converges. Kanada describes three stages: a rapidly changing
start, a long quasi-stationary phase in which the probability of the solution
state keeps growing, and a final stationary state. The paper's Figure 12 gives
a six-queens position with global score 14 out of 15 that no swap improves
globally, a *local maximum* of the global score. The single-catalyst rule
escapes it anyway, because the local score of the right three queens rises
from 2 to 3 while the global score falls from 14 to 13. Chemart's tests check
this exact position, and check that the global score falls during runs in which
no reaction lowers a local score.

**Catalysts trade reactions for effort.** Figures 9 and 10 of the same paper
show that adding catalysts biases the search toward better boards and cuts the
number of reactions, while the matching effort, and so the running time, grows
beyond two catalysts. Without a catalyst the system "does not stop even when a
solution is found": it is a random walk. Chemart's tests confirm both: three
catalysts need well under 60% of the reactions of one, and the no-catalyst rule
runs until its budget is spent.

**Colouring the USA, and frustration.** Kanada's FUZZ-IEEE'95 paper applies
CCM to four-colouring the mainland USA map, starting from a single colour, and
reports "a correct solution in every run". It found about 112 reactions and
4,406 attempts on average for the variable-catalyst rule with frustration.
Chemart's version of the rule, taken from Kanada's own demo program, gives
about 117 reactions and 2,500 attempts. The reaction count is checked in the
tests; the attempt count is not, because the paper ran "a slightly modified
version of the rule" that it does not describe. The paper also shows that the
variable-catalyst rule needs frustration: without it the rule can freeze in a
wrong colouring, which Chemart's tests reproduce.

**Larger graphs.** A 1996 manuscript extends the frustration method and tests
it on the DSJC random graphs used by Johnson et al. to benchmark simulated
annealing, with 125 or 250 vertices and edge probabilities from 0.1 to 0.9.
Annealed CCM colours DSJC125.1 with 6 colours in 0.2 s and with 5 in 13 s, and
DSJC250.5 with 29 colours in about 85 minutes on one CPU. Its times are
comparable to the annealing methods and GSAT in several cases, with failures
under 5% for well-chosen `f0`. Chemart does not reproduce this table. It offers
the same random-graph model, but on a random graph with DSJC125.1's size and
density the variable-catalyst run solves 6 colours in only 2 of 5 seeds (all 5
with 7 colours): the stopping rule, taken from the demos, gives up much earlier
than Kanada's C programs did. Nor does Chemart reproduce the manuscript's
parallel runs. Its claim that parallel reactions need almost no locking, since
clashes act as harmless noise, is noted in the book but not simulated.

**Beyond this entry.** Kanada's SMC'95 paper, which the book cites for the
annealing variant, actually introduces a different escape mechanism, CCM\*,
which composes rules at random ("randomised dynamic tunnelling"). Chemart does
not implement it. The frustration method comes from FUZZ-IEEE'95 and the 1996
manuscript.

## Further reading

- Kanada's CCM pages, with the papers, the demo applets (now in JavaScript)
  and short method notes: <https://www.kanadas.com/CCM/>
- Johnson, D. S., Aragon, C. R., McGeoch, L. A. & Schevon, C. (1991).
  Optimization by simulated annealing: an experimental evaluation. Part II,
  graph coloring and number partitioning. *Operations Research* 39(3),
  378–406. The source of the DSJC graphs and the annealing baselines.
