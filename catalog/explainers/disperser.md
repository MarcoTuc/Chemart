## Introduction

The chemical disperser is a load-balancing algorithm for computer networks,
written as a chemical reaction network. Each computer (node) holds a pile of
jobs waiting to be run. The jobs are treated as molecules, and each network
link as a catalyst that, now and then, picks up one job molecule on its side
and moves it to the node at the other end. A node with many jobs loses them
faster than a node with few, simply because more of its molecules meet a
catalyst per second. Left alone, the network ends with the same number of jobs
on every node: the average of what it started with.

It was proposed by Thomas Meyer and Christian Tschudin at the University of
Basel (Meyer, Yamamoto and Tschudin 2008; Meyer and Tschudin 2009) as part of
their "chemical networking protocols". Their idea was to design network
protocols out of molecule-like packets that react according to the law of mass
action, so that a protocol can be analysed with the tools of chemical
kinetics. The disperser was their first example. They present it as a chemical
alternative to *gossip* protocols for computing averages, such as Push-Sum, in
which each node repeatedly sends half of its value to a random neighbour. The
disperser does the same job, and its convergence proof "keeps on less than
half a page" (Meyer and Tschudin 2009).

Banzhaf and Yamamoto use it (book §17.3.1) to illustrate what they call the
engineering of communication protocols based on chemical dynamics. The point
of the example is not the averaging itself, which diffusion-based load
balancing has done since the late 1980s, but that writing the algorithm as
chemistry makes it safe and provable almost for free. Work is conserved by
every reaction, so a node's load can never go negative. And the reactions
translate directly into differential equations whose convergence and
stability can be checked.

It is a small *simulation model*: a fixed network of reactions on a given
graph, one reaction per directed link, simulated molecule by molecule. It
sits in the book's chapter on computing with artificial chemistries, next to
[Fraglets](fraglets.md), the string-rewriting language in which Meyer and
colleagues implemented the full protocol, and [computing with chemical
organizations](organization-computing.md). Fraglets is a general programming
language whose molecules are code; the disperser keeps only the core reaction,
with molecules that carry no information at all.

## How it works

### Molecules and reactions

There are two kinds of molecule. `X4` is one job waiting at node 4; the number
of `X4` molecules is node 4's load. `C1_4` is a catalyst sitting at node 1 and
responsible for the link from node 1 to node 4. Every undirected link carries
two catalysts, one at each end. The only reaction is

```
C_ij + X_i -> C_ij + X_j      for every directed link (i, j)
```

A job at node `i` meets the catalyst of the link towards `j`, and the job
reappears at `j`: it has been sent over the link. The catalyst is unchanged.
The number of catalysts never changes and is the same on every link, so the
catalysts only set how often each link is used.

The default network is the book's four-node example (figure 17.3): node 1 is
linked to 2, 3 and 4, and 2 is linked to 3. Its eight reactions are printed
below. Take node 1. It has three catalysts, `C1_2`, `C1_3` and `C1_4`, so
each of its jobs has three ways of leaving, while a job at node 4 has only
one, `C4_1 + X4 -> C4_1 + X1`.

### The reactor

The reactions are scheduled by the law of mass action: each happens at a rate
`k × c × x_i`, where `k` is the rate coefficient (the same for every link),
`c` the amount of catalyst and `x_i` the number of jobs at node `i`. The book
and the original work simulate this with Gillespie's stochastic simulation
algorithm (SSA), which draws individual reaction events at random times with
exactly these rates. The nodes are separate compartments, but because each
species name already says which node it belongs to, the whole network can be
simulated as one well-mixed vessel. This is how the book's reference program,
`Disperser.py`, does it, and what Chemart's network is.

### Why it converges

Adding up the flows into and out of node `i` gives one equation per node:

```
dx_i/dt = k c ( Σ_{j linked to i} x_j  −  deg(i) × x_i )
```

where `deg(i)`, the *degree*, is the number of links at node `i`. A node with
more links loses jobs faster, but also receives jobs from more neighbours. The
loads stop changing when each node holds the average of its neighbours, and on
a connected network that happens only when every node holds the same amount.
Since reactions only move jobs, the total never changes, so that amount is the
starting total divided by the number of nodes `n`: `x̄ = Σ x_i(0) / n`.

In matrix form the equation is `dx/dt = −k c L x`, where `L` is the
*Laplacian* of the graph: degrees on the diagonal, −1 for each link. This is
exactly a discrete diffusion process on the graph (the book calls it "a
chemical implementation of a discrete mass-conserving diffusion process on an
amorphous topology", citing Yamamoto et al. 2011). Meyer and Tschudin prove
stability from the eigenvalues of `L`, which are never negative. The smallest
is always 0, belonging to the direction "all nodes equal", which the conserved
total pins down; every other direction decays. The smallest non-zero
eigenvalue therefore sets how fast the last imbalance disappears: the better
connected the graph, the faster it balances.

## Using it

The default call returns the network of figure 17.3 with every coefficient 1,
one unit of each catalyst, and 1,000 jobs on node 4 in `net.initial_state`.
`net.extras["compartments"]` lists, for each node, its job species and the
catalysts located there, for example `'4': ['X4', 'C4_1']`. The network has no
random parts, so the `seed` changes nothing. Chemart generates the network; to
watch it run you simulate it yourself, as below.

#### Reproducing figure 17.3

The book's figure injects 1,000 jobs at node 4, adds 600 at node 2 at
`t = 10`, and removes 300 from node 2 at `t = 20`. This script runs the
Gillespie SSA on the generated network with those events, treating the
concentrations as molecule counts (one molecule per unit), as `Disperser.py`
does:

```python
import numpy as np

import chemart

net = chemart.generate_network("disperser", seed=1)
n = {s: round(v) for s, v in net.initial_state.items()}  # X4 = 1000, every C = 1
k = [r.rate["k"] for r in net.reactions]
changes = {10: ("X2", +600), 20: ("X2", -300)}  # the book's figure 17.3
rng = np.random.default_rng(1)

t, fired = 0.0, 0
for stop in (5, 10, 15, 20, 25, 30):
    while True:
        # propensity of C_ij + X_i: k * (number of C_ij) * (number of X_i)
        a = np.array([ki * np.prod([n[s] for s in r.reactants])
                      for ki, r in zip(k, net.reactions)])
        dt = rng.exponential(1 / a.sum())
        if t + dt > stop:  # memoryless: restart the clock at the stop
            t = stop
            break
        t += dt
        r = net.reactions[rng.choice(len(a), p=a / a.sum())]
        for s in r.reactants:
            n[s] -= 1
        for s in r.products:
            n[s] += 1
        fired += 1
    print(f"t={stop:2d}  " + "  ".join(f"X{i}={n[f'X{i}']:4d}" for i in "1234"))
    if stop in changes:
        s, d = changes[stop]
        n[s] += d
print(fired, "reactions")
```

```
t= 5  X1= 245  X2= 250  X3= 252  X4= 253
t=10  X1= 263  X2= 232  X3= 246  X4= 259
t=15  X1= 391  X2= 395  X3= 434  X4= 380
t=20  X1= 390  X2= 400  X3= 399  X4= 411
t=25  X1= 350  X2= 301  X3= 316  X4= 333
t=30  X1= 342  X2= 319  X3= 338  X4= 301
77241 reactions
```

The rows at `t = 10` and `t = 20` are printed just before each change. The
loads settle around 250 (1,000 / 4), then 400 (1,600 / 4), then 325
(1,300 / 4), with the random scatter of a few dozen jobs that the book's
figure also shows. The run takes about nine seconds; the number of events,
and so the time, grows with the number of jobs.

#### How fast, on which network

Integrating the equations instead (the same `S`-matrix recipe as on the
[Brusselator](brusselator.md) page) shows how the shape of the network sets
the speed. The table gives the time until every node is within 1% of the
average, with 1,000 jobs placed on one node, and the smallest non-zero
Laplacian eigenvalue (the *gap*):

| network | links | gap | time to within 1% |
|---|---|---|---|
| default, 4 nodes | 4 | 1.000 | 5.6 |
| line of 20 nodes, jobs at one end | 19 | 0.025 | 214.9 |
| ring of 20 nodes | 20 | 0.098 | 54.1 |
| star of 20 nodes, jobs on a leaf | 19 | 1.000 | 7.5 |
| complete graph of 20 nodes | 190 | 20.000 | 0.4 |

The graphs were built as lists, for example
`graph=[[i, i + 1] for i in range(1, 20)]` for the line, with
`initial_jobs={"20": 1000}`. On the default network, doubling `k` halves the
time (2.8) and halving `catalyst_concentration` doubles it (11.2): only their
product matters. On a disconnected graph each piece balances separately:
`graph=[[1, 2], [3, 4]]` with all jobs on node 1 ends at 500, 500, 0, 0.
Each of these integrations takes a fraction of a second.

The parameter table below lists the four parameters: the graph, the rate
coefficient, the catalyst amount and the initial loads. Node ids in
`initial_jobs` are strings.

## Results

**Convergence to the average, and back after a disturbance.** The published
result is a proof, not a measurement. Meyer and Tschudin (2009) write the
rate equation of every node, show that its only resting point on a connected
network is the equal share `x̄`, and show that this point is stable for any
topology, because the Jacobian of the system is minus the graph Laplacian.
Stable here means that after jobs are added or removed anywhere, the network
settles again at the new average. The book illustrates both claims with the
stochastic run of figure 17.3: about 250 jobs per node after the first
injection, about 400 after 600 more jobs arrive at node 2, and a lower average
after 300 are removed. It also notes that nodes far from the injection point
(nodes 2 and 3 after the injection at node 4; node 4 after the changes at node
2) converge more slowly. Chemart's test integrates the default network's
equations for 60 time units and checks that every node ends at 250 to within
0.01%. The stochastic run and the injection and removal events are not tested;
the script above reproduces them.

**Comparison with gossip.** Meyer and Tschudin compare the disperser with
Push-Sum. Both rely on a conserved quantity, and they report that the
convergence time is the same for suitably chosen parameters. The difference is
in how work moves: Push-Sum ships half of a node's value at once, the
disperser one job per reaction. Neither protocol survives the loss of
messages, but a lost packet costs the disperser one job, where it costs
Push-Sum half a node's value. The price is many more messages. This
comparison is qualitative in the paper; Chemart does not implement Push-Sum.

**The full protocol and what came after.** The complete disperser protocol was
implemented in Fraglets and is depicted in Meyer, Yamamoto and Tschudin
(2008) (book §17.3.1); the 2009 paper also points there for a variant in which
nodes need no information about their neighbours. Chemart implements only the
core reaction. The book places the disperser in a line of work by the same group
that relied on mass-action scheduling: a packet scheduler, a chemical
congestion-control protocol that behaves like TCP, and a distributed rate
controller based on enzyme kinetics, tested on a real network. None of these
is part of this entry.

**What it leaves out.** The book is explicit that the algorithm ignores
problems of practical load balancing: processors of different speed and
capacity, transmission delays, jobs arriving and completing continuously, and
nodes that steal or ignore jobs. It presents it as an illustration of the
concept of chemical programming for distributed algorithms. The book's own
figure is the only numerical result, and Chemart's network matches the
reference program behind it: same topology, one catalyst per link, all rates
1.

## Further reading

- Kempe, D., Dobra, A. & Gehrke, J. (2003). Gossip-based computation of
  aggregate information. *Proc. 44th IEEE Symposium on Foundations of
  Computer Science*, pp. 482–491. The Push-Sum protocol the disperser is
  compared with.
- PyCellChemistry, the Python package accompanying the book, which contains
  `Disperser.py`: <https://www.cs.mun.ca/~banzhaf/pycellchem-1.0/doc/index.html>
- Yamamoto, L., Miorandi, D., Collet, P. & Banzhaf, W. (2011). Recovery
  properties of distributed cluster head election using reaction–diffusion.
  *Swarm Intelligence* 5(3–4), 225–255 (book [941]). The book's reference
  for the disperser as mass-conserving diffusion.
