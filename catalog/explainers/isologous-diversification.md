## Introduction

Isologous diversification is a theory of cell differentiation proposed by
Kunihiko Kaneko and Tetsuya Yomo in the mid-1990s (Kaneko & Yomo 1994, 1997),
together with the simulation model built to demonstrate it. Its question is
how the cells of an organism, which all carry the same genome and start from
one fertilised egg, come to form distinct cell types (muscle, nerve, blood),
and how a cell's type is then passed on to its daughters. The usual answer is
a programme: genes switch on and off in response to signals or to a position
in the embryo. Kaneko and Yomo set out to show that a prototype of
differentiation can arise "even without implementing a programmed switching
process of genes" (1997, §1.1), from nothing more than cells that grow, divide
and talk to each other through a shared medium.

The picture is this. Every cell contains the same small network of catalytic
reactions, chosen so that the concentrations inside a single cell oscillate.
The cells all draw nutrient from, and leak chemicals into, one well-stirred
medium, so each cell's oscillation is pushed around by all the others. A cell
divides when it has consumed enough, and its two daughters get almost, but not
exactly, half of its chemicals each. While there are few cells they oscillate
in step and divide together: 1, 2, 4, 8 cells. Past some number the coupled
oscillators fall out of step, the tiny differences left by division are
amplified, and the cells settle into a few groups with clearly different
average chemical compositions. These groups are the cell types, and the
daughters of a cell of one type are of the same type. The name stresses the
point: *isologous* ("identical", in contrast with *homologous*, "similar")
units diversify through their interaction, not through errors or mutations.

Behind it lies Kaneko's earlier work on *globally coupled* dynamical systems,
in which many identical chaotic or oscillating elements, each driven by the
average of all the others, split spontaneously into clusters. The cell model
adds what those systems lacked: the number of elements grows, because cells
divide. The 1997 paper also answers a gap in Kauffman's 1969 proposal (see
[random Boolean networks](rbn.md)) that cell types are the different
attractors of one gene network: that picture needs cells to start from
different initial conditions, and does not say where those come from. Here the
cells choose their own initial conditions, through the interaction.

It is a simulation model: ordinary differential equations for the
concentrations inside each cell and in the medium, plus rules for division and
death. Banzhaf and Yamamoto describe it in the chapter on modelling biological
systems (book §18.5, "Cell Differentiation and Multicellularity"), mainly
through the *spatial* version of Furusawa and Kaneko (1998), in which cells sit
on a two-dimensional grid and stick together, and return to it in §18.6 as an
explanation of where positional information could come from. Chemart
implements the earlier, non-spatial model of Kaneko and Yomo (1997), in which
all cells share one medium. Its nearest neighbours in the catalog are the
[French flag model](french-flag.md), where a cell's fate is read off an
imposed morphogen gradient, and [random Boolean networks](rbn.md), where cell
types are attractors of a single isolated network; in isologous
diversification the types exist only because cells interact. The
[cellular Potts evo-devo models](cpm-grn-evodevo.md) of §18.6.1 also produce
cell types on a grid, but through an evolved gene network rather than
unprogrammed dynamics.

## How it works

### The chemicals of one cell

Each cell holds k + 1 chemicals, named `X0` to `Xk` (k = 8 by default).
`X0` is the **source**, a nutrient that the cell takes from the medium. The
others are unspecified cellular chemicals: the authors say they may include
metabolites as well as products of gene expression. The concentration of
chemical m in cell i is written x(m); the concentration of the same chemical
in the medium is X(m). The chemicals are not given any chemistry beyond the
network of reactions that connects them.

### The reaction network

Three kinds of reaction run inside every cell. Here are some from the default
network (seed 1):

```
X0 + X1 -> 2 X1  [mass-action k=1.0]
X1 + X4 -> 2 X4  [mass-action k=1.0 saturation=michaelis-menten saturated_species=X1 x_M=10.0]
X2 -> DF  [mass-action k=0.2]
```

- **Source paths.** Every chemical is made from the nutrient, catalysed by
  itself: `X0 + X1 -> 2 X1` at rate `e0 x(0) x(1)`. The more of `X1` a cell
  has, the faster it turns nutrient into more `X1`.
- **Catalytic paths.** Each chemical is converted into a few others (three by
  default, the `connections` of the network). `X1 + X4 -> 2 X4` is the path
  from `X1` to `X4`, catalysed by its own product `X4`: an autocatalytic path,
  as in the paper's networks. Its rate is `e1 x(4) x(1) / (1 + x(1)/x_M)`.
  The last factor is a *Michaelis-Menten* saturation: when the substrate
  `x(1)` is much smaller than the constant `x_M` the rate is simply
  proportional to `x(4) x(1)`, and when it is much larger the rate stops
  growing with `x(1)`. The paper writes the path from m to l catalysed by j as
  `Con(m, l, j) = 1`.
- **Paths to the division factor.** Some chemicals (four of the eight by
  default, here `X2`, `X4`, `X5`, `X6`) are consumed at rate `gamma x(l)` into
  the **division factor**, `DF`, which stands for everything a cell must make
  before it can divide: membrane lipids, DNA, ATP.

The network is drawn at random from the seed. That matters, because the whole
scenario rests on the cell's chemistry oscillating, and most random networks
do not. The paper found that with few paths per chemical the dynamics settles
to a fixed point; with many, one chemical takes over and the cell reduces to
source → one chemical → division factor (the "winner takes all" state of the
1997 paper's appendix 1); only a medium number of paths gives oscillations,
often switching between chemicals in turn. The authors kept only oscillating
networks. Chemart does the same: it draws up to `network_attempts` networks,
runs each as a single cell in its medium, and keeps the first one that still
oscillates after a transient and feeds the division factor.

### The medium, and how cells affect each other

Cells never touch. They exchange every chemical with the common medium in two
ways: **active transport**, `p (x(1) + ... + x(k)) X(m)`, so that a cell rich
in chemicals takes up more (the authors call the sum `x(1) + ... + x(k)` the
cell's **activity**), and **diffusion** through the membrane,
`D (X(m) - x(m))`. The medium has a volume `V`, measured in cell volumes,
which dilutes what the cells take and give. Nutrient flows into the medium at
rate `f` from a reservoir held at a fixed concentration (the parameter
`nutrient`, 40 by default), and the other chemicals wash out at rate `D_out`. Since all cells feed from the same nutrient, a cell that
takes more leaves less for the others: this competition is the coupling.

### Division and death

A cell divides when the division factor it has made since its birth, the
integral of `gamma x(l)` over its paths to `DF`, passes a threshold `R`. Its
daughters receive `(1/2 + ε)` and `(1/2 - ε)` of every concentration, where
the noise ε is drawn uniformly from a small range (`[-0.001, 0.001]` in the
paper). The volume of a cell is taken as constant except during the short
moment of division, which is why the concentrations are halved. A cell whose
chemicals `x(1) + ... + x(k)` sum to less than a starvation threshold `S`
dies, and its contents are released into the medium.

### What a cell type is

A cell type is not a species. It is the kind of oscillation a cell settles
into, read off the concentrations averaged over time. Two cells may differ at
every instant only because they oscillate out of phase; they are of different
types only if their *averages* differ. Chemart averages each cell's
concentrations since its last division, divides them by their sum to get a
composition, and groups cells whose compositions lie within `type_tolerance`
of each other (a Euclidean distance, as in Furusawa & Kaneko 1998).

### The stages

Kaneko and Yomo summarised their simulations as a sequence of stages, which
Chemart's analysis reports by number:

1. **Synchronous oscillation.** Up to a threshold number of cells, all cells
   are identical, oscillate in step and divide together, so the cell number
   runs 1, 2, 4, 8.
2. **Phase clustering.** Beyond it, the oscillations lose synchrony and cells
   group by phase; their averages are still the same.
3. **Fixed differentiation.** The averages themselves separate into a few
   distinct groups, which differ in composition, activity and oscillation
   period.
4. **Determination.** A cell's type is inherited by its daughters: the average
   composition of a daughter matches her mother's.
5. **Successive differentiation.** Types split again, giving a hierarchy.
   (Chemart counts this as stage 4 repeated.)

The formal specification below summarises the same model: the species are the
k + 1 chemicals plus `DF`, the reactions are the three kinds of path, and the
reactor integrates the equations for all cells and the medium, adding and
removing cells at divisions and deaths.

## Using it

The call above builds the network and also runs the cell society. For seed 1
the first random network drawn fails the screen and the second passes
(`network_attempts: 2`). A single cell with random starting concentrations
then divides synchronously, reaching 8 cells by about t = 4.3; `max_cells=8` stops further division and the run continues to
`t_max = 300`. What happened is in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["cells"], a["divisions"], a["deaths"]        # (8, 7, 0)
a["stage"], a["stage_name"]                    # (1, 'synchronous oscillation of identical cells (cell number 1, 2, 4, 8, ...)')
a["snapshot_spread"], a["average_spread"]      # (5e-06, 5e-06)
a["n_types"], a["types"][0]["composition"]     # (1, [0.107335, 0.0, 0.0, 0.0, 0.001569, 0.0, 0.891096, 0.0])
a["accepted"], a["single_cell_amplitude"]      # (True, 3.833245)
```

This is stage 1: eight cells whose compositions differ by at most 5 × 10⁻⁶.
The composition lists `X1` to `X8` as fractions of the cell's total; here the
cells hold about 89% `X7` and 11% `X1`, with the other chemicals near zero.
`snapshot_spread` compares the cells at the final instant and `average_spread`
their time averages; the stage is 2 when only the first exceeds
`type_tolerance`, and 3 or 4 when the second does. `recursivity` and
`inherited_fraction` compare each mother's averaged composition with her
daughters'; they matter from stage 3 on, and in stage 1 they only reflect that
the composition is still changing from one generation to the next.

Other useful extras: `extras["reaction_network"]["paths"]` lists the drawn
paths as `[m, l, j]` triples, `division_factor_paths` the chemicals that feed
`DF`, and `extras["interaction_law"]` the transport, diffusion and medium
equations with their constants and the final medium state. These couplings
are not reactions of one cell, so they are not in `net.reactions`.

**How often the default works.** The default depends heavily on the seed. Over
seeds 0 to 9 with default settings, the network screen succeeded within its
four attempts for seeds 1, 4, 6 and 9; six runs (seeds 0, 1, 2, 3, 8, 9)
reached 8 cells in stage 1, seed 7 stopped at 4, seed 4 at 2, and seeds 5 and 6
never divided. Seeds 4 and 6 show that passing the screen does not guarantee
division: the screen runs the single cell from different starting
concentrations than the real run, and the cell can end in another state. Seed 1's own
network, started from other random concentrations, divides once and then
settles on `X8`, which does not feed the division factor, so the two cells
never divide again. Raising
`network_attempts` helps the screen: with `network_attempts=16` seeds 2 and 3
found an accepted network (after 12 and 6 draws) and reached 8 cells; seed 0
found one after 6 draws but its cell never divided. A default run takes a few
seconds.

**The role of the unequal split.** With `split_noise=0.0` the daughters are
exact halves, and the eight cells of seed 1 are identical to machine
precision (`snapshot_spread` 0.0 instead of 5 × 10⁻⁶). Cells that start
identical stay identical, so the imbalance is what seeds any later
differentiation; the 1997 paper stresses that its size and mechanism do not
matter, only that some difference exists.

**A single cell.** `max_cells=1` integrates one cell in its medium without
division, the setting of the 1997 paper's Fig. 4 used to inspect a network's
oscillation.

**More cells.** `max_cells` raises the cap. It does not always matter: seed 1
stops at 8 cells even with `max_cells=32`, because the dominant chemical `X7`
does not feed the division factor and the medium's nutrient runs out. With
`max_cells=32` and default settings otherwise, seed 2 reached 16 cells (25 s),
and seeds 0 and 3 reached 32 cells (two and four minutes, on a busy machine).
All three were still in stage 1 at t = 300: the largest difference between
averaged compositions was 0.00005, 0.004 and 0.04, below the default
`type_tolerance` of 0.1. Runs grow slower with the number of cells, since every
cell adds 2k + 3 variables to one stiff system.

**The paper's parameters.** Most kinetic defaults are those of the 1997
paper's main simulation (its Fig. 6). Three are not: the division threshold
`R` is 100 (the paper's Fig. 18 simulations) instead of 2000, the death
threshold `S` is 0.01 instead of 0.05, and the medium volume `V` is 100 (the
value of the 1998 paper) instead of 1000. The paper's two published sets are

```python
fig6  = dict(medium_volume=1000, division_threshold=2000, death_threshold=0.05)
fig18 = dict(medium_volume=1000, nutrient=10, connections=2)   # R = 100, S = 0.01 as by default
```

with `max_cells` and `t_max` raised to the paper's 32 to 64 cells and several
hundred time units. These runs take minutes, not seconds, and are not covered
by the tests. One example: the Fig. 18 set with `max_cells=32, t_max=1000,
network_attempts=16` and seed 1 took 11 draws to find an oscillating network,
reached 32 cells, and ran for about two minutes; at t = 1000 the cells were
still in stage 1 (instantaneous compositions within 0.0004 of each other),
whereas in the paper's run with these parameters differentiation had
begun by 32 cells. The paper's networks (its Fig. 5) cannot be read off the figure,
so the network is still a random draw of the same shape; `paths` pins one by
hand.

**Variants.** `enzyme="quadratic"` replaces the saturating catalytic term with
`e1 x(m) x(j)²`, the enzyme term of Furusawa and Kaneko (1998), which is
exactly the mass action of `Xm + 2 Xj -> Xl + 2 Xj`. Only that term changes;
the rest of the model stays the 1997 one. `autocatalytic=False` draws each
path's catalyst at random instead of using its product.

## Results

### Kaneko and Yomo (1997): the scenario

The founding paper (Bulletin of Mathematical Biology, 1997), which extends a
simpler model of 1994 from rich and poor cells to successive cell types,
simulated networks of k = 8, 16, 32 and 64 chemicals with 2
to 6 paths per chemical, and reports that typical behaviour was common to
them. Its main example uses k = 8 with three autocatalytic paths per chemical.

- **Stages 1 and 2.** Cells divide together up to 8. Around the division from
  8 to 16 cells the oscillations desynchronise and the 8 cells split roughly
  into two phase groups, while their averages stay almost identical.
- **Stage 3.** Averaged compositions start to differ at about t = 280 (16
  cells), two groups are fixed by about t = 400 (32 cells), and a third group
  appears. One group has higher activity, takes up more of the source, and
  oscillates and divides faster. The authors stress that the difference in
  phase is "analogue" and easily reset by division, while the difference in
  averages is "digital": few, well-separated groups.
- **Stage 4.** The return map of a mother's average against her daughter's
  lies on the diagonal after about the 90th division, and the lineage keeps
  its colours for t > 400. The same type also appears from different branches
  of the lineage, as in the lineages of the worm *C. elegans*.
- **Stage 5.** By 64 cells (t = 940) each group splits again. The
  undifferentiated "red" cells give rise to red, green or blue cells, while
  green and blue only reproduce themselves; the authors liken the red cells to
  stem cells.
- **Rare chemicals.** The clearest differences between types, and the first
  differences at the onset of clustering, are in chemicals present at very low
  concentration.
- **Specialisation and tumours.** Cells that concentrate on few chemicals are
  more active and divide faster. With a stronger diffusion coupling
  (D = 0.2, R = 500) one cell of 32 appears at t ≈ 140 with `x(4)` around 3.6
  and the other chemicals near zero; its offspring divide faster, keep the
  trait and take over, and its return map loses recursivity. The authors call
  it a tumour-like cell and predict that tumour cells have reduced chemical
  diversity.
- **Cell death.** With a second network (two paths per chemical) and death
  switched on, the cell number fluctuates aperiodically around 32, with many
  cells dying simultaneously, which the authors compare with programmed cell
  death.
- **Parameters.** Weak nonlinearity delays differentiation (to around 128
  cells), strong nonlinearity brings it forward to around 8. Lowering `R` or
  raising `D` suppresses differentiation. The distribution of cell types is
  robust to noise, and recovers after cells of one type are removed.
- **Transplantation.** Determined cells placed among undifferentiated ones
  keep their type, but a population made mostly of one determined type partly
  dedifferentiates: cell memory lives in each cell, but needs the interaction
  to be kept.

### Furusawa and Kaneko (1998): rules and stability

A second paper with a modified model (20 chemicals, diffusion only, cells that
grow in volume and divide when it doubles) found the same scenario and made it
quantitative. From a single "type-0" cell, two cells change type when there
are 16 cells, and six types appear following a fixed rule: 0 → {0, 1, 2} and
1 → {3, 4, 5}. Types 1 to 4 cannot exist alone: a single cell started in one
of them returns to type 0, so they are sustained only by the other cells. The
rate at which type-0 cells differentiate depends on how many there are, which
holds the population (100 cells, division switched off) near
(n0, n1, n2) = (40, 30, 30); removing type-2 cells triggers more 0 → 2
differentiations until the distribution recovers. Over 100 developments from
random single cells to 300 cells, the number of type-2 cells fell into four
peaks (0, about 100, 150 and 220), that is, distinct stable "colonies". About
5% of random networks oscillated, and about 20% of those differentiated; with
added autocatalytic paths the figures rose to 40% and more than 20%.

### Later work described in the book

The book's account centres on the spatial version (Furusawa & Kaneko 1998, in
*Artificial Life*): cells on a grid take up nutrient, divide and adhere; a
cluster of identical cells forms, cells inside it differentiate once it
exceeds a threshold size, and a ring pattern of types emerges. Small clusters
can detach and found new colonies, giving a life cycle of multicellular
"organisms" that was not programmed. Later refinements (2002) found rings most
often but also stripes, and patterns robust to cell damage. Takagi and Kaneko
(2002) showed that predefined cells are not needed: in a multi-component
extension of the Gray-Scott reaction-diffusion model, self-replicating spots
form and differ in their internal dynamics, some oscillating and some at a
fixed point, like cell types. Furusawa and Kaneko (2006) showed that
differences in nutrient uptake between types create chemical gradients around
the cells, which could serve as positional information, and a 2009 survey
argues that the theory is consistent with experiments on induced pluripotent
stem cells. The book notes that the actual mechanism of cell differentiation
remains debated.

### What Chemart reproduces

The tests check the model rather than the late stages. They verify that every
reaction is one term of the paper's equation (1) and that the simulated rate
equations match the paper's equations (1) to (7), written out term by term,
for both enzyme forms; that the quadratic network integrates exactly as mass
action, and the saturating one reduces to mass action when `x_M` is large;
that two chemicals feeding each other show the winner-takes-all behaviour of
the paper's appendix 1; that the drawn network has the prescribed shape; that
the screened network oscillates as a single cell; that a cell divides exactly
when the integral of equation (8) reaches `R`; that the split conserves the
mother's chemicals with an imbalance no larger than the noise; that a starved
cell dies and releases its chemicals into the medium; and, in slow tests,
stage 1 (8 cells, 7 divisions, one type) and that the cells stay exactly
identical when the split is exactly equal.

Chemart does not reproduce stages 2 to 5, and nothing in the tests asserts
them: the phase clustering, the fixed differentiation into types, their
inheritance, the successive splitting, the role of rare chemicals, the
tumour-like cells, the collective deaths, the recovery of the type
distribution after removal, and the transplantation experiments. The
implementation decisions say these need the paper's 32 to 300 cells, where a
single run costs minutes; the runs above reached 32 cells without leaving
stage 1. The spatial model the book illustrates (grid, adhesion, ring and
stripe patterns, detaching clusters) is not implemented, nor is the full 1998
model.

## Further reading

- Kaneko, K. & Yomo, T. (1994). *Physica D* 75, 89–102. The earlier,
  simpler model with cell division and dynamic clustering that the 1997 paper
  extends.
- Kaneko, K. & Yomo, T. (1999). Isologous diversification for robust
  development of cell society. *Journal of Theoretical Biology* 199(3),
  243–256, doi:10.1006/jtbi.1999.0952.
- Furusawa, C. & Kaneko, K. (2002). Origin of multicellular organisms as an
  inevitable consequence of dynamical systems. *The Anatomical Record*
  268(3), 327–342 (book ref. [302]).
- Furusawa, C. & Kaneko, K. (2006). Morphogenesis, plasticity and
  irreversibility. *International Journal of Developmental Biology* 50,
  223–232 (book ref. [303]).
- Furusawa, C. & Kaneko, K. (2009). Chaotic expression dynamics implies
  pluripotency: when theory and experiment meet. *Biology Direct* 4, 17 (book
  ref. [304]).
- Takagi, H. & Kaneko, K. (2002). Pattern dynamics of a multi-component
  reaction diffusion system: differentiation of replicating spots.
  *International Journal of Bifurcation and Chaos* 12(11), 2579–2598 (book
  ref. [840]).
