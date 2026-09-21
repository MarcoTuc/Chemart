## Introduction

A Self-Organising Assembly System (SOAS) is a factory assembly line that
designs itself. The shop floor holds a pool of robotic modules: robots, axes
(motorised slides that move along one direction), grippers, part feeders,
positioning devices and a human operator. A product order arrives, describing
the tasks needed to build the product. The modules then group themselves into
teams, called *coalitions*, each able to carry out one task, and the teams
together make an assembly line for that product. Nobody picks the modules or
decides where they go.

Frei, Di Marzo Serugendo and Şerbănuţă proposed SOAS in 2010, as part of research on *evolvable assembly systems*, a term
coined by Onori (2002) for assembly systems that adapt as the product and the
assembly process change. Their motivation, stated in the 2010 paper, is agile
manufacturing: companies increasingly make many product variants in low
volumes, and rebuilding and reprogramming an assembly line for each is manual,
work-intensive and error-prone. The running example is an adhesive tape roller
dispenser, assembled on a carrier from a body case and a tape roll (2010), and
later from two body-case parts, a tape roll and a screw (2012).

What makes SOAS an artificial chemistry is how the self-organisation is
written down. The authors use the [Chemical Abstract Machine](cham.md) (CHAM),
a formalism in which a computation is a solution of molecules and rules say
which molecules may react. In SOAS the solution is the shop floor, the
molecules are the modules and the order, and a reaction lets a module join a
coalition when it offers a skill the task needs and has a physical interface
that fits. The rules are executable in Maude, a rewriting-logic language, so
one run produces candidate assembly lines, and the trace of rules that built
each line shows how it satisfies the order. The book puts it this way: the
system "is able to find feasible configurations that can be proved correct".

Banzhaf and Yamamoto describe SOAS in one paragraph of the chapter "Beyond
Chemistry and Biology" (book §20.1), next to [mechanical self-assembly](mechanical-self-assembly.md).
The two differ in kind. Hosokawa's magnetic triangles really collide and stick,
and their yields follow chemical kinetics; SOAS has no kinetics at all. It is a
constructive rule system in which the interesting output is *which* structures
can be built, not how fast. Among the rule formalisms, it is the catalog's
worked application of the [CHAM](cham.md), much as the
[Chemical Casting Model](ccm.md) applies chemical rules to puzzles.

## How it works

### The molecules

There are four kinds of molecule.

- A **manufacturing resource agent** (MRA) is one module. It has a type, a set
  of *skills* it offers (a robot offers several `move` skills, each with a
  direction and a range in millimetres; a gripper offers `open-close` with its
  opening), the skills it *requires* from a partner, and its physical
  *interfaces*. An interface has a shape (circular, square, triangular, and so
  on) and a sign: `+` means the module is held, `-` that it holds another. A
  gripper has a `circular+` interface; a robot that can carry a gripper has a
  `circular-` one.
- The **order**, called the *generic assembly plan* (GAP), is a list of tasks.
  Each task names the parts it handles, where it starts, and the operations it
  needs, for example `pick&place` for "take this part and put it on the
  product".
- A **coalition** is a group of MRAs formed for one task. It records the skills
  its members provide, the operations the task still requires, and the
  interfaces still free. A coalition that requires nothing more is *complete*.
- An **assembly line** assigns complete coalitions to tasks. When every task
  has one, the line can build the product.

Before any reaction, the parts sharpen the tasks. A task that picks up the
body case, which is gripped at two points 60 mm apart by a two-finger gripper,
also requires `grip(gripper-type=2finger,range=60.0)`. A two-finger gripper
offers `grip` with its opening as range. A task that starts at "the feeder"
also requires a feeder of the right part.

### The rules

The papers name six rule types. The first four build coalitions:

1. **interface compatibility**: two modules connect only through a free
   interface of the same shape and opposite signs, and both interfaces are then
   occupied;
2. **composition patterns**: a short list of typical module combinations
   (`move & grip`, `feed & move`, and four more) that a module and a coalition
   must match before they combine;
3. **composite skills**: some skills decompose into simpler ones, so a
   required `pick&place` becomes a `move` plus a vertical `move`;
4. **task-coalition matching**: a module that offers a skill a task requires
   starts a coalition for that task.

A skill offered matches a skill required when it has the same type and every
detail the requirement names, with a range at least as large. The other two
types lay the line out on the floor (type 5) and turn the order into concrete
robot movements, the *layout-specific assembly instructions* (type 6).

### A worked example

The two-part dispenser (`system="tape-roller-gap1"`) needs four tasks: load
the carrier, place the body case, place the tape roll, unload. Its second task
is `pick&place` on the body case, which starts at a named feeder, so no feeder
has to join. Here are two reactions from the Chemart network:

```
gap1 + g1 -> g1 + gap1 + gap1->t(2)[g1]
r3 + gap1->t(2)[g1] -> r3 + gap1->t(2)[g1] + gap1->t(2)[g1,r3]
```

In the first, gripper `g1` meets the order `gap1`. It offers
`grip(gripper-type=2finger,range=70.0)`, and the task requires a grip of
60 mm, so `g1` starts a coalition for task `t(2)` (rule type 4). The new
molecule records what is left to do:

```
task=gap1->t(2); coalition=g1; provided-skills=grip(gripper-type=2finger,range=70.0);
required-ops=move, move(direction=vertical); open-interfaces=circular+
```

The task's `pick&place` has been decomposed (type 3) into a `move` and a
vertical `move`, which a gripper cannot supply. In the second reaction robot
`r3` joins. It offers three rotational moves, one of them vertical, so it
supplies both missing skills; `move & grip` is a composition pattern (type 2);
and the robot's `circular-` interface fits the gripper's `circular+` (type 1).
The new coalition `gap1->t(2)[g1,r3]` has `required-ops=none`: it is complete,
and its remaining free interfaces are the robot's triangular, straight and
diamond ones. The 2010 paper walks through exactly this step (its figures 24
and 26).

In both reactions the reactants reappear among the products: the modules, the
order and the old coalition act as catalysts, because a coalition is only a
proposal, and `g1` could just as well join robot `r1` in another proposal. A module is used up only when an assembly
line adopts a coalition containing it. That is the last rule, **assign**: a
complete coalition fills an open task of a partial line, provided it shares no
module with the coalitions already in the line and every task still open can
still be filled. The human operator stands for any number of people, so
coalitions of the human are assigned all at once and never block anything.

### The reactor

In the papers the rules fire concurrently, in no set order, until nothing more
can react. Maude's search command explores every possible order. Chemart
computes the same end result directly: every coalition that can be built from
the initial solution, and every assembly line that can be built from the
order. The layout (type 5) is then computed for each finished line; it is a
fixed function of the line, not a reaction. Type 6 is not implemented.

## Using it

The default run is the four-part tape roller of the 2012 paper
(`system="tape-roller-gap2"`): 19 modules and a six-task order. It takes about
a second and a half. The 582 species are the order, the modules, 165
coalitions and the partial and complete lines; no reaction has a rate. The
results are in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["coalitions"]          # {'t(1)': 1, 't(2)': 44, 't(3)': 44, 't(4)': 44, 't(5)': 31, 't(6)': 1}
a["n_assembly_lines"]    # 48
a["complete_coalitions"]["t(5)"]   # ['f4+g4+r3', 'f4+g4+r1', 'f4+g4+r2', 'a1+a4+a7+f4+g4']
a["assembly_lines"][0]["assignment"]
# {'t(1)': 'unassigned-human', 't(2)': 'f1+g1+r3', 't(3)': 'f2+g2+r1',
#  't(4)': 'f3+g3+r2', 't(5)': 'a1+a4+a7+f4+g4', 't(6)': 'unassigned-human'}
```

Module names follow the authors' specification: `r` robots, `a` axes, `f`
feeders, `g` grippers, `pd` positioning devices, and `unassigned-human` the
operator. Each entry of `assembly_lines` also holds the line's species id and
its `layout`, a list of conveyors and stations with positions in millimetres.
`unassignable_tasks` lists tasks that no complete coalition can fill.
`extras["reaction_rules"]` labels every reaction with its rule; the default
network has 37 `initiate`, 232 `join`, 1 `assign-generic` and 1,056 `assign`
reactions.

**Removing a module.** `unavailable` takes modules off the floor, as in the
papers' resilience scenarios:

```python
net = chemart.generate_network("soas", unavailable=["g4"])
net.extras["analysis"]["unassignable_tasks"], net.extras["analysis"]["n_assembly_lines"]
# (['t(5)'], 0)       no screw-driver gripper: the screw task cannot be done
net = chemart.generate_network("soas", system="tape-roller-gap1", unavailable=["r1"])
net.extras["analysis"]["n_assembly_lines"]      # 6  (12 with r1)
```

**Your own shop floor.** With `system="custom"`, `mras`, `gap` and `parts` take
the modules, the order and the parts as text, one per line; the notation is in
the parameter table below and in the module docstring. Here one robot, two
grippers and a feeder face a one-part order. The narrow gripper opens 30 mm,
less than the 50 mm the box needs, so it is never used:

```python
mras = """\
rb robot | move(subtype=linear,direction=horizontal,range=200.0), move(subtype=linear,direction=vertical,range=100.0) | | circular-, triangular-
gw gripper/2finger | open-close(range=80.0) | | circular+
gn gripper/2finger | open-close(range=30.0) | | circular+
fd feeder | feed(subtype=feeds(box)) | | triangular+
unassigned-human human | load, unload | |
"""
gap = """\
demo 1
t(1) | other | tray | (0,0,0) | load
t(2) | pick&place | b | feeder | pick&place
t(3) | other | tray+b | | unload
"""
parts = """\
b box | 2finger | (0,0,5) (0,50,5)
tray carrier | |
"""
net = chemart.generate_network("soas", system="custom", mras=mras, gap=gap, parts=parts)
[x["line"] for x in net.extras["analysis"]["assembly_lines"]]
# ['demo{t(1)=unassigned-human;t(2)=fd+gw+rb;t(3)=unassigned-human}']
```

The robot needs a `triangular-` interface to hold the feeder; without it the
pick-and-place task has no complete coalition. `patterns` replaces the six
composition patterns, and `floor` sets the shop-floor size for the layout.

## Results

**The assembly system designs itself.** The 2012 paper runs the complete
specification, all six rule types, on the four-part dispenser. Each solution
assigns a coalition to every task and places the coalitions on the floor,
linked by conveyors: the line runs east, turns north with two corner conveyors
at the end of the floor, and comes back west, where the screw is inserted and
an operator unloads the product (its appendix B). Chemart reproduces this
layout position by position for the first solution; the test compares every
conveyor and station.

**48 lines.** The same paper says in section 5.5 that "Maude generates 47
possible solutions in just a few seconds", and a few lines later that "Maude
generates 48 solutions for this GAP", describing a 25th solution in which
robots `r3`, `r1` and `r2` swap tasks. Running the authors' published Maude
specification gives 48, and so does Chemart. The count has a simple reading:
four "workers" (robots `r1`, `r2`, `r3`, and a robot built from the three axes
`a1+a4+a7`) are permuted over the four pick-and-place tasks, 24 ways, times 2
ways of sharing the two body-case feeders. Chemart's tests check the 48, the
first and 25th solutions, the 165 coalitions (30 complete) behind them, and,
for the two-part dispenser, 46 coalitions and 12 lines, all counts taken from
the Maude specification.

**Only compatible modules combine.** A gripper is held by a robot or an axis
on a circular interface. The three-axis worker forms because the middle axis
`a4` needs a base axis, which `a1` supplies, and holds `a7`, which holds the
gripper. The 55 mm gripper `g2` can never take the body case, which needs
60 mm, so it always goes to the tape roll; positioning devices offer nothing
the tasks need and never join. The tests check each of these restrictions, and
that removing the `feed & move` pattern leaves every feeder task of the
four-part order without a coalition.

**Correct by construction.** In every line each task has a complete coalition,
no module works at two stations, and each feeder feeds the part its task
takes. Chemart's tests check all three for every one of the 48 lines. The
authors' argument for correctness is the trace: any solution Maude gives "can
first be traced (to check that each rule applied correctly) and validated"
(2012). The 2010 paper traces one run of the two-part order: "t1 by human1, t2
by coalition r3 - g1, t3 by r1 - g2, and t4 by human2". Chemart contains this
line and the reactions and open interfaces of the paper's coalition figures;
it uses the 2012 modules, where one generic operator replaces the two humans.

**Incomplete termination.** The 2010 paper checks its CHAM against
Wermelinger's (1998) three properties for a correct CHAM specification. The
first is that the CHAM terminates, and it notes that it may terminate "because
not all the required skills are provided in sufficient number, which means
that one or several tasks remain open". Chemart shows both ways this can
happen. Without the screw-driver gripper `g4` the screw task has no coalition
at all. Without robot `r1` every task still has coalitions, but only three
workers remain for four pick-and-place tasks, so no line exists. The `assign`
rule looks ahead, so a partial line is never extended into a dead end.

**Resilience.** The 2010 paper argues that if a module becomes unavailable
during design, the process starts afresh with what is left. Chemart's
`unavailable` parameter runs this: the two-part order still yields 6 lines
without `r1`, and the tests check it.

**What the papers do not claim.** The book says the assembly line is
"optimized for a given manufacturing order"; the 2012 paper says the model
"only searches for a viable solution, without giving any importance to
performance characteristics or efficiency", and the 2010 paper calls itself a
feasibility study. Book reference [297], Frei and Whitacre (2012), discusses
degeneracy and networked buffering as routes to evolvability in agile
manufacturing systems; Chemart does not model it. The book adds that the work
is "still at a prospective research stage and not ready for real industrial
deployment".

## Further reading

- Onori, M. (2002). Evolvable assembly systems: a new paradigm? *Proceedings
  of the 33rd International Symposium on Robotics (ISR)*, 617–621. The origin
  of the term, as cited by the 2010 paper.
- Wermelinger, M. (1998). Towards a chemical model for software architecture
  reconfiguration. *IEE Proceedings - Software* 145(5), 130–136. The three
  properties used to check the SOAS CHAM.
- Frei, R. & Whitacre, J. (2012). Degeneracy and networked buffering:
  principles for supporting emergent evolvability in agile manufacturing
  systems. *Natural Computing* 11(3), 417–430. doi:10.1007/s11047-011-9295-4
