## Introduction

EVOLVE is a family of computer ecosystems built by the biologist and computer
scientist Michael Conrad and his students between 1969 and the late 1990s. Each
model is a small simulated world of organisms. Each organism has a genome, gets
matter and energy from its surroundings, reproduces with mutation and dies, and
evolution is left to run. Conrad was not trying to optimise anything. He wanted
to see which ecological and evolutionary patterns come out of such a world,
where no fitness score is handed out and an organism survives only if it pays
its own way.

The line starts with Conrad's 1969 Stanford thesis and a 1970 paper with Howard
Pattee. In that model organisms live in a one-dimensional world, collect
materials, reproduce by *conjugation* (exchanging genetic material, as some
bacteria do) and die, and their materials are recycled so that the total amount
is conserved. The book reports that its usual outcome was a homogeneous, very
efficient population in which evolution stopped. The EVOLVE series proper
followed: EVOLVE II (Conrad and Strizich, 1985), EVOLVE III (Rizki and Conrad,
1985) and EVOLVE IV (Brewster and Conrad, 1998). In these, organisms live in a
two-dimensional world and obey *energy and matter conservation*. They harvest
light, secrete chemicals, and spend energy on metabolism and reproduction.
Their enzymes are strings, and what an enzyme does is looked up in a table by
how well the string matches, with a middle *critical section* counting most,
like the binding site of a real protein.

The result the book singles out is ecological. Seed the world with
*autotrophs*, organisms that make their own organic matter from inorganic
matter and light, as plants do. Then *scavengers*, which live by breaking down
the remains the autotrophs leave, arise by mutation, go extinct and come back.
A new way of making a living, an *ecological niche*, has emerged that nobody
programmed. A 2021 review of digital evolution (Vostinar et al.) calls EVOLVE
"one of the earliest examples of digital evolution".

Chemart's entry is a **simulation model**, and a **reconstruction**. None of
the EVOLVE papers is openly available, so the model is built from the book's
one-paragraph description (§8.2.3, "Some Early Models", in the chapter on
coevolution and ecologies) and the papers' published abstracts. The string
matching, the trophic economy and every number are Chemart's reading of that
paragraph, not Conrad's code. Its nearest neighbours in the catalog are the
other ecosystems of §8.2.3. [Dorin and Korb's ecosystem](dorin-korb-ecosystem.md)
also keeps a closed matter inventory and an energy budget, but builds organisms
out of bonded atoms rather than genomes. [Ecolab](ecolab.md) evolves the
coefficients of population equations and has no individuals at all.
[Urdar](urdar.md) makes the metabolism itself a computation. The
assembly-language ecosystems [Tierra](tierra.md) and [Avida](avida.md) evolve
programs rather than enzyme strings. Conrad's other entry,
[the lock-and-key enzymatic processor](conrad-enzymatic.md), uses the same
shape-matching idea for computation instead of ecology.

## How it works

### The world

The world is a square grid that wraps around at the edges (a *torus*), 10 × 10
cells by default. A cell holds at most one organism and two stocks of matter,
counted in whole units:

- `mineral`: free inorganic matter;
- `organic`: organic matter, meaning secretions and dead remains. Bodies are
  made of it: a living organism holds `body` units (2 by default).

Matter is never created or destroyed; it only moves between these two stocks
and the organisms' bodies. Energy comes from one source, sunlight, which falls
on every cell each step. An organism can reach the matter in its own cell and
in the four cells that touch it.

### Genomes, enzymes and the table

A genome is a string of genes written in the letters `a`, `b`, `c`, `d`. Each
gene is one *sensitivity* letter followed by an eight-letter *protein string*.
The sensitivity letter sets how many letters of the protein a mutation of that
gene changes: `a` means 1, `d` means 4. This follows EVOLVE II, whose abstract
states that "the magnitude of phenotypic change resulting from mutation is
itself a property of the gene".

A protein's function is fetched from a table of four reference strings:

```
light      abcdabcd   harvest sunlight
fix        bdbdbdbd   turn mineral matter into organic matter (autotrophy)
respire    cbadcbad   break organic matter back down to mineral (scavenging)
replicate  dddddddd   the machinery of reproduction
```

The protein is compared letter by letter with each reference. Letters 3 to 5
form the critical section, and a match there counts three times as much as a
match elsewhere. The score is scaled to between 0 and 1. The best-scoring entry
gives the protein its function, and the score is how efficiently the enzyme
does it. A protein whose best score is below 0.55 has no function at all. The
book gives only the principle, a table and a weighted critical section; the
four functions, the reference strings and the threshold are Chemart's.

Here is a founder genome from the default run, split into its six genes:

```
c aacdabcd   -> light      0.93
c bdbdbdbd   -> fix        1.00
c dbdddddd   -> replicate  0.93
c cdaababb   -> none       (best score 0.36)
c ccbacbdc   -> none       (0.36)
c ccbbbcba   -> none       (0.50)
```

The light gene differs from `abcdabcd` in its second letter, outside the
critical section, so it scores (3 × 3 + 4) / 14 = 0.93. If that one mismatch
had fallen inside the critical section the score would be (3 × 2 + 5) / 14 =
0.79. The founder has no respire gene, so it is a pure autotroph. Every founder
is built this way, so scavenging can only arise by mutation.

### One step of life

Each step, every living organism takes its turn, in random order:

1. **Harvest light.** It gains its light efficiency × the incident light
   (4 units by default), rounded down. The rest is lost.
2. **Scavenge**, if it has a respire enzyme and organic matter is within
   reach. It turns one organic unit back into a mineral unit and recovers part
   of the energy stored in organic matter.
3. **Fix**, if it has a fix enzyme and enough energy. It spends `fix_cost`
   (3) energy to turn one mineral unit into an organic unit, which it
   secretes into its cell. That energy is banked in the organic stock, where a
   scavenger can later recover it. An autotroph fixes only while less than one
   body's worth of organic matter is within reach, so it builds a store and
   then stops.
4. **Pay maintenance** (1 energy). An organism that cannot pay starves, and
   one older than `max_age` (40 steps) dies of old age. Its body goes back to
   the organic stock, together with the energy it still held.
5. **Reproduce**, if it has enough energy, a free neighbouring cell and a
   body's worth of organic matter within reach. It then succeeds with
   probability equal to its replicate efficiency. The offspring's genome is a
   mutated copy of the parent's. With probability 0.1, if a neighbour is
   present, the reproduction is *conjugative*: before mutation, the parent's
   genome is cut at a gene boundary and completed with the neighbour's, a nod to
   Conrad and Pattee's conjugation.

After all organisms have acted, a quarter of each cell's matter (rounded
down) spreads equally to its four neighbours.

### What a run produces

Each of these acts is written down as a reaction. The organism appears on both
sides as a catalyst: it carries out the reaction but is not used up by it.
Reactions from the default run:

```
ge30eb23c + mineral -> ge30eb23c + organic                              fix
g668b6a94 + organic -> g668b6a94 + mineral                              respire
gd689176b + 2 organic -> gd689176b + g3a08f8dd                          reproduction
gbcd53d11 + 2 organic + g2d905df8 -> gbcd53d11 + g9abbc715 + g2d905df8  conjugation
g0ea189f7 -> 2 organic                                                  death
```

A species name like `ge30eb23c` is `g` plus a hash of the genome. Every
distinct genome is a separate species, so the species set grows as mutants
appear, which is what *constructive* means in the table above. Each distinct
reaction is listed once, with how often it fired. Every reaction balances in
matter: a genotype species stands for a body of 2 organic units, so in
`G + 2 organic -> G + O` the two units become the offspring's body.

This is how scavengers appear. Here is a genotype from the default run whose
old light gene mutated into a respire gene:

```
c abadddcd   -> respire    0.57   (was the light gene)
a bdbabdbc   -> fix        0.71
d ddddddbd   -> replicate  0.93
```

It can no longer harvest light, but it can live off organic remains. Because
it can also fix, Chemart classes it as a *mixotroph*: an organism that both
fixes and scavenges.

The formal specification below summarises this: the species are the genotypes
plus the two matter stocks, the reactions are the five kinds of event, and the
reactor is the grid.

## Using it

The default run is not a published experiment; no EVOLVE parameters are
available. It is a small world set up so that the book's phenomenon can
appear in a couple of seconds: 8 autotroph founders, 120 steps, and little
enough matter (3 mineral units per cell, 316 units in all) that recycling
matters. It takes under two seconds.

What happened is in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["event_counts"]
# {'conjugation': 26, 'death': 207, 'failed_replication': 49, 'fix': 759,
#  'reproduction': 273, 'respire': 500}
a["genotypes_seen"], a["survivors"], a["final_guilds"]
# (178, 100, {'autotroph': 84, 'mixotroph': 5, 'none': 11})
a["scavenger_episodes"]
# [{'from': 24, 'to': 120, 'peak': 10}]
a["history"]["scavengers"][::10]
# [0, 0, 0, 1, 1, 3, 6, 10, 10, 8, 8, 5, 5]
```

The eight founders fill the grid by about step 30. The first organism able to
scavenge appears at step 24 and the scavengers stay until the end, peaking at
10. A *guild* is an organism's way of life, read from its enzymes: `autotroph`
(fixes), `scavenger` (respires), `mixotroph` (both), or `none` (neither: it
lives on light and builds its offspring from organic matter others made).
`history["scavengers"]` counts every organism with a respire enzyme, mixotrophs
included. `history` also holds, for every step, the population, the matter
stocks, the number of genotypes and the mean enzyme efficiencies.
`net.extras["phenotypes"]` gives each genotype's efficiencies, mean
sensitivity and guild.

The energy accounts are in `net.extras["energies"]["ledger"]`. In this run
42,528 units of light fell on the grid and 38,158 were harvested. Two balances
hold exactly: what organisms took in equals what they spent, deposited at death
or still hold (52,882 on each side), and the energy stored in organic matter
equals what was put in minus what was taken out. Most of the energy put into
organic matter, the store scavengers draw on, is what dead organisms leave
behind (28,553 units); fixing adds only 2,277.

**Watching the niche come and go.** The book's full claim, that scavengers go
extinct and reappear, needs longer runs. With `steps=600` (about a second per
run), seeds 0 to 6 give:

```python
net = chemart.generate_network("evolve-series", seed=4, steps=600)
net.extras["analysis"]["scavenger_episodes"]
# [{'from': 137, 'to': 353, 'peak': 8}, {'from': 382, 'to': 479, 'peak': 3},
#  {'from': 505, 'to': 600, 'peak': 10}]
```

Seeds 1, 4 and 5 show extinction and return (in seeds 1 and 5 the gap is only
3 and 6 steps). Seeds 0, 2, 3 and 6 keep scavengers from their first
appearance to the end. At the default 120 steps, seed 4 produces no scavengers
at all.

**Selection on the match.** Seed the world half with a perfect autotroph and
half with a copy whose light and fix critical sections have been scrambled
(efficiencies 0). Pass them as `founder_genomes=[good, bad]`, where `good` is
`"cabcdabcdcbdbdbdbdcddddddddcaaaaaaaacaaaaaaaacaaaaaaaa"` and `bad` is
`"cabdabbcdcbdcacdbdcddddddddcaaaaaaaacaaaaaaaacaaaaaaaa"`. In seeds 1 to 6 the
scrambled genome leaves no survivors, and the mean light + fix efficiency rises
from 1.0 to between 1.70 and 1.87.

**Mutation rate.** `mutation_rate` is the chance per gene per reproduction that
the gene is mutated. Over seeds 1 to 6 the mean fix efficiency at step 120 is
0.82 at the default 0.08, 0.94 at 0.02 and 0.99 at 0. With `mutation_rate=0`
no scavenger ever appears.

**Other settings.** `light_input=0` switches the sun off. `light_period` and
`light_swing` make the light vary as a sine wave, the ingredient of EVOLVE
III's constant-versus-variable experiment. `founder_genomes` also lets you
transplant an evolved genome, taken from `net.species[i].structure`, into a new
run. `width`, `height` and `steps` set the world's size and the run's length.

## Results

### What the published EVOLVE models found

Only the abstracts and the book's paragraph could be read, so this is what they
state, without the papers' numbers.

- **Early model (Conrad and Pattee, 1970).** The book reports that the usual
  outcome was a homogeneous population of organisms that used materials very
  efficiently, an ecosystem that was collectively stable, and an end to
  evolution.
- **EVOLVE II (Conrad and Strizich, 1985).** Organisms compete for a limited
  food supply in a changing environment. The abstract reports three findings.
  Two lineages with "distinctly different survival strategies" evolved and
  coexisted. Organisms "developed a resistance to phenotypic change in response
  to mutation in slowly varying environments". And traits that favour the
  individual at the expense of reproduction could still change under mutation,
  which the authors read as gene structures that benefit the lineage rather
  than the individual offspring.
- **EVOLVE III (Rizki and Conrad, 1985).** Population, organism and genetic
  structure are modelled as separate, replaceable levels. In the experiment the
  abstract describes, populations raised in a constant environment usually beat
  populations raised in a variable one when both were moved into a variable
  environment early in their development, and the reverse held later. The
  authors note that this agrees with experiments on laboratory microcosms.
- **Symbiosis in EVOLVE III (O'Callaghan and Conrad, 1992).** Parasitic feeding
  was added. The model then showed obligate and facultative parasitism,
  transient mutualism, speciation and mass die-offs, and total biomass rose
  with symbiotic activity. Stable mutualism was impossible, because organisms
  could not feed on other organisms and on the environment at the same time.
  Vostinar et al. (2021) review this as early work on symbiosis in digital
  evolution. They also describe EVOLVE as having light intensity and
  temperature in the environment and genomes of up to 40 genes, each up to 200
  bases long.
- **EVOLVE IV (Brewster and Conrad, 1998, 1999).** Organisms interact by
  exchanging metabolites and by changing their environment. The 1999 abstract
  reports that "niche formation occurs in the model".
- **The book (§8.2.3).** Starting from autotrophs alone, scavenger populations
  "emerge, go extinct, and reappear", showing the emergence of a niche for
  decomposing the autotrophs' remains.

### What Chemart reproduces

Each item below is one of the entry's recorded phenomena, with what the tests
in `tests/chemistries/test_evolve_series.py` check.

- **The scavenger niche (the book's result): partly.** Scavenging arises by
  mutation from pure autotrophs in the default run, and the tests check this.
  Extinction and return need 600 steps and appear in 3 of the 7 seeds tried.
  A slow test checks it for one seed (seed 4, three episodes).
- **Coexisting strategies (EVOLVE II): in a loose sense.** The tests check that
  autotrophs and organisms with a respire enzyme are both alive at the end of
  the default run. This is Chemart's reading of "two species"; the paper's
  lineages cannot be compared.
- **Exact matter conservation: yes, by construction.** Every reaction balances
  and the world total never changes. The tests check both.
- **Exact energy accounting: yes, by construction.** Light is the only input,
  and the two balances above hold to the unit. The tests check them. The sources
  say only "energy and matter conservation"; the ledger is Chemart's.
- **No light, no ecosystem: yes.** With `light_input=0` each founder can fix
  once from its starting energy, nothing reproduces, and all eight are dead by
  step 2. The tests check that nothing reproduces and nobody survives.
- **Trophic structure: yes, by construction.** Organic matter comes only from
  fixing and death, and goes back to mineral only by scavenging.
- **Selection on the lock-and-key match: yes.** The scrambled founder leaves
  no descendants, and the efficiency rises well above its starting average.
  The tests check this for three seeds.
- **Mutation–selection balance, not improvement: yes, as a negative result.**
  From near-perfect founders the enzyme match decays at the default mutation
  rate (fix efficiency 0.98 → 0.75 in the default run) and is held near its
  start at 0.02. The tests check this. Over 600 steps the decay goes much
  further: fix efficiency ends between 0.22 and 0.65 in seeds 0 to 6, and
  organisms of guild `none` make up 16 to 69 of the survivors.
- **Resistance to phenotypic change (EVOLVE II): no.** The mean sensitivity
  falls below the founders' 3 in every default run (to 2.80 on average over
  seeds 1 to 6), and the tests check that fall. But runs started at other
  values show it is mostly mutation bias. A mutated sensitivity letter is drawn
  at random, with mean 2.5, so the population drifts toward about 2.5 from any
  start. Founders at 1 rise to 2.35 and founders at 4 fall to 2.62 in 600
  steps. The tests start only at 3, so they cannot tell these apart. Making the
  light vary (`light_period` 20 or 100) barely changes the result (2.79 and
  2.75 at step 120).

### What Chemart does not reproduce

- The EVOLVE III competition experiment (constant- versus variable-raised
  populations). A variable light and genome transplants are available, but the
  abstract gives no parameters and does not say what "stage of development"
  means.
- The one-dimensional 1969/1970 model and its state-transition organisms.
- Parasitism and symbiosis (O'Callaghan and Conrad, 1992) and EVOLVE IV's
  metabolite exchange and environment modification.
- Any published number. The table, the costs, the world size and the
  inventory are all Chemart's choices, recorded under *Implementation
  decisions*.

## Further reading

- Vostinar, A. E., Skocelas, K. G., Lalejini, A. & Zaman, L. (2021). Symbiosis
  in digital evolution: past, present, and future. *Frontiers in Ecology and
  Evolution* 9, 739047. <https://doi.org/10.3389/fevo.2021.739047>. Section
  5.1 places EVOLVE in the history of digital evolution.
- O'Callaghan, J. & Conrad, M. (1992). Symbiotic interactions in the EVOLVE
  III ecosystem model. *BioSystems* 26(4), 199–209.
  <https://doi.org/10.1016/0303-2647(92)90025-T>
