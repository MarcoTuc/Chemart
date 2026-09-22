# What each chemistry is: given, generator or gas

Every catalog entry records a `type` (see `catalog/SCHEMA.md`). This table is the reasoning behind each value.

- **given**: the chemistry *is* a reaction network, written down. It may be a menu of named variants (Michaelis-Menten's `form`, flow-ac's `example`) or a rule set the user supplies (P systems, ARMS). You choose its rates and initial state, then simulate it.
- **generator**: an algorithm computes the network from the chemistry's arguments, by random draws (Kauffman sets, random catalytic networks) or by closing a set of rules over molecules (gamma, CHAM, ToyChem). Once built, the network is treated as a given one; its measures change with the arguments.
- **gas**: a Turing gas. Molecules carry structure (terms, programs, strings) and a procedure makes them react, so the soup and its network evolve in chemical-evolutionary time (AlChemy, BFF, Combinatory Chemistry).

How a chemistry can be run is a separate question, answered by its module: `generate` returns one network, `evolve` runs a process frame by frame. A gas usually has `evolve`, and some also have `generate` (the closure of their rules); a given network may have `evolve` when the chemistry defines a reactor for it, such as a lattice.

*Borderline* marks entries where another type is defensible; the note says why this one was chosen.

| chemistry | type | why |
|---|---|---|
| `acgp` | generator | A random program is drawn and read as a network: each instruction r3 <- r1 op r2 becomes r1 + r2 -> r3. |
| `aevol` (archived: pruned) | gas | Genomes evolve on a lattice; the network records the expressions, births and replacements of one run. |
| `alchemy` | gas | Applying one lambda-term to another computes the product; the soup of terms evolves (its closure is a second face). |
| `analog-function-crn` | given | A menu of designed networks (X -> X + Y, 2Y -> 0, ...) whose steady states compute functions; `function` picks one. |
| `arms` | given *(borderline)* | The chemistry is its rule list over atomic symbols ('3 a -> c'), named or supplied; drawing a random subset of rules would make it a generator. |
| `arn` (archived: pruned) | generator *(borderline)* | A random genome is decoded into a regulatory network of proteins; the ODE then runs on it. |
| `automata-reaction` | gas | A 32-bit word acts as an automaton on another to compute the product; the soup evolves. |
| `autopoiesis-vmu` | given *(borderline)* | Three reaction types over a fixed set of species (substrate, catalyst, link states); the lattice is a spatial reactor for them. |
| `avida` (archived: artificial-life) | gas | Self-replicating programs run on a grid. |
| `bagley-farmer` | generator *(borderline)* | Every polymer up to max_length with all splits and random catalytic links; the optional threshold prunes the network by its dynamics. |
| `bff` | gas | Pairs of self-modifying programs run as one tape and rewrite each other; the soup evolves towards self-replicators. |
| `bigan-conservative-crn` | generator | A random network drawn from stated distributions: species, reactions and free energies. |
| `bnc-cell` (archived: pruned) | generator | All molecules allowed by the bond rule up to max_length, and every recombination between them. |
| `bondable-ca` | generator | Whether atoms bond is computed from their CA dynamics; the closure over atoms generates molecules and bonds. |
| `brane-calculi` | generator | Reduction rules rewrite nested membrane terms; the closure generates the configurations. |
| `brusselator` | given | Four reactions, written down. |
| `ca-embedded-particles` | gas *(borderline)* | The particles and their collisions are read out of a cellular automaton run (its generate face is the published catalogue, a written-down network). |
| `ccm` | gas *(borderline)* | A search procedure: rules rewrite atom states when a local order test passes; the species set is fixed, the trace evolves. |
| `cham` | generator | Heating, cooling and reaction rules rewrite process-algebra terms; the closure generates the network. |
| `chameleon` | given | Three reactions over three colours, written down. |
| `chemoton` | given | Gánti's hand-designed network of metabolism, template and membrane; N only sets the template length. |
| `combinator-chemistry` | gas | Applying one combinator to another and reducing computes the products; the soup evolves (its closure is a second face). |
| `combinatory-chemistry` | gas | Reducing S, K, I expressions computes products and consumes reactants; the soup of expressions evolves from free atoms. |
| `conrad-enzymatic` | given *(borderline)* | The binding network of one wet experiment; shape complementarity fixes it once from the ligand and enzyme set. |
| `corewar` (archived: artificial-life) | gas | Programs execute in a shared core. |
| `coreworld` (archived: artificial-life) | gas | Redcode runs with resources. |
| `cpm-grn-evodevo` (archived: pruned) | gas *(borderline)* | Cell types appear as the GRN and the Cellular Potts dynamics run. |
| `dimerization` (archived: pruned) | given | A + B <-> C, written down. |
| `disperser` | given | One catalysed transfer per link of a supplied computer network. |
| `dna-automaton` | generator *(borderline)* | The processing stages are computed by running the automaton on the inputs; the wet protocol itself is fixed. |
| `dna-hpp` | generator | Hybridisation and ligation build ever longer strands; the closure generates them. |
| `dorin-korb-ecosystem` | gas | Atoms bond and unbond on a grid; molecules are the clusters a run produces. |
| `ecolab` (archived: pruned) | generator *(borderline)* | A Lotka-Volterra ecosystem whose species set grows by mutation during the run. |
| `energy-gated-collision` (archived: pruned) | given | A reactor algorithm applied to an explicit small set of species and reactions. |
| `evolve-series` (archived: pruned) | gas | Genomes evolve in an ecosystem. |
| `farmer-immune` | generator | Random antibody bitstrings whose matchings fix every reaction and rate. |
| `flow-ac` | given | The book's fixed networks placed in a flow; `example` picks one, space is the reactor. |
| `fraglets` | gas | Fraglets are programs that transform fraglets; running them evolves the soup. |
| `french-flag` (archived: pruned) | given | Each cell's fate reaction is fixed by the morphogen threshold at its position. |
| `gamma` | generator | Conditional rules rewrite the values in the multiset; the closure generates the network. |
| `gard` | generator | N_G lipid types with a catalytic matrix drawn at random. |
| `hbcb-psd` (archived: pruned) | generator | A hierarchy of classes and the reactions between neighbouring levels, built by closure from the parameters. |
| `high-order-chem` | gas | Rule molecules act on data molecules to produce new data; the soup evolves. |
| `hill-kinetics` | given | Cooperative binding written as fixed reactions; `form` picks the variant. |
| `ikegami-hashimoto` | gas | Machines read tapes and write new machines and tapes; the soup evolves. |
| `isologous-diversification` (archived: pruned) | generator | A random catalytic network inside every cell; cells divide and differentiate. |
| `jain-krishna` | generator *(borderline)* | A random catalytic graph; its rewiring over graph updates is an evolution of the network by an external rule. |
| `kappa-calculus` | generator | Rules over site-graph patterns; the closure generates complexes and reactions. |
| `kauffman-autocatalytic-sets` | generator | All polymers up to a length, all condensations and cleavages, catalysis drawn at random. |
| `l-systems` | generator | Productions rewrite words in parallel; the derivation closure generates new words. |
| `laing-molecular-machines` | gas | Machines run on tapes and produce new strings. |
| `logistic-chemistry` (archived: pruned) | given | Two reactions, written down. |
| `lotka-volterra` (archived: pruned) | given | Three reactions, written down. |
| `matrix-chemistry` | gas | A string folded into a matrix acts on another string to compute a third; the soup evolves (its closure is a second face). |
| `mccaskill-polymer-tm` | gas | Polymers process polymers as tapes and write new strings. |
| `mcs-bl` | gas | Strings act on strings by pattern matching and rewriting, producing new strings. |
| `mechanical-self-assembly` | given | Nine assembly reactions between monomers and clusters, written down. |
| `metabolic-robot-controller` | generator | A random mass-balanced reaction graph used as a controller. |
| `mgs` | generator | Rules rewrite paths in topological collections; the closure generates new collections. |
| `michaelis-menten` | given | E + S <-> ES -> E + P, written down; `form` picks the variant. |
| `molecular-tsp` | gas | Machines rewrite candidate tours into new tours; the soup of tours evolves. |
| `music-ac` | gas | Recombination rules join melodic fragments into new objects; the soup evolves. |
| `n-economy` (archived: pruned) | given | Production reactions with integer stoichiometry, listed. |
| `nac` | gas | Local rewiring of a graph of molecules; the clusters are what a run produces. |
| `naming-game-ac` (archived: pruned) | given | Fixed reactions between a meaning, words and codes. |
| `nk-landscape` (archived: pruned) | generator | Every genotype of length N with random fitness contributions: replication and mutation reactions built by formula. |
| `nuclear-reaction-networks` | given | Measured nuclear reactions such as the pp chain and the CNO cycle, written down (without rates). |
| `okamoto-switch` | given | A fixed flip-flop of seven reactions. |
| `ono-ikegami-protocell` | given *(borderline)* | Six reactions among five particle types; the lattice and its membrane physics are a spatial reactor for them. |
| `oregonator` | given | Five reactions, written down. |
| `organization-computing` (archived: pruned) | given | Designed gate reactions (XOR, maximal independent set) over presence/absence species. |
| `p-systems` | given | A named or supplied rule set over atomic symbols in regions: the rules, placed in their membranes, are the network. |
| `prime-number-chemistry` | gas | Division of integers produces new integers; the soup evolves. |
| `proof-ac` | gas | Resolution of clauses produces new clauses; the soup evolves (its closure is a second face). |
| `quasispecies` (archived: pruned) | generator | Every genotype of length L, with replication and mutation reactions built by formula. |
| `raf` | generator | Kauffman's binary polymer model up to length n, with catalysis drawn at random, to analyse for RAF sets. |
| `random-catalytic-networks` | generator | A random catalytic network over structureless species. |
| `rbn` | generator *(borderline)* | A random Boolean network drawn once, each truth-table row one reaction (RBN World, built from the same networks, is the gas rbn-world). |
| `rbn-world` | gas | Atoms are random Boolean networks with bonding sites; whether two bond is computed from their attractors, and the soup of molecules evolves. |
| `reflexive-ac` | gas | Finite state machines compose with machines into new machines. |
| `replication-death` (archived: pruned) | given | Birth and death, written down. |
| `replicator-equation` (archived: pruned) | given | Reactions fixed by a supplied payoff matrix. |
| `repressilator` | given | The three-gene repression ring, written down. |
| `rna-folding-ac` | gas | Each sequence is folded and its structure read as the reaction it catalyses; the products are computed. |
| `sac` | gas | Strings rewrite and cut strings, producing new strings. |
| `selection-equation` (archived: pruned) | given | Replication and dilution reactions, written down. |
| `self-propelled-droplets` | given | The hydrolysis reactions of one wet experiment, written down. |
| `smn` | generator | A random enzyme genome; its compound set is closed under the enzymes' ligation, cleavage and recombination. |
| `soas` | generator | CHAM rules combine resources and skills; the closure generates coalitions. |
| `social-communication-ac` (archived: pruned) | gas *(borderline)* | Agents choose each reply from memories that change during the run. |
| `squirm3` | gas | Atoms bond and change state on a grid, producing molecules such as replicated templates. |
| `sr-loops` | gas *(borderline)* | Self-replicating loops emerge in a cellular automaton; the network is read from the run. |
| `srsim` | gas *(borderline)* | Rules bind and change spatial molecules; complexes are what a run forms (only the spatial run is implemented). |
| `stringmol` | gas | Strings bind and execute on each other, producing new strings. |
| `swarm-chemistry` (archived: artificial-life) | given *(borderline)* | No reactions: a fixed set of particle recipes and an interaction law. |
| `synthon` | generator *(borderline)* | Reaction classes rewrite molecule graphs; the closure, or a kinetic sample of it, generates the network from the starting molecules. |
| `tierra` (archived: artificial-life) | gas | Machine-code creatures replicate and mutate. |
| `tominaga-stacked-strings` | gas | Recombination rules on stacked strings produce new objects. |
| `toychem` | generator | Graph rewriting rules on structural formulas generate the molecules and reactions. |
| `typogenetics` | gas | Strands code enzymes that act on strands, producing new strands. |
| `urdar` (archived: pruned) | gas | CA-rule organisms rewrite bitstring metabolites and reproduce. |
