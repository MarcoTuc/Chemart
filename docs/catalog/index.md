# Catalog

**72** chemistries, one page each, generated from the catalog so these pages cannot drift from the library. A further 28 archived entries are listed [at the end](#archive).

| | |
|---|---|
| network given (written down, built by a formula, or drawn at random) | 29 |
| network generated (the output of the chemistry's algorithm) | 43 |
| constructive (open, growing species set) | 43 |
| carry their own rate constants or rate law | 29 |
| carry energetics or thermodynamic consistency | 6 |
| declare a conservation law | 27 |
| define space | 13 |
| define compartments | 7 |

By kind: 57 generator, 7 formalism, 4 wet, 2 analysis, 2 framework.

Columns: **network** is *given* when the chemistry is a reaction network Chemart instantiates from its parameters, and *generated* when the network is the output of the chemistry's algorithm (see `catalog/NETWORKS.md`); **grows** is whether the species set is open and expands at run time; **S/R** is the species and reaction count at *default* parameters with `seed=1`, which for most chemistries scales up considerably; **beyond topology** lists only the capabilities that distinguish entries, since every one of them supplies topology and stoichiometry. **fidelity** says how close the implementation is to a published specification — see [Fidelity and trust](../trust.md).

Each chemistry's own page has the full capability list, the complete attribution, its parameters and its provenance.

## application

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Algorithmic Chemistry GP (ACGP)](acgp.md) | Banzhaf & Lasarczyk, 2004… | given | · | book+decisions | 12/100 | — |
| [Analog computation of algebraic functions with concentrations](analog-function-crn.md) | Hjelmfelt et al. | given | · | book+decisions | 2/2 | kinetics |
| [Brusselator](brusselator.md) | Prigogine & Lefever, 1968 | given | · | book | 6/4 | kinetics |
| [Chemical Casting Model (CCM)](ccm.md) | Kanada, 1992-1996 | generated | · | reconstructed | 64/187 | — |
| [Chemical disperser (load balancing)](disperser.md) | Meyer & Tschudin, 2009 | given | · | book+decisions | 12/8 | kinetics, compartments |
| [Fraglets](fraglets.md) | Tschudin, 2003 | generated | yes | reconstructed | 7/4 | compartments |
| [Metabolic / artificial biochemical network robot controllers](metabolic-robot-controller.md) | Ziegler & Banzhaf, 2001 | given | · | reconstructed | 8/11 | kinetics, conservation, flow |
| [Molecular Traveling Salesman](molecular-tsp.md) | Banzhaf, 1990 | generated | yes | reconstructed | 84/79 | — |
| [Algorithmic chemistry for music composition](music-ac.md) | Miura & Tominaga, 2006 (h… | generated | yes | reconstructed | 772/719 | conservation |
| [Okamoto's biochemical switch](okamoto-switch.md) | Okamoto, Sakai & Hayashi… | given | · | book+decisions | 8/7 | kinetics |
| [Artificial chemistry as a proof search system](proof-ac.md) | Busch & Banzhaf, 2003 | generated | yes | reconstructed | 100/402 | — |

## automata

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Automata reaction (32-bit binary string chemistry)](automata-reaction.md) | Dittrich & Banzhaf, 1998 | generated | yes | reconstructed | 4818/9872 | kinetics, flow |
| [BFF (self-modifying Brainfuck soup)](bff.md) | Agüera y Arcas, Alakuijal… | generated | yes | reconstructed | 837/646 | space |
| [Bondable Cellular Automata (BCA)](bondable-ca.md) | Hatcher, Banzhaf & Yu, 20… | generated | yes | book+decisions | 23/11 | conservation |
| [Embedded particles in cellular automata](ca-embedded-particles.md) | Hanson & Crutchfield 1992 | generated | · | reconstructed | 9/7 | space |
| [Machine-tape chemistry](ikegami-hashimoto.md) | Ikegami & Hashimoto, 1995 | generated | yes | reconstructed | 164/1482 | kinetics, flow |
| [Laing's artificial molecular machines](laing-molecular-machines.md) | Laing, 1972-1977 | generated | yes | book+decisions | 17/15 | flow |
| [Polymers as Turing machines / pattern processing chemistry](mccaskill-polymer-tm.md) | McCaskill, 1988 | generated | yes | reconstructed | 7/12 | flow |
| [Self-replicating loops in cellular automata](sr-loops.md) | von Neumann 1966 | generated | yes | reconstructed | 1/1 | space |
| [Typogenetics](typogenetics.md) | Hofstadter, 1979 | generated | yes | reconstructed | 4/4 | flow |

## bio-inspired

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Conrad's lock-and-key enzymatic / self-assembly processor](conrad-enzymatic.md) | Conrad, 1985-1992 | given | · | reconstructed | 12/12 | conservation |
| [Bitstring immune system model (idiotypic network)](farmer-immune.md) | Farmer, Packard & Perelso… | given | yes | book+decisions | 12/88 | kinetics |
| [Molecular Classifier System (MCS.bl)](mcs-bl.md) | Decraene, Mitchell & McMu… | generated | yes | reconstructed | 4/6 | kinetics, flow |
| [SAC (string-based artificial chemistry with cells)](sac.md) | Suzuki & Ono, 2002-2003 | generated | yes | reconstructed | 410/410 | — |
| [Stringmol](stringmol.md) | Hickinbotham, Clark, Step… | generated | yes | reconstructed | 1/1 | flow |

## core

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Colored chameleon chemistry](chameleon.md) | Winkler, 2009 puzzle (Com… | given | · | book | 3/3 | kinetics, conservation |
| [High-order chemistry (rules as molecules)](high-order-chem.md) | Yamamoto, 2014 (PyCellChe… | generated | yes | book+decisions | 198/154 | — |
| [Matrix chemistry](matrix-chemistry.md) | Banzhaf, 1993 | generated | yes | book+decisions | 23/375 | kinetics, flow |
| [Prime number (number-division) chemistry](prime-number-chemistry.md) | Banzhaf, Dittrich & Rauhe… | generated | yes | book+decisions | 201/159 | — |

## evolutionary-dynamics

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Jain-Krishna autocatalytic set model](jain-krishna.md) | Jain & Krishna, 1998-2002 | generated | yes | book+decisions | 100/22 | kinetics, flow |
| [Random catalytic reaction networks](random-catalytic-networks.md) | Stadler, Fontana & Miller… | given | · | book+decisions | 10/42 | kinetics, flow |
| [Random Boolean Networks (RBN) and RBN World](rbn.md) | Kauffman, 1969 | given | · | reconstructed | 20/40 | conservation |

## network

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Conservative random chemical reaction networks](bigan-conservative-crn.md) | Bigan, Steyaert & Douady… | given | · | reconstructed | 10/46 | kinetics, rate law, energies, thermo, conservation, flow |
| [Network Artificial Chemistry (NAC)](nac.md) | Suzuki, 2004-2009 | generated | yes | reconstructed | 61/64 | conservation, space |
| [ToyChem (graph-based toy model of chemistry)](toychem.md) | Benko, Flamm & Stadler, 2… | generated | yes | reconstructed | 11/16 | kinetics, rate law, energies, conservation |

## non-chemical

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Mechanical self-assembly (Hosokawa)](mechanical-self-assembly.md) | Hosokawa, Shimoyama & Miu… | given | · | book+decisions | 6/9 | kinetics, conservation |
| [Nuclear reaction networks](nuclear-reaction-networks.md) | — | given | · | book+decisions | 11/6 | conservation |
| [Self-Organizing Assembly Systems (SOAS)](soas.md) | Frei, Di Marzo Serugendo… | generated | yes | reconstructed | 582/1326 | — |

## origin-of-life

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Varela-Maturana-Uribe autopoiesis model](autopoiesis-vmu.md) | Varela, Maturana & Uribe… | given | · | reconstructed | 8/33 | space |
| [Bagley & Farmer autocatalytic metabolism](bagley-farmer.md) | Bagley & Farmer, 1992 | given | yes | reconstructed | 125/696 | kinetics, thermo, conservation, flow |
| [Chemoton](chemoton.md) | Ganti, ~1952/1971 | given | · | reconstructed | 38/40 | kinetics, compartments |
| [GARD (Graded Autocatalysis Replication Domain)](gard.md) | Segre, Lancet et al., 199… | given | · | reconstructed | 200/10200 | kinetics, compartments |
| [Kauffman autocatalytic sets (binary polymer model)](kauffman-autocatalytic-sets.md) | Kauffman, 1986 | given | yes | book+decisions | 62/1608 | conservation, flow |
| [Ono & Ikegami autopoietic protocells](ono-ikegami-protocell.md) | Ono & Ikegami, 1999-2003 | given | · | reconstructed | 5/6 | energies, space |
| [RAF sets (reflexively autocatalytic, F-generated)](raf.md) | Hordijk & Steel, 2004-2015 | given | · | reconstructed | 254/746 | conservation |

## rewriting

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [AlChemy (lambda-calculus chemistry)](alchemy.md) | Fontana, 1991 | generated | yes | reconstructed | 365/1156 | kinetics, flow |
| [ARMS (Abstract Rewriting System on Multisets)](arms.md) | Suzuki & Tanaka, 1997 | given | · | reconstructed | 7/6 | — |
| [Brane calculi](brane-calculi.md) | Cardelli, 2004 | generated | yes | reconstructed | 4/3 | compartments |
| [Chemical Abstract Machine (CHAM)](cham.md) | Berry & Boudol, 1990 | generated | yes | reconstructed | 18/30 | compartments |
| [Combinator chemistry](combinator-chemistry.md) | Speroni di Fenizio, 2000 | generated | yes | reconstructed | 107/3268 | conservation, flow |
| [Combinatory Chemistry](combinatory-chemistry.md) | Kruszewski & Mikolov, 2020 | generated | yes | reconstructed | 591/1001 | conservation |
| [Gamma / gamma-calculus](gamma.md) | Banatre & Le Metayer, 198… | generated | yes | reconstructed | 4/10 | — |
| [Kappa calculus](kappa-calculus.md) | Danos & Laneve, 2004 | generated | yes | reconstructed | 10/11 | kinetics, conservation |
| [L-systems](l-systems.md) | Lindenmayer, 1968 | generated | yes | reconstructed | 6/5 | — |
| [MGS](mgs.md) | Giavitto & Michel, 2001 | generated | yes | reconstructed | 21/25 | conservation, space |
| [P systems (membrane computing)](p-systems.md) | Paun, 1998 | given | · | reconstructed | 8/4 | compartments |
| [Reflexive artificial chemistry](reflexive-ac.md) | Salzberg, 2007 | generated | yes | reconstructed | 6/7 | flow |

## spatial

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Dorin & Korb virtual ecosystem chemistry](dorin-korb-ecosystem.md) | Dorin & Korb, 2007 | generated | yes | reconstructed | 26/46 | energies, conservation, space |
| [Flow artificial chemistry](flow-ac.md) | Kreyssig & Dittrich, 2011 | given | · | reconstructed | 4/8 | kinetics, space |
| [Squirm3](squirm3.md) | Hutton, 2002-2007 | generated | yes | reconstructed | 74/67 | conservation, space |

## systems-biology

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Hill kinetics (cooperative binding)](hill-kinetics.md) | — | given | · | book+decisions | 4/3 | kinetics, rate law |
| [Michaelis-Menten enzyme kinetics](michaelis-menten.md) | — | given | · | book+decisions | 4/3 | kinetics, rate law |
| [RNA-folding ribozyme artificial chemistry](rna-folding-ac.md) | Ullrich & Flamm, 2008 | generated | yes | reconstructed | 21/8 | energies, conservation |
| [String Metabolic Network (SMN)](smn.md) | Ono, Fujiwara & Yuta, 2005 | generated | yes | book+decisions | 21/24 | conservation, flow |
| [SRSim (rule-based spatial simulator)](srsim.md) | Gruenert & Dittrich, 2010… | generated | yes | reconstructed | 13/19 | kinetics, conservation, space |
| [Synthon artificial chemistry](synthon.md) | Lenaerts & Bersini, 2009 | generated | yes | reconstructed | 31/161 | kinetics, conservation |
| [Tominaga's stackable-string chemistry](tominaga-stacked-strings.md) | Tominaga et al., 2004-2009 | generated | yes | reconstructed | 47/59 | — |

## wet

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [DNA automaton (Benenson-Shapiro)](dna-automaton.md) | Shapiro 1998-1999 (bluepr… | generated | · | reconstructed | 40/44 | kinetics, conservation |
| [Adleman's DNA Hamiltonian Path computation](dna-hpp.md) | Adleman, 1994 | generated | yes | reconstructed | 797/769 | conservation |
| [Oregonator (Belousov-Zhabotinsky)](oregonator.md) | Field & Noyes, 1974 | given | · | book+decisions | 7/5 | kinetics, space |
| [Repressilator](repressilator.md) | Elowitz & Leibler, 2000 | given | · | book | 12/18 | kinetics |
| [Self-propelled oil droplets](self-propelled-droplets.md) | Hanczyc et al., 2007-2011 | given | · | reconstructed | 6/2 | conservation, space |

## Archive

These entries are not part of the chemistry catalog. `list_chemistries()` and the LLM tools leave them out, but each keeps its specification, generator, tests and page, and `generate_network(id)` still runs it.

### artificial-life

5 entries: artificial life rather than artificial chemistry.

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Avida](avida.md) | Adami & Brown, 1994 | generated | yes | reconstructed | 178/248 | space |
| [Core War / Redcode](corewar.md) | Dewdney, 1984 | generated | yes | reconstructed | 192/771 | conservation, space |
| [Coreworld (VENUS)](coreworld.md) | Rasmussen, Knudsen, Feldb… | generated | yes | reconstructed | 10/740 | conservation, space |
| [Swarm Chemistry](swarm-chemistry.md) | Sayama, 2009-2011 | given | · | reconstructed | 2/0 | rate law, space |
| [Tierra](tierra.md) | Ray, 1991 | generated | yes | reconstructed | 63/131 | — |

### pruned

23 entries: set aside from the catalog.

| chemistry | origin | network | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Aevol](aevol.md) | Knibbe, Beslon et al., 20… | generated | yes | reconstructed | 134/377 | space |
| [Artificial Regulatory Network (ARN)](arn.md) | Banzhaf, 2003 | given | · | reconstructed | 7/98 | kinetics, flow |
| [BNC (bond-number chemistry) cell model](bnc-cell.md) | Hintze & Adami, 2008 | given | · | reconstructed | 87/11814 | kinetics, conservation, space, compartments |
| [Cellular Potts + GRN evo-devo models](cpm-grn-evodevo.md) | Hogeweg, 2000 | generated | · | reconstructed | 10/17 | energies, space, compartments |
| [Reversible dimerization](dimerization.md) | textbook | given | · | book | 3/2 | kinetics |
| [Ecolab](ecolab.md) | Standish, 1994-2004 | generated | yes | reconstructed | 54/686 | kinetics |
| [Arrhenius-gated collision algorithm](energy-gated-collision.md) | Banzhaf & Yamamoto, 2015… | given | yes | book+decisions | 4/2 | kinetics, rate law, energies, thermo |
| [EVOLVE virtual ecosystems](evolve-series.md) | Conrad & Pattee, 1969-1970 | generated | yes | book+decisions | 180/450 | energies, conservation, space |
| [French flag model (positional information)](french-flag.md) | Wolpert, 1969 | given | · | book+decisions | 60/30 | space |
| [HBCB / PSD degradation-and-reuse chemistry](hbcb-psd.md) | Oohashi, Ueno, Maekawa, K… | given | · | book+decisions | 16/24 | energies, conservation |
| [Isologous diversification](isologous-diversification.md) | Kaneko & Yomo, 1994-1999 | given | · | reconstructed | 10/36 | kinetics, compartments |
| [Logistic growth ('replicate and fight')](logistic-chemistry.md) | — | given | · | book+decisions | 1/2 | kinetics |
| [Lotka-Volterra](lotka-volterra.md) | — | given | · | book+decisions | 2/3 | kinetics |
| [N-economy (natural number economy)](n-economy.md) | Straatman, Banzhaf et al. | given | yes | book+decisions | 7/5 | — |
| [Naming game as an artificial chemistry](naming-game-ac.md) | De Beule, Hovig & Benson… | given | · | book+decisions | 7/21 | kinetics |
| [Kauffman NK model](nk-landscape.md) | Kauffman, 1993 | given | · | book+decisions | 64/448 | kinetics, flow |
| [Computing with chemical organizations](organization-computing.md) | Matsumaru, Speroni di Fen… | given | · | book+decisions | 6/7 | — |
| [Quasispecies equation](quasispecies.md) | Eigen, 1971 | given | · | book+decisions | 64/4032 | kinetics, flow |
| [Replication and death](replication-death.md) | textbook evolutionary dyn… | given | · | book+decisions | 1/2 | kinetics |
| [Replicator equation (evolutionary game dynamics)](replicator-equation.md) | — | given | · | book+decisions | 3/6 | kinetics, flow |
| [Selection equation under a dilution flow](selection-equation.md) | — | given | · | book+decisions | 3/3 | kinetics, rate law, flow |
| [Artificial chemistry of social communication](social-communication-ac.md) | Dittrich, Kron & Banzhaf… | generated | · | reconstructed | 10/57 | — |
| [Urdar](urdar.md) | Gerlee & Lundh, 2010 | generated | yes | reconstructed | 2847/2453 | energies |
