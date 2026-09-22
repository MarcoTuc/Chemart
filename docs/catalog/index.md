# Catalog

**72** chemistries, one page each, generated from the catalog so these pages cannot drift from the library. A further 28 archived entries are listed [at the end](#archive).

| | |
|---|---|
| given — a reaction network written down; you choose its rates and initial state | 19 |
| generator — an algorithm computes the network from the chemistry's arguments | 24 |
| gas — a Turing gas: structured molecules react by a procedure and the soup evolves | 29 |
| constructive (open, growing species set) | 43 |
| carry their own rate constants or rate law | 29 |
| carry energetics or thermodynamic consistency | 6 |
| declare a conservation law | 27 |
| define space | 13 |
| define compartments | 7 |

By kind: 57 generator, 7 formalism, 4 wet, 2 analysis, 2 framework.

The catalog is grouped by what each chemistry is: a **given** network, a **generator** of networks, or a Turing **gas** (the reasoning for every entry is in `catalog/TYPES.md`). Columns: **grows** is whether the species set is open and expands at run time; **S/R** is the species and reaction count at *default* parameters with `seed=1`, which for most chemistries scales up considerably; **beyond topology** lists only the capabilities that distinguish entries, since every one of them supplies topology and stoichiometry. **fidelity** says how close the implementation is to a published specification — see [Fidelity and trust](../trust.md).

Each chemistry's own page has the full capability list, the complete attribution, its parameters and its provenance.

## given

19 chemistries: a reaction network written down, possibly as a menu of named variants or from user-supplied rules; you choose its rates and initial state.

| chemistry | origin | family | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Analog computation of algebraic functions with concentrations](analog-function-crn.md) | Hjelmfelt et al. | application | · | book+decisions | 2/2 | kinetics |
| [ARMS (Abstract Rewriting System on Multisets)](arms.md) | Suzuki & Tanaka, 1997 | rewriting | · | reconstructed | 7/6 | — |
| [Varela-Maturana-Uribe autopoiesis model](autopoiesis-vmu.md) | Varela, Maturana & Uribe… | origin-of-life | · | reconstructed | 8/33 | space |
| [Brusselator](brusselator.md) | Prigogine & Lefever, 1968 | application | · | book | 6/4 | kinetics |
| [Colored chameleon chemistry](chameleon.md) | Winkler, 2009 puzzle (Com… | core | · | book | 3/3 | kinetics, conservation |
| [Chemoton](chemoton.md) | Ganti, ~1952/1971 | origin-of-life | · | reconstructed | 38/40 | kinetics, compartments |
| [Conrad's lock-and-key enzymatic / self-assembly processor](conrad-enzymatic.md) | Conrad, 1985-1992 | bio-inspired | · | reconstructed | 12/12 | conservation |
| [Chemical disperser (load balancing)](disperser.md) | Meyer & Tschudin, 2009 | application | · | book+decisions | 12/8 | kinetics, compartments |
| [Flow artificial chemistry](flow-ac.md) | Kreyssig & Dittrich, 2011 | spatial | · | reconstructed | 4/8 | kinetics, space |
| [Hill kinetics (cooperative binding)](hill-kinetics.md) | — | systems-biology | · | book+decisions | 4/3 | kinetics, rate law |
| [Mechanical self-assembly (Hosokawa)](mechanical-self-assembly.md) | Hosokawa, Shimoyama & Miu… | non-chemical | · | book+decisions | 6/9 | kinetics, conservation |
| [Michaelis-Menten enzyme kinetics](michaelis-menten.md) | — | systems-biology | · | book+decisions | 4/3 | kinetics, rate law |
| [Nuclear reaction networks](nuclear-reaction-networks.md) | — | non-chemical | · | book+decisions | 11/6 | conservation |
| [Okamoto's biochemical switch](okamoto-switch.md) | Okamoto, Sakai & Hayashi… | application | · | book+decisions | 8/7 | kinetics |
| [Ono & Ikegami autopoietic protocells](ono-ikegami-protocell.md) | Ono & Ikegami, 1999-2003 | origin-of-life | · | reconstructed | 5/6 | energies, space |
| [Oregonator (Belousov-Zhabotinsky)](oregonator.md) | Field & Noyes, 1974 | wet | · | book+decisions | 7/5 | kinetics, space |
| [P systems (membrane computing)](p-systems.md) | Paun, 1998 | rewriting | · | reconstructed | 8/4 | compartments |
| [Repressilator](repressilator.md) | Elowitz & Leibler, 2000 | wet | · | book | 12/18 | kinetics |
| [Self-propelled oil droplets](self-propelled-droplets.md) | Hanczyc et al., 2007-2011 | wet | · | reconstructed | 6/2 | conservation, space |

## generator

24 chemistries: an algorithm computes the network from its arguments (random draws or a closure); once built, it is treated as a given network.

| chemistry | origin | family | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Algorithmic Chemistry GP (ACGP)](acgp.md) | Banzhaf & Lasarczyk, 2004… | application | · | book+decisions | 12/100 | — |
| [Bagley & Farmer autocatalytic metabolism](bagley-farmer.md) | Bagley & Farmer, 1992 | origin-of-life | yes | reconstructed | 125/696 | kinetics, thermo, conservation, flow |
| [Conservative random chemical reaction networks](bigan-conservative-crn.md) | Bigan, Steyaert & Douady… | network | · | reconstructed | 10/46 | kinetics, rate law, energies, thermo, conservation, flow |
| [Bondable Cellular Automata (BCA)](bondable-ca.md) | Hatcher, Banzhaf & Yu, 20… | automata | yes | book+decisions | 23/11 | conservation |
| [Brane calculi](brane-calculi.md) | Cardelli, 2004 | rewriting | yes | reconstructed | 4/3 | compartments |
| [Chemical Abstract Machine (CHAM)](cham.md) | Berry & Boudol, 1990 | rewriting | yes | reconstructed | 18/30 | compartments |
| [DNA automaton (Benenson-Shapiro)](dna-automaton.md) | Shapiro 1998-1999 (bluepr… | wet | · | reconstructed | 40/44 | kinetics, conservation |
| [Adleman's DNA Hamiltonian Path computation](dna-hpp.md) | Adleman, 1994 | wet | yes | reconstructed | 797/769 | conservation |
| [Bitstring immune system model (idiotypic network)](farmer-immune.md) | Farmer, Packard & Perelso… | bio-inspired | yes | book+decisions | 12/88 | kinetics |
| [Gamma / gamma-calculus](gamma.md) | Banatre & Le Metayer, 198… | rewriting | yes | reconstructed | 4/10 | — |
| [GARD (Graded Autocatalysis Replication Domain)](gard.md) | Segre, Lancet et al., 199… | origin-of-life | · | reconstructed | 200/10200 | kinetics, compartments |
| [Jain-Krishna autocatalytic set model](jain-krishna.md) | Jain & Krishna, 1998-2002 | evolutionary-dynamics | yes | book+decisions | 100/22 | kinetics, flow |
| [Kappa calculus](kappa-calculus.md) | Danos & Laneve, 2004 | rewriting | yes | reconstructed | 10/11 | kinetics, conservation |
| [Kauffman autocatalytic sets (binary polymer model)](kauffman-autocatalytic-sets.md) | Kauffman, 1986 | origin-of-life | yes | book+decisions | 62/1608 | conservation, flow |
| [L-systems](l-systems.md) | Lindenmayer, 1968 | rewriting | yes | reconstructed | 6/5 | — |
| [Metabolic / artificial biochemical network robot controllers](metabolic-robot-controller.md) | Ziegler & Banzhaf, 2001 | application | · | reconstructed | 8/11 | kinetics, conservation, flow |
| [MGS](mgs.md) | Giavitto & Michel, 2001 | rewriting | yes | reconstructed | 21/25 | conservation, space |
| [RAF sets (reflexively autocatalytic, F-generated)](raf.md) | Hordijk & Steel, 2004-2015 | origin-of-life | · | reconstructed | 254/746 | conservation |
| [Random catalytic reaction networks](random-catalytic-networks.md) | Stadler, Fontana & Miller… | evolutionary-dynamics | · | book+decisions | 10/42 | kinetics, flow |
| [Random Boolean Networks (RBN) and RBN World](rbn.md) | Kauffman, 1969 | evolutionary-dynamics | · | reconstructed | 20/40 | conservation |
| [String Metabolic Network (SMN)](smn.md) | Ono, Fujiwara & Yuta, 2005 | systems-biology | yes | book+decisions | 21/24 | conservation, flow |
| [Self-Organizing Assembly Systems (SOAS)](soas.md) | Frei, Di Marzo Serugendo… | non-chemical | yes | reconstructed | 582/1326 | — |
| [Synthon artificial chemistry](synthon.md) | Lenaerts & Bersini, 2009 | systems-biology | yes | reconstructed | 31/161 | kinetics, conservation |
| [ToyChem (graph-based toy model of chemistry)](toychem.md) | Benko, Flamm & Stadler, 2… | network | yes | reconstructed | 11/16 | kinetics, rate law, energies, conservation |

## gas

29 chemistries: a Turing gas: molecules carry structure and a procedure makes them react, so the soup and its network evolve in chemical-evolutionary time.

| chemistry | origin | family | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [AlChemy (lambda-calculus chemistry)](alchemy.md) | Fontana, 1991 | rewriting | yes | reconstructed | 365/1156 | kinetics, flow |
| [Automata reaction (32-bit binary string chemistry)](automata-reaction.md) | Dittrich & Banzhaf, 1998 | automata | yes | reconstructed | 4818/9872 | kinetics, flow |
| [BFF (self-modifying Brainfuck soup)](bff.md) | Agüera y Arcas, Alakuijal… | automata | yes | reconstructed | 837/646 | space |
| [Embedded particles in cellular automata](ca-embedded-particles.md) | Hanson & Crutchfield 1992 | automata | · | reconstructed | 9/7 | space |
| [Chemical Casting Model (CCM)](ccm.md) | Kanada, 1992-1996 | application | · | reconstructed | 64/187 | — |
| [Combinator chemistry](combinator-chemistry.md) | Speroni di Fenizio, 2000 | rewriting | yes | reconstructed | 107/3268 | conservation, flow |
| [Combinatory Chemistry](combinatory-chemistry.md) | Kruszewski & Mikolov, 2020 | rewriting | yes | reconstructed | 591/1001 | conservation |
| [Dorin & Korb virtual ecosystem chemistry](dorin-korb-ecosystem.md) | Dorin & Korb, 2007 | spatial | yes | reconstructed | 26/46 | energies, conservation, space |
| [Fraglets](fraglets.md) | Tschudin, 2003 | application | yes | reconstructed | 7/4 | compartments |
| [High-order chemistry (rules as molecules)](high-order-chem.md) | Yamamoto, 2014 (PyCellChe… | core | yes | book+decisions | 198/154 | — |
| [Machine-tape chemistry](ikegami-hashimoto.md) | Ikegami & Hashimoto, 1995 | automata | yes | reconstructed | 164/1482 | kinetics, flow |
| [Laing's artificial molecular machines](laing-molecular-machines.md) | Laing, 1972-1977 | automata | yes | book+decisions | 17/15 | flow |
| [Matrix chemistry](matrix-chemistry.md) | Banzhaf, 1993 | core | yes | book+decisions | 23/375 | kinetics, flow |
| [Polymers as Turing machines / pattern processing chemistry](mccaskill-polymer-tm.md) | McCaskill, 1988 | automata | yes | reconstructed | 7/12 | flow |
| [Molecular Classifier System (MCS.bl)](mcs-bl.md) | Decraene, Mitchell & McMu… | bio-inspired | yes | reconstructed | 4/6 | kinetics, flow |
| [Molecular Traveling Salesman](molecular-tsp.md) | Banzhaf, 1990 | application | yes | reconstructed | 84/79 | — |
| [Algorithmic chemistry for music composition](music-ac.md) | Miura & Tominaga, 2006 (h… | application | yes | reconstructed | 772/719 | conservation |
| [Network Artificial Chemistry (NAC)](nac.md) | Suzuki, 2004-2009 | network | yes | reconstructed | 61/64 | conservation, space |
| [Prime number (number-division) chemistry](prime-number-chemistry.md) | Banzhaf, Dittrich & Rauhe… | core | yes | book+decisions | 201/159 | — |
| [Artificial chemistry as a proof search system](proof-ac.md) | Busch & Banzhaf, 2003 | application | yes | reconstructed | 100/402 | — |
| [Reflexive artificial chemistry](reflexive-ac.md) | Salzberg, 2007 | rewriting | yes | reconstructed | 6/7 | flow |
| [RNA-folding ribozyme artificial chemistry](rna-folding-ac.md) | Ullrich & Flamm, 2008 | systems-biology | yes | reconstructed | 21/8 | energies, conservation |
| [SAC (string-based artificial chemistry with cells)](sac.md) | Suzuki & Ono, 2002-2003 | bio-inspired | yes | reconstructed | 410/410 | — |
| [Squirm3](squirm3.md) | Hutton, 2002-2007 | spatial | yes | reconstructed | 74/67 | conservation, space |
| [Self-replicating loops in cellular automata](sr-loops.md) | von Neumann 1966 | automata | yes | reconstructed | 1/1 | space |
| [SRSim (rule-based spatial simulator)](srsim.md) | Gruenert & Dittrich, 2010… | systems-biology | yes | reconstructed | 13/19 | kinetics, conservation, space |
| [Stringmol](stringmol.md) | Hickinbotham, Clark, Step… | bio-inspired | yes | reconstructed | 1/1 | flow |
| [Tominaga's stackable-string chemistry](tominaga-stacked-strings.md) | Tominaga et al., 2004-2009 | systems-biology | yes | reconstructed | 47/59 | — |
| [Typogenetics](typogenetics.md) | Hofstadter, 1979 | automata | yes | reconstructed | 4/4 | flow |

## Archive

These entries are not part of the chemistry catalog. `list_chemistries()` and the LLM tools leave them out, but each keeps its specification, generator, tests and page, and `generate_network(id)` still runs it.

### artificial-life

5 entries: artificial life rather than artificial chemistry.

| chemistry | origin | family | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Avida](avida.md) | Adami & Brown, 1994 | automata | yes | reconstructed | 178/248 | space |
| [Core War / Redcode](corewar.md) | Dewdney, 1984 | automata | yes | reconstructed | 192/771 | conservation, space |
| [Coreworld (VENUS)](coreworld.md) | Rasmussen, Knudsen, Feldb… | automata | yes | reconstructed | 10/740 | conservation, space |
| [Swarm Chemistry](swarm-chemistry.md) | Sayama, 2009-2011 | spatial | · | reconstructed | 2/0 | rate law, space |
| [Tierra](tierra.md) | Ray, 1991 | automata | yes | reconstructed | 63/131 | — |

### pruned

23 entries: set aside from the catalog.

| chemistry | origin | family | grows | fidelity | S/R | beyond topology |
|---|---|---|:--:|---|--:|---|
| [Aevol](aevol.md) | Knibbe, Beslon et al., 20… | systems-biology | yes | reconstructed | 134/377 | space |
| [Artificial Regulatory Network (ARN)](arn.md) | Banzhaf, 2003 | network | · | reconstructed | 7/98 | kinetics, flow |
| [BNC (bond-number chemistry) cell model](bnc-cell.md) | Hintze & Adami, 2008 | systems-biology | · | reconstructed | 87/11814 | kinetics, conservation, space, compartments |
| [Cellular Potts + GRN evo-devo models](cpm-grn-evodevo.md) | Hogeweg, 2000 | systems-biology | · | reconstructed | 10/17 | energies, space, compartments |
| [Reversible dimerization](dimerization.md) | textbook | core | · | book | 3/2 | kinetics |
| [Ecolab](ecolab.md) | Standish, 1994-2004 | evolutionary-dynamics | yes | reconstructed | 54/686 | kinetics |
| [Arrhenius-gated collision algorithm](energy-gated-collision.md) | Banzhaf & Yamamoto, 2015… | systems-biology | yes | book+decisions | 4/2 | kinetics, rate law, energies, thermo |
| [EVOLVE virtual ecosystems](evolve-series.md) | Conrad & Pattee, 1969-1970 | evolutionary-dynamics | yes | book+decisions | 180/450 | energies, conservation, space |
| [French flag model (positional information)](french-flag.md) | Wolpert, 1969 | systems-biology | · | book+decisions | 60/30 | space |
| [HBCB / PSD degradation-and-reuse chemistry](hbcb-psd.md) | Oohashi, Ueno, Maekawa, K… | systems-biology | · | book+decisions | 16/24 | energies, conservation |
| [Isologous diversification](isologous-diversification.md) | Kaneko & Yomo, 1994-1999 | systems-biology | · | reconstructed | 10/36 | kinetics, compartments |
| [Logistic growth ('replicate and fight')](logistic-chemistry.md) | — | evolutionary-dynamics | · | book+decisions | 1/2 | kinetics |
| [Lotka-Volterra](lotka-volterra.md) | — | evolutionary-dynamics | · | book+decisions | 2/3 | kinetics |
| [N-economy (natural number economy)](n-economy.md) | Straatman, Banzhaf et al. | non-chemical | yes | book+decisions | 7/5 | — |
| [Naming game as an artificial chemistry](naming-game-ac.md) | De Beule, Hovig & Benson… | application | · | book+decisions | 7/21 | kinetics |
| [Kauffman NK model](nk-landscape.md) | Kauffman, 1993 | evolutionary-dynamics | · | book+decisions | 64/448 | kinetics, flow |
| [Computing with chemical organizations](organization-computing.md) | Matsumaru, Speroni di Fen… | application | · | book+decisions | 6/7 | — |
| [Quasispecies equation](quasispecies.md) | Eigen, 1971 | evolutionary-dynamics | · | book+decisions | 64/4032 | kinetics, flow |
| [Replication and death](replication-death.md) | textbook evolutionary dyn… | evolutionary-dynamics | · | book+decisions | 1/2 | kinetics |
| [Replicator equation (evolutionary game dynamics)](replicator-equation.md) | — | evolutionary-dynamics | · | book+decisions | 3/6 | kinetics, flow |
| [Selection equation under a dilution flow](selection-equation.md) | — | evolutionary-dynamics | · | book+decisions | 3/3 | kinetics, rate law, flow |
| [Artificial chemistry of social communication](social-communication-ac.md) | Dittrich, Kron & Banzhaf… | non-chemical | · | reconstructed | 10/57 | — |
| [Urdar](urdar.md) | Gerlee & Lundh, 2010 | evolutionary-dynamics | yes | reconstructed | 2847/2453 | energies |
