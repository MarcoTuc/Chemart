# Findings from the explainer writers

While writing each chemistry's documentation page, the writer agents checked the entry against the book, the source papers, the generator, the tests and their own runs. They were allowed to fix only clear factual errors in the YAML prose (those fixes are in the git history). Everything else they noticed is listed here, unverified: bugs, tests that do not check what they claim, catalog text that disagrees with the code or the sources, and errata in the book and papers.

Treat each item as a lead to check, not as an established fact.

## Chemistry catalog

### acgp

- tests/chemistries/test_acgp.py does not exist. The page footer (generator) links to it, but the only acgp test is test_acgp_program_is_the_reaction_multiset in tests/chemistries/test_w2_book.py.
- The reaction drops the operator, so different instructions can collapse into the same reaction. With seed=1 at defaults there are 100 instructions, 99 distinct instructions and 96 distinct reactions, e.g. 'r11 = r5 / r4' and 'r11 = r5 - r4' both become 'r5 + r4 -> r11'. Repeated reactions in the network are therefore not always repeated instructions. The decision text ('repeated instructions stay repeated') is true but can mislead.
- provides lists 'catalysts' only because an instruction that writes back into one of its sources (r10 = r10 * r9) becomes r10 + r9 -> r10. Chemically this is an artefact of the reading, not catalysis.
- The A notes describe start()/next() random execution, but Chemart has no executor. The network carries no dynamics at all.
- The operator set {+,-,*,/} is narrower than the papers' sets: add, sub, div, mult, pow, and, or, not in the 2005 chapter, and AND/OR/NAND/NOR for parity. There are no constant registers either.
- Book §16.6 calls the UCI task 'thyroid cancer'. The source paper classifies thyroid function (normal vs hyper/hypo).
- Book §17.1 credits resilience to instruction deletion to [575], not ACGP. The old YAML phenomenon had given ACGP that credit (now fixed).

### alchemy

- YAML alias 'Turing gas' could not be checked. It is not in the book section, the Fontana & Buss 1994 paper or Mathis et al. 2024; it presumably comes from Fontana 1991, which I did not read. I left it in and did not use it on the page.
- The Level 0 test (tests/chemistries/test_alchemy.py::test_level0_reactor_collapses_to_copiers) checks only seeds 0 and 1. Over seeds 0-19 (M=100, 20,000 collisions, no filter) only 7 of 20 runs ended with a single species (always λx1.x1). The rest ended with 2 to 41 species, some of them a projector family. So 'collapse to a copier' is not the typical outcome with Chemart's default random generator. This is consistent with Mathis et al. 2024 but not with the paper's 'in many instances', and no test checks the statistic.
- The book's 'syntactic convergence of survivors' phenomenon has no test of its own. It is covered only indirectly, by the Level 1 test in which the survivors form a single projector family.
- The paper's ODE result for the A10 family (a single stable fixed point with equal rates, a limit cycle when the basic cycle runs 5:1 faster) is not tested. test_rates_give_the_flow_reactor_equation checks only that the right-hand side matches eq. 18.
- Mathis et al. 2024 (arXiv:2408.12137v2, sec. 4.1) contains a garbled sentence: 'Even replacing 50% of the expressions with the identity function did cause collapse'. From context it should almost certainly read 'did not'. I did not quote it.

### analog-function-crn

- YAML `origin` reads "Hjelmfelt et al.; Deckard & Sauro; Dittrich et al.", and this line is printed in the page header. The book's refs [141] and [511] are Buisman, ten Eikelder, Hilbers & Liekens (2009) and Liekens et al. (2008), not Dittrich (book bibliography, lines 22278 and 23482). `origin` is not one of the fields I may edit, so I left it unchanged. Suggested value: "Buisman et al.; Deckard & Sauro; Liekens et al.; Hjelmfelt et al."
- The page footer names tests/chemistries/test_analog_function_crn.py, but that file does not exist. The only test is test_analog_functions_compute_at_steady_state in tests/chemistries/test_w1_dynamics.py.
- The YAML has no `sources` or `phenomena` fields, so the page has no 'Sources used for the implementation' list. The four papers are cited in full only under the explainer's Further reading.
- The YAML R scheme describes only the sqrt network, although the generator offers five functions. It is incomplete rather than wrong: the other networks appear only in `decisions`.
- The YAML intuition says 'you must wait for equilibrium'. These are irreversible networks, so the answer is a non-equilibrium steady state, which is what the book calls it. The intuition is not printed on a page that has an explainer.
- Small slip in the book's §17.4.3: after listing four shortcomings it says 'This second shortcoming severely limits the scalability', but it means the fourth, composition.

### arms

- Chemart does not reproduce the gradual termination curve of [831] fig. 15. With all 380 two-symbol rules, 0/20 runs halt at every p strictly between 0 and 1 (p = 0.05 to 0.95), and 20/20 halt at exactly p = 0 and p = 1. The paper has most runs terminating below p = 0.1 and above 0.85. Random subsets of 12 or 30 rules halt in roughly the same fraction at every intermediate p (about 31/40 for 12 rules, 15/40 falling to 5/40 for 30 rules), which does not match either. test_termination_at_extreme_heating_probability tests only p = 0, 0.5 and 1, so it cannot catch this.
- [831] states 30976 rules for its two-symbol rule set, but its own eq. 4 gives 400 for n = 2, m = 5, so the paper's actual rule set is uncertain (the YAML decision already records this). This may be the cause of the fig. 15 mismatch.
- The Brusselator system (seed 1, steps=10000) does not oscillate. It halts at step 9288 with a full bag {B:212, D:2288, E:2356, Y:144}, because all four rules keep the bag's size and inputs fill it to 5000. This is consistent with the YAML decision, but the system therefore shows none of the paper's Brusselator behaviour.
- In the Ru1 system with uniform random selection, all 20 seeds (0-19) halt, because no rule consumes c. So the default ru1 run shows none of the cycles the paper shows for Ru1 under fixed rule orders (figs. 4-5, not reproducible, as the decisions say).
- The YAML top-level `notes` field is an export list for the module (run, Rule, ...), not reader-facing content. It was printed as the 'How it works' text on the old page without an explainer; the explainer now replaces it.
- The absolute numbers of Chemart's 'revisits' (about 936) and 'kinds_of_periods' (about 190) are not comparable with the paper's 'number of generated cycles' (about 450 plateau, fig. 16), because the paper does not define its counting precisely.

### automata-reaction

- YAML top-level notes (not printed on the page) claim it 'reproduces the qualitative results of AlChemy at a fraction of the compute' and call it 'the best cost/benefit constructive chemistry in the book'. Neither is in the book: it says only that the observations 'support the findings of Fontana et al.' on syntactic and semantic closure.
- Phenomenon 'lazy -> eager replicator transitions' comes only from book 10.6.1 (footnote 2), with no source cited. The 1998 paper only says the soup evolves 'towards active replication' (Fig. 4). No test checks it.
- Phenomena 'syntactic and semantic closure', 'M=10^4 exploration ended by active replicators' and 'M=10^5-10^6 self-evolution' are not checked by soup runs in the tests. Only static reaction tables are tested (Figs. 4 and 9) plus the filter mechanics. The claim that 'products resemble parents' is covered only by Table 2, not by any statistical test.
- test_passive_and_active_replicators_are_common_and_rare checks only 0.15 < passive < 0.45 and active < 0.002, a much looser bound than the ~0.005% active fraction given in the decisions.
- test_small_soup_converges_to_an_organization checks only the Hamming distance between the top two survivors. It does not check that the survivors form a closed set, although the phenomenon claims a closed, self-maintaining organization. (In my seed-0 run the four survivors were in fact closed: status=complete.)
- The later ALife VI paper (Dittrich, Ziegler & Banzhaf 1998) states filter f1 as (s1 != s2 and s1 != s3), which differs from eq. 3 of the Artificial Life paper (s1 != s3 and s2 != s3). Chemart implements eq. 3. Only a source discrepancy, noted on the page.
- Paper typo: section 5.3 says the process after gen. 40 'is similar to the extinction process shown in Fig. 3'. The extinction figure is Fig. 2.
- Paper-scale evolution runs (M=10^5-10^6, 1000-7000 generations, 10^8 or more collisions) are impractical: the pure-Python machine does about 10^4 collisions/s. M=10^4 for 60 generations already takes 45-55 s.

### autopoiesis-vmu

- **Fixed (2026-09-22):** Generator bug (chemart/chemistries/autopoiesis_vmu.py, World.try_bond, lines ~365-367): reactants/products are built as Counter({species_at(i): 1, species_at(j): 1}). When both links are the same species the duplicate dict key collapses to count 1, so the observed network records 'L1S -> L2S', 'L1 -> L2', 'L0 -> L1', 'L0S -> L1S' instead of '2 L1S -> 2 L2S' etc. (seen in the default seed=1 run: 'L1S -> L2S (x7)'). The lattice is unaffected; the recorded stoichiometry is wrong.
- Geometry artefact: with Moore-neighbourhood catalysis and a diagonally connected membrane, a catalyst in the inner corner of the Von Kamp twelve-link ring can turn a substrate diagonally OUTSIDE the ring into a link (verified: 15x15, seed 0, no decay, a link appears at (9,5) outside the ring at step 14 with the catalyst at (8,6)). Interior and exterior links can also bond across the corner. Links therefore leak out of cells.
- Long-run behaviour does not match McMullin & Varela (1997): in their set-up (15x15, 12-link cell, disintegration 0.001, 3000 steps, seeds 0-9) Chemart's catalyst is enclosed a median 122 steps with bond inhibition and 376 without; first opening at median step 72 vs 117. The paper reports cells lasting about 1000 steps with repeated repair only WITH inhibition. Chemart reproduces the single-rupture repair effect (14/20 vs 4/20 in my run) but not a lifetime benefit.
- test_boundary_is_repaired_after_one_of_its_links_decays: its whole-run half asserts only ruptures > 0 and repairs > 0, and the 'repairs' metric counts any re-enclosure of the catalyst by any closed chain. That also holds without bond inhibition (e.g. 18 'repairs' for seed 1 at decay 0.01), so this part does not check the claimed phenomenon.
- Parameters: the defaults are taken from Von Kamp (2002) table 1, the standard set of the extended SCL-DIV system, not original SCL. McMullin & Varela (1997) ran original SCL with disintegrationProbability 0.001 (default here 0.01). The original FORTRAN program EXP29.FOR let a free link bond only if at most one doubly bonded link was in its Moore neighbourhood (McMullin 1997), so its threshold corresponds to chain_inhibit_count=2, while SCL's (and the default) is 1.
- YAML decisions/sources say McMullin (1997) 'Computational Autopoiesis: The Original Algorithm' and the full McMullin & Varela (1997) paper could not be retrieved. Full texts now exist in /home/marco/.claude/jobs/42535c4e/tmp/agents/autopoiesis-vmu/ (bmcm9701.txt, bmcm-ecal97.tex/.pdf, and McMullin 2004 alj.txt). I did not change the YAML because this is a statement about what the implementation used.
- Disintegration does not conserve matter: when a decaying link has no adjacent hole, its second substrate (and any absorbed one) is lost (default run: 899 substrate units at t=0, 847 substrates + 18 links x 2 = 883 at t=120). The decisions document the 'when there is one' rule but not the loss.
- mode='reactions' is marked status complete but lists bonding only between unloaded links (L0/L1). Bonding involving loaded links (e.g. 'L0 + L0S -> L1 + L1S', observed) and disintegrations that change partners' bond states are not among the defined reactions.
- chain_inhibit_count param range text says Von Kamp's count is 'turned off at the neighbourhood size'. Von Kamp's value 6 turns it off on his hexagonal lattice. On Chemart's 8-cell Moore neighbourhood the max value 8 only nearly disables inhibition.

### bagley-farmer

- The YAML lists provides [catalysts, thermodynamic-consistency, ...], but the default network's net.summary() provides only flow, initial-state, mass-conservation, rate-constants, stoichiometry and topology. With saturation on, a catalyst is consumed into E_bound and not returned in the same reaction, so it is not detected as a catalyst. The at-a-glance table (from the YAML) and the summary printed on the same page disagree.
- Deterministic metadynamics often does not reach a fixed graph. With the exact Figure 3 constants (nu=1e4, ku=1e4, max_length=4, p=0.08, seed 2, threshold 0.01) the outcome is 'cycle-merged' after 10 rounds. The default network with threshold=0.01 also gives 'cycle-merged' for seeds 1, 2 and 3. Only the test's softened setting (nu=100, ku=1000) gives 'fixed'. So the YAML phenomenon 'the deterministic metadynamics always reaches a unique fixed point' is covered by one test case and does not hold in general in Chemart.
- The YAML decision on Table 2 says only the (b) ranges for nu and ku look swapped. Other printed ranges also exclude their own optima: in (a), ku's range [10^1, 10^3.5] excludes 5.00e4, and in (b), kr's range [10^3, 10^3.2] excludes 2.70. A simple nu/ku swap does not explain all of them.
- YAML decision: 'The paper gives no value of p'. That is true of the emergence paper [49], but the companion paper [50] (Figure 3 caption) uses p = 4.5e-3.
- With the default saturation=True, a near-equilibrium run (delta=0.01) holds 1.11 of the 2.0 units of monomer mass in bound pools, so free concentrations are 0.40-0.83 of the uncatalysed ones. Catalysis therefore changes the (free + bound) profile even near equilibrium. This follows from the paper's saturation scheme, not from a bug, but it is not stated anywhere in the YAML.
- The book's bibliography entry [260] gives Physica D volume '2(1-3)'. This looks like a typo for 22. I left the volume out of Further reading.
- The paper scan prints the eq. 12 fixed point as 'm0 = N_f/K delta'. From dm/dt = N_f delta - K m it should be N_f delta / K. Chemart uses the correct form.
- The at-a-glance row 'constructive: yes — the species set grows at run time' does not describe the default generator: it builds a fixed, complete network up to max_length. Only threshold > 0 changes the species set, and then only by selecting among the pre-built polymers.

### bigan-conservative-crn

- Book §11.3.2 says the authors stop adding reactions once only one mass vector is admissible. The paper (arXiv:1303.7439 §2.1) keeps adding until no remaining reaction is orthogonal to that vector (maximum size), and the code follows the paper. The page points out this book error.
- The generated page footer names tests/chemistries/test_bigan_conservative_crn.py, which does not exist. The tests are in tests/chemistries/test_w2_papers.py.
- The tests check only structure (a positive, unique mass vector; a loose size range over 8 seeds, 10..110 with mean 20..70; the detailed-balance ratio; the saturating K and the inflow). Neither YAML phenomenon (one chemical dominating at high density; directed transformation under flux) is tested, because Chemart has no ODE integrator for this network.
- The YAML `provides` lists rate-law and thermodynamic-consistency. The default network's net.provides has neither: rate-law appears only with kinetics='saturating', and nothing ever derives thermodynamic-consistency.
- The YAML intuition says backward rates are set 'so detailed balance holds'. The paper (§3.1) says this holds at equilibrium only for mass-action kinetics, not for saturating kinetics. The overclaim is minor and the intuition is not printed on the page, so I left it unedited.
- The default nutrient=5 raises ValueError for N<6. Users must also set nutrient when they choose a small N.
- In the YAML, A.dilution is 'external nutrient inflow drives the system out of equilibrium', and the page prints it under 'How the population is bounded'. An inflow does not bound anything: under flux one concentration grows without limit.
- My runs (seeds 0-19, flux 1 M/s, saturating) show that at the default flux some networks are in the paper's saturated regime (the nutrient A5 accumulates: seeds 6, 8, 19). On small networks the product depends on which chemical is fed (seed 1: A3, A6 or A8), unlike the paper's reference network. The 63-65 pair networks (seeds 2, 12, 13) behave as the paper describes. This is consistent with the paper's remark that maximum size is presumably what makes the behaviour robust, but it qualifies its 'independent of nutrient' claim.

### bondable-ca

- The paper was accessible after all. The earlier job's bca.pdf is a Cloudflare HTML page, but an open copy is on Banzhaf's page. The YAML decisions still say that neither [373] nor [372] was accessible. I did not edit the decisions field because it is outside the fields I may change.
- The default network's decays all come from the settling cap. With settle_iterations=256 (the default), r90+r165 has not settled when the cap is reached, and the bond run of 3 in that state triggers dissociation. With settle_iterations=2048 or 65536 the same pair settles after 1830 ticks with a run of 10, and the default seeds give only 3 associations and no dissociation (11 species instead of 23). The 11-reaction network that the test pins therefore depends on where the transient is cut off.
- Association compares snapshot polarity, not settled mean polarity. Chemistry.associate calls Molecule.polarity(), which reads the one stored configuration (the state where the cycle closed). The YAML R text says 'settled mean polarities'. Rule 110 has a cycle mean of +0.67 but a stored-state polarity of -2, so it is treated as negative. Rule 18 has a mean of -5.14 and a stored value of -8. analysis.molecule_polarity reports these snapshot values.
- Canonicalisation over chain reversal does not match the dynamics. Molecule.reversed() reverses the atom order and negates the offsets but does not mirror the cells. So a molecule and its reverse get the same species code while step_molecule runs them differently. Verified: [r90, r165, +0] and its reverse have equal canonical codes but settle to different states. This is why the default network has both 'r90 + r165 -> r165-r90' and 'r165 + r90 -> r165-r90#2'.
- In truncated runs, analysis.unsettled_molecules lists molecules that are not in the network. It is built from chem.unsettled, which includes every registered molecule. Example: polarity_rule='any' reports 1166 unsettled against 200 species.
- coupling='xor' closes over molecules that already have a bond of run 0 (for example r165-r90 with bond_strengths [0], which then associates with r0 to give r0-r165-r90). Dissociation is a separate reaction in the closure, so molecules whose bond is already broken still react.
- test_settling_times_fit_the_figure_shading only checks that settling is at most 2048 and has a maximum of 28. Chemart's settling iteration is the transient length. The figure's shading, like the paper's, spreads rules across buckets up to 2048, which suggests the time for a running mean to settle. So the shading is not reproduced, although the test name suggests it is.
- The paper contradicts itself on the bonded-pair run length: 10,000 iterations in the text, 100,000 in Table 2.
- The book's bibliography prints ECAL 2011 as '2001' (already noted in the YAML decisions).

### brane-calculi

- The bitonality check (check_wellformed in tests/chemistries/test_brane_calculi.py) runs only on the viral-infection, viral-reproduction and a custom Pino network. Its docstring and the YAML phenomenon suggest bitonality is verified generally, but eat-me, seek-and-store and viral-replication are never checked. The page says this plainly.
- viral-reproduction reaches a third reaction the paper does not show: the same cell membrane engulfs the newly budded virus again, and it then stays stuck because Z' has no endosome. The YAML decisions record this and the test asserts it. It is not a bug, but the terminal state is not the paper's end state (membrane[Z'], virus).
- viral-replication's closure is explored breadth first, so max_nucap_copies grows very slowly with max_species: 2 at 1,000, 3 at 5,000, 4 at 20,000 (about 12 s). The truncated network mostly shows different orderings of the same few copies, not the 'any number n, m' of copies the paper claims.
- The paper's two printed derivations in sec. 3.2 (Mate, Drip) contain typographical errors. The YAML decisions already record them; they are errors in the source, not in Chemart.

### brusselator

- The YAML's `book` field lists "19.3", but chapter 19 of the book never mentions the Brusselator. The only treatment is in §17.4.2. The book's own §9.4 says "chapter 19 discusses the Brusselator model in more detail", which is an error in the book, and the 19.3 entry seems to come from it. I did not fix the `book` field because it is not on the list of fields I may edit. The explainer points out the discrepancy under Further reading.
- No tests/chemistries/test_brusselator.py exists. The only test of the entry is test_brusselator_oscillates_only_above_hopf in tests/chemistries/test_w1_dynamics.py. The footer that the generator writes on docs/catalog/brusselator.md still links to the missing file.
- The YAML notes call this entry "a good test that Chemart handles arity-3 reactions and the associated l_i! factor in the k->c conversion". No test checks it: test_core.py::test_k_to_c covers only 1- and 2-molecule reactions. I ran k_to_c by hand on 2X+Y with omega=200 and it gives 5e-05 = 2/omega^2, which is correct.
- A lists an `ssa` reactor, but Chemart has no simulator. tests/chemistries/odes.py describes itself as test scaffolding, so users must write their own ODE or SSA code. The explainer includes both.
- The rate constants are numbered differently in the original paper. Prigogine & Lefever (1968) put k2 on 2X+Y->3X and k3 on B+X->Y+D; the book, and so Chemart, swap them. This is harmless but confusing when comparing with the paper, so the page mentions it.

### ca-embedded-particles

- Observed-reaction under-counting (generator bug-like behaviour). An event opened before t_c keeps absorbing later nearby changes, because `last` is re-extended while anything changes within 14 cells, and emit() then drops the whole event as pre-condensation. Real post-condensation collisions are lost with it. Reproduce with chemart.evolve('ca-embedded-particles', seed=7, lattice=75, steps=150): alpha -> gamma + mu at step 7 and mu + beta -> delta at step 9 (t_c = 5) were merged into an event that started at step 1 at site 20.5 and were discarded. Only `delta + gamma -> ∅` (x2) and `beta + mu -> delta` (x1) are reported. The page documents this.
- The decisions say the measured t_c is 6 at N = 149 for the default IC. The current code gives condensation_time = 5 for seed=1 (the default run).
- The decisions and the test docstrings compare Chemart's t_c with 'the published average of about 12'. That figure is for phi_dens5 (≈ phi_100) and uses the published transducer definition, not Chemart's 14-cell width criterion, so the comparison is loose. The test only checks 0 < t_c <= 30.
- The A.reactor is `lattice-2d`, and the generator glosses it as 'a two-dimensional grid'. The model is a 1-D ring. The catalog has no 1-D lattice reactor type, so the formal specification misleads. The page adds a clarifying sentence.
- The P_599 of phi_par^a came out at 0.795 on 200 ICs (rng seed 1), against the published 0.740: about 1.9 sampling SDs high. GKL P_149 was 0.795 on 1000 ICs against 0.816 (about 1.6 SD). Both pass the slow tests' tolerances (0.10 and 0.08), but they may be worth a larger-sample check.
- Default-run measured_velocities are unreliable (delta -1.907 over 3 tracks, eta +0.333 over 1 track, against the published -3 and +3), because they are averaged over few, short tracks. probe() gives the clean values.
- Crutchfield, Mitchell & Das (1998) contradict themselves slightly. The Table 1 caption gives the best known r = 3 rules as P_149 ≈ 0.85, while their section 9 reports Juillé & Pollack's coevolved rule at 0.863 ± 0.001. The page reports both.
- The book's §10.7.2 heading is 'Self-Replication in Cellular Automata from an AC Viewpoint'. Embedded particles is an unnumbered subheading inside it, which the YAML's `book: 10.7.2` does not make visible.

### cham

- catalog/chemistries/arms.yaml and mgs.yaml also use the reactor tag maximally-parallel. I did not check whether it is accurate for them; worth checking, given it was wrong for cham.
- YAML decision on sources: the book's [422] (Inverardi & Wolf 1995) is closed access, so the page cannot describe its chams. iw.html in the earlier job's download folder is an empty file (0 bytes).
- tests/chemistries/test_cham.py test_custom_rules_and_machines asserts that tccs and ccs-minus give identical networks on the default program. That is true, but it means the default network does not exercise any TCCS-only rule. The TCCS rules are covered only by the external-sum test and test_tccs_rule_schemata.
- Minor naming clash, not a bug: the paper's §5 'γ-calculus' (Boudol) is a different calculus from the Gamma-family γ-calculus (Banâtre et al., book [55]) that the book mentions in §9.2. The page points this out.

### chameleon

- The brief names tests/chemistries/test_chameleon.py, and the page footer links to it, but the file does not exist. The only test is test_chameleon_counts_and_invariants in tests/chemistries/test_w1_structure.py. It checks the default counts (2700/0/2700) and that all three conservation laws, the two modular ones included, hold against the stoichiometric matrix. No test checks the dynamics (relaxation to 1/3 each).
- The generator only emits the network. Chemart has no chameleon simulator and no helper for rescaling time between the multiset run and the ODE. Reproducing figures 2.5/2.6 needs user code (chemart.soup.soup plus a react function, shown on the page).
- chemart.soup.soup swap-removes the reactants and appends the products. The book instead puts the products 'into the same places where the reactants were found', which the YAML's A notes also state. With uniform random draws the dynamics are the same, but the A note describes behaviour that no Chemart code implements.
- The book's condition for when a monochrome colony is possible ('the difference between the number of chameleons of different colors is a multiple of 3') is loosely worded. The necessary condition is that the two counts other than the target colour are congruent mod 3. The page states that precise necessary condition.

### chemoton

- Nothing in Chemart or its tests runs the Chemoton's cell cycle. tests/chemistries/test_chemoton.py only checks that the rate equations match eqs. 1-14 at one random state with V=40, above the threshold, where tests/odes.py ignores the threshold gate. None of the published dynamics (replication times 0.455/0.65/2.4) is tested.
- A fixed-step Euler run (dt=1e-4) of the generated network, with the paper's volume rescaling Q=S^1.5 and division at S=2, gives generation times of 0.17 ([X]=100), about 0.25 ([X]=10) and about 0.92 ([X]=1). The paper gives 0.455, 0.65 and 2.4. The ratios agree (1.45 vs 1.43, 5.4 vs 5.3) but the absolute times are about 2.6x shorter. Peak [V] at N=25 is about 37 in the run against about 85 in the paper, and peak [A1] at [X]=100 is 1.3 against 2.5; peak [A1] at [X]=1 matches (20.3 vs 20). I found no cause. The YAML phenomenon 'replication time 0.455 ...' is therefore not reproduced in absolute terms.
- Fernando & Di Paolo (2004) disagree with themselves on the gate. The text says propagation also runs only when [V'] > [V']*, while the appendix gates only k6 (initiation). The generator and the YAML decision follow the appendix. Gating propagation too made no difference to the period in my runs (still 0.17).
- Also in the paper: the text says the time step is 0.0001, but the Figure 2/3 axes say 'Time Steps of 0.00001 units'. Eq. 10 appears to have a typo (a stray '- k9 T*' term). The tests follow the physically consistent form k9r*T - k9*T*·R.
- The only record of the Q=S^1.5 volume rule and division is text in extras.compartments, so a generic ODE integrator given this network would ignore both the threshold and division.
- At [X]=1 the Euler run settles into a period-2 cycle, alternating generation times of about 0.72 and 1.11. The paper does not mention this.

### combinator-chemistry

- Three of the eight YAML phenomena are not checked by any test: the level-1 catalytic collapse to a single self-copier (L0), the level-3 soup that settles and moves between organisations, and 'universal copiers crash diversity'. The tests check only the ladder closures, that atoms are conserved and the population varies, and the rule g + y -> 4y.
- My own runs of the level-1 catalytic soup (fig. 6.4 limits, M=100, 5000 generations): seeds 4 and 5 collapse to one molecule (WK, I), but seeds 1-3 freeze first with 7, 42 and 2 kinds of molecule between which every collision is elastic because of the size limits. So the L0 collapse is only partly reproduced.
- The default soup (method='soup', 20 generations) shows only the first crash: the population falls from 103 to about 5 molecules. Organisations appear only after thousands of generations. With seed 2, 3000 generations and 600 atoms per type, one molecule takes over and the pool runs out of R. Such a run takes 35 s to 2.5 min, so the default generations value shows nothing of the published behaviour.
- The thesis contradicts itself on level 2. Sec. 6.1.3 and 6.4 say the reaction stayed catalytic, but table 6.2 marks level 2 as Catalytic=No and Inflow=Yes. Sec. 6.1.3 also says that in level 3 'none [of the reactants] was being used up', which contradicts sec. 6.5 ('molecules are being used up when they react').
- ECAL 2001 says the 2000 paper's runs lasted 1000 generations, but the 2000 paper's figures show 300 generations (the YAML range uses 300).
- The thesis bibliography gives the ECAL 2001 paper a different title ('Metabolic and stable organisations') from the paper itself and the YAML ('Stability of metabolic and balanced organisations'). It also lists the ALife VII paper twice, as 2000a and 2000b, and credits the first combinator reproduction of AlChemy to 2000b, whereas the master thesis is 1999 ('Building life without cheating').
- The YAML source line for S2000 says type B is 'beta/gamma/delta'. In the paper, type B is generated by four molecules, including a different alpha, C(K(WC))(K(WC)). This is minor and I left it unchanged.

### conrad-enzymatic

- tests/chemistries/test_conrad_enzymatic.py::test_the_device_classifies_all_135_presentations does not reproduce Figure 9. It checks the constants stored in PUBLISHED_CLASSIFICATION, then thresholds a hand-set response {a:0, b:1, c:0}. The pass is guaranteed by construction.
- test_the_five_realisable_operations_appear_and_nxor_never_does uses an invented response shape (1+4u)exp(-u), not the measured surface. It shows that the Table 1 formulas are consistent with a unimodal response, not that MDH realises these operations.
- The network carries no rates, and the free and ion-bound conformers run the same catalytic cycle. The feature the experiment is about, that the ions change catalytic activity nonmonotonically, is not represented in the network at all. Only the Table 1 analysis layer computes anything.
- The YAML declares `provides` including sequence-structure-function, and the at-a-glance table prints it. The generated network's provides (net.summary) lists only initial-state, mass-conservation, stoichiometry and topology.
- The phenomenon 'NXOR is realisable at none, because a unimodal response has no interior minimum' presents Chemart's reasoning as the literature's. The paper only says NXOR is 'not implementable with the convention that a high response is considered to be an active output'.
- The YAML intuition says the enzyme classifies patterns 'in a milieu of two ion species'. The 135-presentation device result used MgCl2 alone. The intuition is not printed on the page, so it was left unchanged.
- Bibliography mismatch for [194], Conrad 1992, Computer 25(11): the book gives pages 11–22, the YAML refs give 11–20. Not verified which is right.
- The YAML sources say Conrad 1985 (CACM) was 'not obtained'. The earlier job's conrad1985_cacm.pdf is in fact a 5 KB HTML paywall page, not the paper. So that statement is correct, but the file name is misleading.
- The A reactor is listed as well-stirred-multiset, a pot of discrete molecules meeting at random. The experiment is a deterministic batch assay read once, and no multiset simulation exists or can run without rates.

### disperser

- The generated page footer names tests/chemistries/test_disperser.py, but that file does not exist. The only disperser test is test_disperser_converges_to_network_average in tests/chemistries/test_w1_dynamics.py (the generator assumes the file name).
- The only test integrates the ODE for 60 time units and checks X_i = 250 (rel 1e-4). It does not check the book's actual demonstration: stochastic SSA dynamics, and settling again after the +600 injection and -300 removal (averages 400 and 325).
- The YAML decision says the figure 17.3 topology is 'only drawn' and was inferred from the text. The reference program Disperser.py (PyCellChemistry 1.0) actually hard-codes topo = [[1,2],[1,3],[1,4],[2,3]]. So the choice is confirmed by the source, and the decision could say so. Not edited: decisions are outside what I may change.
- Meyer & Tschudin 2009 say the fixpoint is asymptotically stable 'for arbitrary network topologies' because 'any Laplacian has positive real eigenvalues'. That is loose: every Laplacian has a zero eigenvalue, and the equal-share result holds only on a connected graph. A disconnected graph balances each component separately; I checked with graph=[[1,2],[3,4]], which gives 500/500/0/0. The page states the connected-graph condition.
- The YAML origin says 'Meyer & Tschudin, 2009', but the book cites [576, 580], and [580] (Meyer, Yamamoto & Tschudin 2008) already contains the full Fraglets disperser protocol. Left unchanged (origin is not a field I may edit).

### dna-automaton

- conditions='fast' and conditions='parallel' only set concentrations and the cleavage rate constant. They do not select program A2 and input I3, which both published experiments used. With the defaults you get A1 on I8 at A2/I3 concentrations; the software concentrations are per molecule, so they carry over, but the network is not the published one. The page tells readers to pass program='A2', inputs=['I3'] themselves.
- Book §19.3.1 says the DNA automaton 'relies on two enzymes ... FokI ... and ligase'. That is the 2001 Nature design. The book's figs. 19.14-19.16 and its detailed reference [100] are the ligase-free 2003 design. The page mentions this mix-up.
- The YAML declares reactor [ode, ssa], but the assembly and hybridisation reactions never carry a rate, because only Kd values are published. Even under fast/parallel only the cleavage steps have rates, so the network cannot be simulated as given. Chemart also has no general ODE/SSA simulator. The page says so plainly.
- test_software_reuse_of_fig_3c checks only that T8 is the most-used rule (6 vs 3 vs 3 uses per input molecule). The published per-molecule counts (29, 21, 54) come from 18 h of recycling and are not reproduced, which is right because the model has no time course. The test name could make it sound like more is checked.
- Tominaga's stacked-string entry models a version of this automaton that has output-detection molecules (probably closer to the 2001 design). The two catalog entries overlap on the same experiment. Not investigated further.

### dna-hpp

- tools/gen_catalog_pages.py renders a parameter default of 0 as an empty cell: the dna-hpp Parameters table shows v_in's default as `` instead of `0` (docs/catalog/dna-hpp.md, the v_in row). The cause is probably a falsy-value check in the default formatter. Other entries whose int defaults are 0 will have the same problem.
- The book (19.3.1) gives the edge encoding as complemented halves with in/out swapped. That is not Adleman's scheme. The YAML decisions already record this and the explainer mentions it; it is noted here only for completeness.
- The generator applies an exact gel, so it cannot represent the 120 bp contamination the paper actually saw after step 3. The simplification is honest, but no test or decision says it explicitly. The explainer now does.
- The scaling claim (Hartmanis: 200 cities needs more DNA than the Earth weighs) is stored as a text string in extras.analysis.experiment. It is not computed, and no test checks it, although it is listed as a phenomenon.

### dorin-korb-ecosystem

- tests/chemistries/test_dorin_korb_ecosystem.py::test_no_sugar_is_made_without_sunlight: its comment says 'only sunlight can pay the high A-B bond energy, so in the dark no sugar bond ever exists'. It only checks seed 1 at the default density. The model does not guarantee this: dense dark runs make sugar from pooled bond energy, both K-catalysed and uncatalysed.
- The figure 6 decomposer body (_decomposer in chemart/chemistries/dorin_korb_ecosystem.py) anchors its O spacer to both a C and the ECC enzyme. That gives O two bonds, although its derived valence is 1 (the appendix: 'one only since O has run out of electrons'). Anchors bypass the valence check. The paper's figure could not be seen in the text extraction, so how the paper resolves this is unclear.
- In the default run all 13 sugar bonds were made by free-floating K atoms in the soup, none by the autotroph body's own chlorophyll. The 'organism' story is carried by the soup more than by the seeded bodies. This is not a bug, but readers may assume otherwise.
- Seeded body bonds are static anchors (a decision in the YAML), so the bodies never face the paper's 'natural decay' that an organic structure must pay to resist. Only bonds added during a run can decay.

### farmer-immune

- The page footer (generator) links tests/chemistries/test_farmer_immune.py, which does not exist. The farmer-immune tests are in tests/chemistries/test_w2_book.py.
- The YAML has constructive: true, so the at-a-glance table says 'the species set grows at run time'. That is true of the model but not of Chemart: the generator returns one fixed snapshot and has no metadynamics. The explainer says this explicitly.
- Decision 3 says 'The book gives no parameter values'. That is true of the book, but the 1986 paper's Fig. 3 uses l_e = l_p = 8 and s = 6, which are exactly Chemart's defaults. The decision could cite this.
- At the defaults (k1 = 1, constant k2 = 0.5) every antibody dies out: the interaction terms cancel in the total, so the total decays as exp(-k2 t). With k1 < 1 and constant k2 the ODE blows up in finite time (seed 1, M=0, k1=0.5: the total passes 1e15 at t of about 0.37). The paper's preferred scheme varies k2 to hold the total constant, and neither the generator nor any parameter offers it. So no default or parameter setting gives a regime where anything persists without code outside Chemart.
- Book 11.2.2 error: it says strings bind when 's > min(l_e, l_p)'. The paper says s < min(l_e, l_p) (the YAML decision already fixes the behaviour; the paper confirms it).
- Book bibliography [262] lists only 'J. D. Farmer and N. H. Packard'. The Physica D 22:187-204 paper is by Farmer, Packard and Perelson.
- The book says the authors 'did not report any concrete simulation results'. The 1986 paper has no figures or numbers, but it describes in words the results of preliminary simulations: k1=1 kills everything, k1<1 favours loops, and memory and forgetting.
- The book's text after eq. 11.5 refers to 'equation 11.2.2' (a typo for 11.4).
- In the formal specification, the A notes print both ODEs on one line ('... - k2 x_i y_i' = ...') because of the YAML folded scalar, which makes them hard to read. This is formatting, not a factual error, so it was left as is.
- There is no ODE for the antigens in the 1986 paper (it only says one is needed). Eq. 11.6 presumably comes from the 1987 Annals paper [261], which I could not access, so eq. 11.6 is unverified against a primary source.
- Generation is pure-Python O(N^2) matching: N=100 takes about 4 s.

### flow-ac

- chemart/chemistries/flow_ac.py membrane field: the paper's PDF (bbox positions) typesets the outward push as 0.0005 e^(5r)/r (x, y), a fraction with r in the denominator. Chemart records '0.0005 exp(5 r) (x, y)', and tests/chemistries/test_flow_ac.py::test_space_records_the_field checks for 'exp(5 r)'. The page warns readers about this; the code was not changed.
- YAML decision 3 calls V1 'an inward spiral'. As a continuous vector field it is one (eigenvalues -0.293 +/- 0.707i), but the paper's discrete update x -> x + V1(x) is exactly ((x-y)/sqrt2, (x+y)/sqrt2), a rigid 45-degree rotation with period 8 that preserves all distances. Under this reading the swirl example never mixes the population. Either the reading of the lost square roots is wrong, or the paper's patches come from local reactions rather than from transport. Not fixed, because decisions may not be edited.
- The YAML R.scheme says 'Book figure 11.23 network'. The book prints no reactions; the network is only in the paper, and the figure is reproduced from it. This is misleading but not strictly wrong, so it was left as is.
- YAML decision 1 says the E. coli example 'is not reproduced in the paper'. It means the paper does not print the network (92 species, 198 reactions); the wording is ambiguous.
- The paper is inconsistent for the swirl run: it states both 2500 starting molecules and n = 1600. This is already recorded in the decisions.

### fraglets

- The generated A line reads `compartments` ("nested compartments (membranes)"), but Fraglets nodes are flat vessels joined by segments, with no nesting. Also, the default method is `closure`, which is not listed among A's reactors.
- The YAML `provides` includes `sequence-structure-function`, but net.summary() of the default network lists only catalysts, compartments, initial-state, stoichiometry and topology.
- The YAML A note and the book (§17.3.1) say the original interpreter ran 'maximally parallel'. Tschudin 2003 (sec. III.A) says instead that when several actions are possible the system 'randomly picks one action', and adds that parallel execution is possible. I left this unchanged and credited the claim to the book.
- The 2003 paper's instruction set also has `new` and `matchS`. Chemart has neither, and neither YAML nor decisions mention them. Their absence is what rules out the paper's flow-control-with-credits protocol, not only the reasons the phenomena line gives. The 2007 tutorial marks both as 'never implemented?'.
- The quine test checks that the quine keeps forking in SSA. Because transformations are instantaneous, a settled state never holds the whole quine, and with no dilution or capacity the self-repair behaviour of the quine papers cannot be shown.

### gamma

- Book §9.2 says Gamma 'was introduced in 1990 by Banâtre and Métayer'. Gamma15 (abstract), RULE04 and the York abstract all date it to 1986 (INRIA report RR0566), and the YAML origin already says 1986. The page gives both dates.
- YAML top-level notes call Gamma 'the ancestor of CHAM, HOCL, P systems and Fraglets'. CHAM is supported (book §9.3, Gamma15 §4.1), and RULE04 lists P systems among later developments of the idea. HOCL and Fraglets are not in any source I read. The notes are not printed on the page, so I left them as they are.
- Several files that the earlier job saved as papers are really HTML error or bot-check pages, not PDFs: scp1990.pdf, scp90.pdf, cacm93.pdf, fifteen.pdf, gamma10.pdf, rr5743.pdf and rr5743b.pdf in /home/marco/.claude/jobs/a14e4701/tmp/agents/gamma/. The 1990 SCP paper and the 1993 CACM paper are still unread.
- test_sort_orders_values_by_index checks the index/value invariant only on the reactions of the default input's network. The random-sequence test checks only the final order, not the invariant.
- With default settings, program='primes' goes over its own default max_states=5000, so analysis.results and deterministic are None. It also spends about 1.8 s in that search, against about 0.25 s with max_states=0. The YAML range note says so, but the default run cannot show that the result does not depend on reaction order.
- Book eq. 9.8 (sort condition, which sorts in decreasing order as printed) and eq. 9.9 (the division rule, not Gamma's rem) are already covered in the decisions. I confirmed both with runs: the printed 9.9 rule on 2..12 reaches 6 different stable multisets, none of them the primes.

### gard

- Chemart's rate law does not match the published GARD equation. Markovitch & Krasnogor 2018, Eq. 1 (checked on the rendered PDF), is dn_i/dt = (k_f rho_i N - k_b n_i)(1 + sum_j beta_ij n_j/N). chemart/chemistries/gard.py encodes k_f rho (1-N/N_max) + sum_j beta_ij n_j rho (1-N/N_max) - k_b n_i instead. Four differences: (1) basal joining is missing the factor N; (2) the catalysed joining rate constant is beta_ij, not k_f*beta_ij, so it runs 100x faster at the default k_f, and it is not divided by N; (3) leaving is not catalysed, although the paper and the book both catalyse it; (4) the crowding factor (1 - N/N_max) appears in neither source.
- The second YAML decision says the 2018 equation 'accelerates only joining'. As printed in the paper, the catalytic factor multiplies both joining and leaving. The book (§6.2.4, 'in both directions') agrees with the paper. The decision about thermodynamic consistency rests on the same misreading.
- tests/chemistries/test_gard.py::test_rate_equation_without_crowding labels its expected value 'paper' (0.01*0.01 + beta@n*0.01 - 1e-4*n). That is the generator's own formula, not the paper's Eq. 1, so the test cannot catch the mismatch above.
- catalysed_leaving=True adds leaving at rate k_b*beta_ij, not k_b*beta_ij/N, so even with it on the network does not reproduce the paper's catalysed leaving.
- The initial_state mixes units: outside species L_i are set to the concentration rho = 0.01, while inside species A_i are molecule counts, initialised to 0. No starting assembly is given.
- A.notes says the original papers used Gillespie SSA. PNAS 2000 actually used fixed time steps (0.05 s) with Poisson-distributed count changes. The 2018 paper used Gillespie. The book attributes SSA to [764, 765], which I did not read.
- The sources entry for the 2018 paper still lists it as the source of the 'GARD rate equation', but the implementation does not follow that equation.

### high-order-chem

- Reference NumberChemHO.isprime returns False for n = 2 (the `n % 2 == 0` test runs before the `n < 3` test), so the reference's prime counts leave out 2. Chemart's is_prime counts 2 as prime, so extras.analysis.prime_fraction is slightly higher than NumberChemHO's own trace whenever 2s are present. This is not listed in the decisions.
- The book is inconsistent about the matrix chemistry file: the module list calls it MatrixChem.py, the appendix text calls it MatrixChemistry.py. This is cosmetic and does not affect this entry.

### hill-kinetics

- Book figure 18.5 is internally inconsistent: the text says both bottom panels start from G(0)=P(0)=1, C(0)=0, but the caption puts the right panel at equilibrium free protein P≈1.5. From G0=P0=1 the free protein ends at 0.55-0.62 for n=1..4 (my runs), so the right panel must use another initial protein, which the book does not give. Chemart cannot reproduce that panel as stated.
- tests/chemistries/test_w1_dynamics.py::test_hill_elementary_equilibrium_matches_hill_function uses P0=1000, where P^4/(K^4+P^4) is within about 1e-12 of 1. It checks that the gene saturates, not the shape of the Hill curve: any n>=2 would pass. A scarce-gene test (G0=0.001, P0 of 0.5 and 1.5) would check the curve itself.
- No test covers the repression variant, the lumped form (the 'hill' rate dict), or figure 18.5. There is no tests/chemistries/test_hill_kinetics.py; the brief assumed one existed.
- In the lumped form the regulator P is a species that no reaction changes, so it stays constant. That is intended, but the network gives no flag (for example extras['buffered']) saying P is an input, as brusselator does for its buffered species.

### ikegami-hashimoto

- Translation start point: the papers' prose (translate 'from a first site of the reading frame', head/tail bits first) and the book ('starting at the head binding position') disagree with the papers' Figure 1, which Chemart follows (read from the source, table bits first). One side effect: in Chemart a tape read from a given source always translates into the same machine. The reading machine matters only through what it writes. That weakens the papers' stated premise that 'a tape encodes several machines depending on which machine reads the tape'. This is covered in the decisions and now on the page, but it is a real modelling difference.
- In random-start runs Chemart almost never ends in the papers' typical minimal loop M1002/T1, although the paper says a typical configuration 'will fall into' it. In 35 runs (seeds 0-4, noise 0.02-0.13, 3000 generations, noise off at 2000) the pairs left were Mbff7/T3f, Mffff/T7f and M0000/T0, never M1002. Mbff7 also dominates the default run. These self-replicators are not among the paper's five.
- Core networks are not reproduced. In the same 35 runs, large surviving sets (24-59 machines, <muA> 0.49-0.66) appear mainly for one starting set (seed 0) at every noise level from 0.02 to 0.1, rather than in a noise window. None of the 35 runs oscillates after noise-off. No test covers core networks, Figures 3-6 or the Figure 2b oscillations.
- Self-replicator count: from at least one starting state there are 18 exact self-replicating machine/tape pairs. The YAML decision says 15, which is right when counted from state 1 only; three more (M2cc5, M4008, M744e) replicate only from state 0. The decision is not wrong as worded, but it is incomplete.
- test_noise_brings_the_published_parasites checks that M3006/M1222 appear with T5/T3 and that distinct machines exceed 5. It does not check that they invade or that populations oscillate, which is what the phenomenon claims.
- The papers do not give N or the initial populations, so the default N=1000 is a guess. This probably affects which regimes appear.

### jain-krishna

- Random phase is biased (generator bug). When lambda1 = 0, attractor() applies (C + I) only m times starting from equal populations. That gives every node a non-zero population, and the minimum always falls on nodes with no incoming link, so only those are ever replaced. In the paper, every node that is not at the end of a longest chain has X = 0 and can be replaced. Measured at m=100, p=0.0025 over 8 seeds, averaging the second half of the random phase: Chemart holds 44-69 links, while the same loop with the paper's rule holds 16-28, close to the random-graph value of about 25. The paper says links stay at the random value in this phase. Arrival time of the first ACS was not clearly different. Also, attractor_support lists all m species whenever lambda1 = 0.
- YAML decision 2 is wrong: it says iterating (C + I)^m 'concentrates on the ends of the longest chains'. After m iterations it does not (see above). I did not edit it, because decisions are outside the fields I may change.
- Tied eigenvalues (generator bug). When two loops with lambda1 = 1 sit one upstream of the other, C has a defective eigenvalue, and np.linalg.eig returns an ill-conditioned eigenvector. Seed 1, m=100, p=0.0025: a 3-cycle upstream of a 2-cycle is present from update 2001. The correct long-run support is just the downstream 2-cycle, but Chemart reports 98 populated nodes at 2001-2008, then 2 at 2009, 98 at 2010 and 2 at 2011, on the same pair of loops. This is the 'downstream' core-shift case of the Handbook's core-transforming theorem, so crash timings at lambda1 = 1 (the usual regime) are unreliable. Two disjoint loops with equal lambda1 would also get an arbitrary eigenvector.
- tests/chemistries/test_jain_krishna.py does not exist. The tests are in tests/chemistries/test_w2_book.py, but the generated page footer links the non-existent file. The tests check only the fixed point on a 4-node graph, the reaction form, and lambda1 >= 1 after 1500 updates. They do not check the random-phase selection rule, arrival or growth times, or the degenerate case.
- The at-a-glance row renders constructive: true as 'yes — the species set grows at run time'. That is false here: the species count m is fixed, and novelty comes from rewiring a node (the book's point in 15.2.2). The generator's wording for 'constructive' does not fit this entry.
- The YAML has no `sources` field, so the page has no 'Sources used for the implementation' list.
- A.notes says 'Fast: integrate to the fixed point', but the code computes the eigenvector directly and never integrates. This is harmless, and I left it unchanged.

### kappa-calculus

- The YAML decision 'Named models' says the EGFR rules are 'labelled R1..R22 in file order', but the generator labels them R1..R23 because there are 23 rules. I did not edit it, because decisions are outside the fields I may change.
- The reference guide's prose on stance ND gives both '1p' and '2p' activity k/2, but its equations give k. The YAML already records this and Chemart follows the equations; I mention it on the page.
- tests/chemistries/odes.py (the integrate helper that test_abc_dynamics_follow_figure_3 uses) lives in tests. Chemart has no public simulator, so the page includes its own solve_ivp snippet, as the Brusselator page does.
- The copy of the reference guide at the earlier job's path is dated 'August 19, 2026'. The v4 guide is a living document, so section numbers the YAML cites may drift.

### kauffman-autocatalytic-sets

- The at-a-glance table says 'constructive: yes, the species set grows at run time', but the generator enumerates every string up to max_length at once and the species set never grows. It is the papers' model that grows. The explainer points this out, but I did not change the flag (not a field I was allowed to edit).
- `provides` lists `flow` and `initial-state`, but net.summary() reports only catalysts, mass-conservation, stoichiometry and topology. No inflow is implemented, and food_set is only recorded in extras, so it does not affect the network.
- A.dilution ('food set inflow keeps the system out of equilibrium') describes the book's chemical setting, not anything Chemart implements. The explainer says so.
- The page footer points to tests/chemistries/test_kauffman_autocatalytic_sets.py, which does not exist. The tests are in tests/chemistries/test_w2_book.py (test_kauffman_polymer_counts_and_conservation, test_kauffman_catalysis_probability). They check only the reaction counts, the conservation laws and P=0/P=1. Neither phenomenon (the P_crit transition, the alphabet-size trend) is tested.
- The YAML has no `sources` field, so the page has no 'Sources used for the implementation' list.
- Book bibliography entry [260] gives Farmer, Kauffman & Packard as Physica D 2(1-3):50-67. The volume is 22 (checked against the paper and ADS).
- Farmer, Kauffman & Packard (1986) use P in two ways. The main text multiplies the number of allowed reactions by P and gives each one a random enzyme. Appendix A treats P as the probability that a given peptide catalyses a given reaction. Chemart follows the appendix and the book. The paper also includes only catalysed reactions, while Chemart also includes the uncatalysed pairs.
- The default food set {a, b, aa, bb} (book fig. 6.12) is not a complete firing disk: it lacks ab and ba. Yet the food_set parameter's text calls it 'the firing disk', and eq. 6.3 assumes a complete disk.

### l-systems

- Book §9.8 says parallel rewriting is 'much closer to the asynchrony of the growth of real cells'. ABOP §1.1 gives the opposite reason: many cell divisions happen at the same time. The book also misspells 'Prusinkiewiecz'. The page follows ABOP.
- YAML `provides` lists catalysts and sequence-structure-function, but the default network's summary reports only initial-state, stoichiometry and topology. Catalysts appear only with reading='symbols' (a -> a + b), and sequence-structure-function is never reported by a generated network.
- The `seed` has no effect: generate() never uses rng, because stochastic systems are enumerated rather than sampled. This is fine, but no document states it (the page now does).
- koch-snowflake records no published_iterations in extras.turtle, because ABOP prints no derivation length for it. That is consistent with the YAML decision, but a reader of extras might expect the field.
- With the default max_length=100000, plant-a with iterations=7 comes back truncated after derivation length 6. This is expected, but the parameter table's range ('4-7 for the plants') does not warn that higher iterations hit the length budget.

### laing-molecular-machines

- The YAML calls this 'the first molecular machine chemistry' (intuition and notes). No source read says 'first': the book does not, and Freitas & Merkle 3.8 only calls Laing 'one of the earliest proponents' of kinematicized cellular automata. Not changed, because no source contradicts it outright, but the claim is unsupported.
- The decisions and the generator docstring say H and D come from 'Laing's 1975 constituents' via Freitas & Merkle 4.8. Freitas & Merkle cite that version only as their ref [557], and their reference list returned 404, so the 1975 date cannot be checked. The page says 'the rigid-constituent version' and gives no year.
- Freitas & Merkle state that the tape is made of two-state molecules, but Laing's 1977 version (Sipper fig. 12) has a three-primitive tape: N (null), 0, 1. Chemart implements only the two-state tape. This is consistent with the decisions but is not written in them.
- Closure cost grows fast with max_length: 10 takes about 6 s, 11 about 13 s, and the earlier job's attempt at 16 hung, so I stopped it. The closure checks every ordered pair of species, including tape+tape pairs that can never react, so the cost grows with the square of the species count. The parameter table (max 10000) gives no warning.

### matrix-chemistry

- Book §3.3 says the N = 9 system has 'n_R = 261,121 reactions', which is 511² (ordered pairs including self-reactions); the book's own eq. 3.14 (excluding self-reactions) gives 260,610. The YAML's max_species range text repeats 261,121. Chemart's full N = 9 closure has 237,121 distinct non-elastic multiset reactions.
- Page numbers disagree: the preprint of [63] says Computers and Mathematics with Applications 26, pp. 109-118, while the book's bibliography and the YAML source say 26:1-8. The preprint of [66] says Complex Systems 8, pp. 205-215; the book says 215-225. I did not check which is right; the page uses the book's numbers.
- The YAML phenomena 'focusing of all mass into one self-replicator (fig. 3.3)' and 'ecosystem-like attractors (fig. 3.4)' have no dynamic test. test_13_2_1 checks only that every reaction among s2..s7 produces s1, not the ODE focusing. My runs reproduce figs 3.4 (right), 3.8, 13.1 and 13.2 with the rate equations, and figs 3.3 (right) and 3.4 (left) with the soup, but no test does.
- Banzhaf 1993 fig. 4 shows s1, s2, s4, s8 as the most common strings and says stochastic runs barely depend on initial conditions. That ODE had decay terms. Without decay, Chemart's 4-bit ODE makes s15 the most common string. This follows the book's eq. 3.24 and is documented in the decisions, but it means Chemart does not reproduce the original paper's dynamics.
- Book 12.5.2 calls table 12.2 'nontopological horizontal folding', but it equals table 3.3 (canonical folding) plus the destructor. This is already noted in the YAML decisions.
- The earlier job's downloads BF00203123.pdf and BF00203124.pdf (the Biological Cybernetics papers [64] and [65]) are 3 KB error pages, not the papers.
- Chemart has no public ODE integrator, so the explainer defines a small SciPy helper for eq. 3.24, as other explainers do.

### mccaskill-polymer-tm

- Decision 7 (reaction scheme) says the book's s1 + s2 -> s1 + s3 'is followed only in the report's optional crossed-tape variant'. But the report (scan p.15, report p.13) says that in the crossed-tape variant 'The resulting strings from both tapes could be returned to the population in addition to the initial strings', and the review says 'The resulting two strings are returned to the population'. So the book's scheme, in which the tape is consumed, matches neither processor. I did not edit the decisions, because they are outside what I may change.
- YAML origin credits 'hardware realisations NGEN/POLYP with ... Breyer, Ackermann'. By the entry's own decisions, Breyer et al. (1998) simulates other (3SR) chemistries in NGEN, not this one, so listing Breyer and Ackermann as realising this chemistry is doubtful. Left unchanged.
- The book puts the arms-race and parasite-extinction observation right after the FPGA hardware paragraph, as if it came from the hardware. In the report it is the 1988 software run (section 4). The page says so.
- Book §10.5.3 says the reactor 'consisted of a two-dimensional domain'. The 1988 report uses a linear domain of 2^M positions; the review says 2D came only in later work. This is already in the A notes.
- Soup runtime does not scale with steps: the default soup (seed 1, 5,000 steps, population 200) took about 3 s, while 100,000 steps with population 1,000 took about 2 s. Probably some early processors run to max_steps=1000 over and over. This is performance only, not a correctness bug.

### mcs-bl

- The YAML's A.reactor lists `compartments`, so the at-a-glance table and the formal specification show a compartment reactor, but Chemart does not model cells or compartments for this entry (decisions say so). A reader of the spec alone could think the cell level is implemented.
- The book's prose (11.1.1) names the dominant cell types c1, c2, c3, c4, but Table 11.2 labels them c0 to c3; it also says 'c1 consists of s1..s4, c2 of s1, s2, s5, s6', shifted by one. The YAML decisions already record this.
- Thesis 5.2 and ALife 2008 give the number of length-4 strings as '4^8 (65,536)'; with 8 symbols it is 8^4 = 4,096 (already in the decisions, confirmed by enumeration: only *$:$ self-replicates).
- In the ACS 2011 paper the Moran-process argument quotes takeovers of 1,400 to 3,150 reproductions as 'compatible' with a neutral expectation of 945 (SD 506, spread about 2 SD); 3,150 lies outside 945 ± 2 SD. The page reports the authors' judgment, not an independent one.
- Thesis 5.4 (specificity: longer-tag replicases take over more often) is runnable with the soup but not tested; the soup tests only check that fired reactions are a subset of the closure, not any published dynamics (the copier's decline and the elongation catastrophe are shown only by my runs, not by tests).
- The default soup of 10,000 collisions takes about 2.5 s; the thesis-scale runs (5 million collisions) are far beyond a quick call.

### mechanical-self-assembly

- tests/chemistries/test_mechanical_self_assembly.py does not exist, but the generated page footer links to it. The only test is test_mechanical_self_assembly_conserves_monomers in tests/chemistries/test_w1_structure.py. It checks that there are 9 reactions and that monomers are conserved. No test checks the dynamics, or the stranding of x4/x5 shown in fig. 20.3a.
- The book says parts both attach and detach when the box is agitated. Eq. 20.1, and so the generator, has association reactions only, and no decision records that detachment is left out.
- The book cites the Artificial Life IV proceedings (pp. 172-180). The widely cited version is the journal article, Artificial Life 1(4):413-427 (1994), doi 10.1162/artl.1994.1.4.413. The YAML has no `sources` or `refs` entry for either.
- The YAML intuition and notes call this 'proof' that kinetics is substrate-independent. The book only says the system 'can be described in a similar way' with kinetic equations, and that the comparison with experiment was statistically weak. These fields are not printed on the page, so I left them unchanged.
- According to Ipparthi et al., Hosokawa's model follows a probability for each possible system state; their own version of it grows like n^7 variables. Chemart uses deterministic mass-action ODEs on mean counts. This matches the book's short description (eq. 20.2) but may not match the paper's actual formulation. I could not verify this because the paper could not be obtained.

### metabolic-robot-controller

- Book §16.1.3 says 'the evolution of the artificial reaction network took place in real hardware' and that 'fundamentally the same network (up to scale factors) would evolve in a simulated and the real world'. Ziegler & Banzhaf (2001, sec. 5 and 7) say the evolution ran in a simulator and the evolved graph was then run 'completely unchanged' on a real Khepera, with only time-scale parameters changed and the sensor inflow amplified. The page reports the paper's version and notes the difference.
- Book §16.1.3 and the YAML phenomenon 'combined genetic + metabolic + signalling networks controlling chaotic dynamical systems (Lones et al.)': the models in Lones et al. (2014) are an artificial genetic network, an artificial metabolic network and a coupled network in which the genetic network controls the metabolic one. Signalling networks appear only in the paper's survey of related work. I did not check the 2010 EuroGP paper, so I left the YAML unchanged.
- The book's Figure 16.4 caption gives the Khepera six front, two side and two back proximity sensors (10 in total). The 2001 paper's Fig. 10 and eq. 29 use eight sensors: six at the front and sides, and two at the back.
- The generated network's only outflow is the actuator substance a. The 2001 lattice reactor's random dilution, which the YAML's A.dilution text mentions, is not recorded in net.outflow. So in any ODE simulation, dead-end substances (e.g. s1 at seed=1) and inputs that are never consumed (c when only the front sensors fire) grow without bound.
- Generating large graphs gets slow and can fail. Once the material-balance weight vector is unique, a random reaction is kept only if it is orthogonal to that vector, and few draws are. n_substances=200, n_reactions=200 takes about 29 s; n_substances=200, n_reactions=400 raises ValueError after 50,000 draws with 238 reactions found. The YAML allows n_reactions up to 1000 and n_substances up to 500.
- A type-7 node can have the same species on both sides (seed=1 r5: '2 s2 -> e + s2'). It is not a no-op, and the generator allows it, but the paper's four types do not clearly intend one participant to act as a pseudo-catalyst.

### mgs

- The book's primary reference [320] (ENTCS 59(4):286-304, 2001, doi:10.1016/S1571-0661(04)00293-2) has carried a CC BY-NC-ND licence since 2013, according to Crossref. The decision that says it 'could not be read' reflects a failed automated download, not a closed paper. A manual download could be used to check the reconstructed language.
- A.reactor lists 'graph-rewrite', which the generated page glosses as 'molecules are graphs and reactions rewrite them'. In MGS the molecules are not graphs. The collection (the reactor's space) is a graph, and rules rewrite paths in it. The tag misdescribes the model.
- The generated Parameters table shows the default of the bool param `torus` (false) as an empty cell (``). This is a gen_catalog_pages.py display issue, probably for every false boolean default.
- Cohen 2003 says the sieve leaves 'the prime integers less than n'. With n in the set, the primes up to and including n remain. The YAML phenomena ('up to n') and the tests (primes_upto(n) inclusive) are right; only the paper's wording is loose.
- The species-level closure for Turn gives 252 reactions, including states that never occur in a real run (the same value at several sites). This is documented as a decision, but a reader who takes the network as 'what MGS does' will be misled. The page now says so.

### michaelis-menten

- tests/chemistries/test_michaelis_menten.py does not exist (the brief names it). The only test is test_michaelis_menten_abridged_matches_elementary_when_enzyme_is_scarce in tests/chemistries/test_w1_dynamics.py. It compares the final P at t=200 (P about 1.96, roughly a fifth converted) to 2%, with E0=0.01, S0=10, ka=10. It checks the quasi-steady-state limit, but not the saturation curve, the value of km over a range of [S], or the form=abridged ka=0 error.
- Book §18.2.1 justifies setting the net rate of formation of ES to zero 'due to mass conservation'. That is really the quasi-steady-state assumption, not conservation. The page describes it as the quasi-steady-state assumption.
- Book eq. 18.4 in the extracted text reads k_m = (k_a + k_b)/k_a, with the prime of ka' lost. This is already recorded in the YAML decisions. My runs confirm numerically that (ka_rev + kb)/ka is right: the elementary rate matches vm·S/(km+S) with km=2 and not with km=ka_rev/ka=1.
- The YAML-level `provides` lists rate-law, but the default elementary network's summary does not provide rate-law; only form=abridged does. This may be intended (union over forms), but the at-a-glance table does not say so.
- The book cites Atkins & de Paula [42] for the model, not Michaelis & Menten (1913).

### molecular-tsp

- The E-machine may not match the paper. Chemart swaps two cities, copying PyCellChemistry exchangeOperator. The paper says the E-machine 'fixes two (randomly chosen) cities and inverts their order in the tour' and calls it an 'exchange-2' operation, which can also be read as reversing the whole stretch between the two cities (a 2-opt move). With the swap, E alone gets stuck on the ring. Seeds 0-2, N=20, 6000 generations: the mean stays 1.26-1.52x the optimum. The paper's Table 1a has E alone reaching the 10% criterion in 3846 generations, and the paper says the ring has 'no local minima'. In a quick check outside the generator, a single 20-city ring tour improved by random swaps reached the polygon in 4 of 8 trials, and by random stretch reversals in 8 of 8. This fits the reversal reading, but the paper's Fig. 2b example (1 2 3 6 5 4 7) is the same under both readings.
- YAML phenomenon 'ring toy problem ... the ring has no local minima' is the paper's claim, but Chemart's own swap E-machine does have local minima on the ring (see above); the tests do not check the claim.
- YAML phenomenon 'E-machine contributes most early; C and I later; recombination later still (paper fig. 5)' is not checked by any test: machine_successes keeps totals only. Recovering the time course by differencing same-seed runs (seed 0, N=30, random layout) shows recombination taking over at the end, but E leads C and I only weakly early on.
- Paper Table 2b says t_R = 1 fails with a variance breakdown. Chemart's test test_recombination_collapses_variance_on_random_cities checks only that the overlap rises (above 0.8 vs below 0.5 after 400 generations), and it asserts that frequent recombination gives SHORTER best tours. That runs opposite to the paper's failure at t_R = 1. In 4000-generation runs (N=20, 5 seeds), t_R = 1 was as good as any other setting on 4 of 5 maps and froze worse on 1.
- Paper Table 1a counts generations, and the generation length c = floor(M / sum t_j) depends on which machines are active (E alone 9 cycles per generation, E+C+I 3, all four 2). So the table's machine comparisons are not per equal number of operations. This is not a Chemart bug, but it matters when reading Table 1a and the YAML phenomenon built on it.
- YAML phenomenon 'raising t_R from 1/1000 to 1 cuts the number of generations' cites tables 1b and 2b, but in Table 2b t_R = 1 fails (variance breakdown). It holds for Table 1b only.

### music-ac

- YAML decision 'Run length' says one phrase is 'reached in about 0.2 s'. On this machine the default run takes about 2.1-3.0 s, and phrases=5 takes 2.4 s. Decisions are not mine to edit.
- Rule (9) can cut out a phrase while a Dummy placeholder is still in it, so a phrase's `notes` list can contain 'Dummy'. Example: seed=5, cadences='with-s-t', steps=3000, phrases=0 gives a phrase ending in 'Am(S) E4 D4 C4 D4 E4 B4 E4 Dummy'. This matches the paper's rule and is not a bug, but assert_phrase_is_well_formed in the tests does not check for it. The YAML phenomenon 'a phrase terminated by Start and Stop is a final product' is only true for phrases with no avoid note and no pending Dummy.
- Long runs saturate. With seed=1, steps=20000, phrases=0, all 1,400 Degree tokens are used, 167 bars are made but only 12 phrases come out. Most later collisions are a futile cycle: rule (8) joins a non-tonic bar and rule (10) detaches it again (8,897 and 8,795 firings). This does not contradict the paper, but it caps what raising `steps` can give.
- The Setomoto first name is 'Masafumi' on Tominaga's lab page. The book and the YAML give only 'M.', so this is not an error.

### nac

- YAML `provides` lists `sequence-structure-function` (and `space`), but the network the generator builds reports only initial-state, mass-conservation, space, stoichiometry, topology. The sequence-structure-function claim belongs to the active layer (node-chain folding), which is not implemented, so the claim looks unjustified.
- YAML origin is 'Suzuki, 2004-2009', but the sources go only to 2008. I found no 2009 NAC paper. A 2011 Springer chapter by Suzuki exists (doi:10.1007/978-3-642-15102-6_3). I did not change the origin because I could not show it was wrong.
- YAML sources call Suzuki (2006, Aust. J. Chem. 59:869-873) 'the journal version of [827]'. But [827] is Suzuki & Ono, 'Statistical mechanical rewiring in network artificial chemistry' (ECAL 2005 workshop CD-ROM), while the 2006 paper has one author and a different title (checked on Crossref). The link between the two is not shown by any source.
- The demixing test (test_hydrophilic_and_hydrophobic_nodes_demix) and the YAML phenomenon credit the passive-rule demixing to Suzuki (2008) and book figure 11.14(b). That figure and the 2008 abstract come from the revised model with programmed agents (centrosome, hydrogen, van der Waals), not from the passive rule. The test also checks 'all but one hydrophilic node in one cluster' on the default seed only. At n_nodes=40, mean_degree=4, 4000 steps, seed 1, 7 of the 20 hydrophilic nodes end up isolated. The rule cannot reconnect an isolated node, and clusters never merge.
- The earlier job's folder has many files named *.pdf (2004.09_ACA, AL9, 2005.03_SICE, etc.). They are HTML error pages of about 4.6 kB. Only the JSAI 2004 slides are a real PDF.
- In the rewiring method, moves that return an isomorphic cluster are counted in attempts['moved'] but not recorded as reactions (default run: 175 moved, 113 recorded firings). This is by design in the code, but neither the YAML nor the decisions document it. I explained it on the page.

### nuclear-reaction-networks

- tests/chemistries/test_nuclear_reaction_networks.py does not exist, but the generated page footer names it. The only tests are in tests/chemistries/test_w1_structure.py: conservation for all three networks, and that the default network equals book eq. 20.4.
- The YAML's A says reactor: ode, but every reaction is generated with rate None and Chemart has no simulator for this entry. The at-a-glance and formal specification suggest the network can be simulated, and it cannot.
- The 'mass-conservation' capability is baryon number (nucleon count), not rest mass, which nuclear reactions do not conserve. The label could mislead, so the page explains it.
- The intuition says the catalysts 'would vanish from a net stoichiometric matrix'. They vanish only from the summed overall reaction; each one still appears in the per-reaction stoichiometric matrix. Loose wording, not an outright error, so the YAML was left unchanged.
- The big-bang-nucleosynthesis network is forward-only with 12 reactions. The book's source (Coc et al. 2012) pairs every reaction with its reverse and compares against a 13-reaction core network (Coc & Vangioni 2010), so Chemart's 12 reactions are not either published network.
- No test checks closure or self-maintenance of the CNO set, although phenomenon 1 claims it.

### okamoto-switch

- YAML `origin` says 'Okamoto, Sakai & Hayashi, 1987-1993', but the cited Okamoto papers run from 1980 ([632]) to 1990 ([631]); PubMed turns up nothing from 1993. I did not change it because `origin` is not one of the fields I was allowed to edit.
- YAML intuition and phenomena say A and B 'flip within seconds of the inputs crossing'. That is the book's figure, but Chemart's network flips 18 s after the crossing at the defaults and 73 s after at k_conv=0.0015. The flip always comes at about twice the crossing time: it happens exactly when the integral of k_in(I1-I2) dt returns to zero (predicted 35.8/71.5/143.0 s, simulated 35.8/71.5/143.1 s). The decision already notes the lag at the defaults but not that it grows as conversion slows. The intuition is not printed on the page, so a reader never sees it.
- The phenomenon 'bistability' is not what the printed six-reaction scheme does. With held inputs (k_conv=0) there is no steady state at all: the substrate of the larger input grows without limit (X1 = 600 at 30 s, 1200 at 60 s). The switch works through a stored surplus, which is a memory of the input history, not through two stable steady states.
- test_okamoto_switch_flips_after_inputs_cross checks only that the flip comes after the crossing and that a slower conversion delays it. It does not check the book's claim that the flip comes 'a few seconds' after the crossing, and it would pass with any lag.
- Possible mix-up in book §17.4.1: it says [635] is the McCulloch-Pitts neuron 'enhanced with memory storage' and that logic gates were 'briefly evoked in [634]'. By the PubMed abstracts, [634] (Biol. Cybern. 1988) is the paper about the McCulloch-Pitts equation and the 'mnemonic mechanism', and [635] (1989) is about setting the switching time by a pulse. I checked abstracts only, not the full texts.
- YAML `provides` lists `catalysts`, but no species in the network is a true catalyst. A and B are consumed and regenerated only through the two-reaction cycle.
- The book prints the input reactions as I1 -> X1 and I2 -> X3. The generator writes them as I1 -> I1 + X1 and I2 -> I2 + X3, and the YAML R text keeps the book's form. A decision documents this, so it is consistent, but the formal specification on the page does not match the generated reactions.

### ono-ikegami-protocell

- tools/gen_catalog_pages.py param_rows: a parameter whose default is 0.0 appears in the parameter table as an empty `` cell (seen for M_fraction, default 0.0). The falsy default seems to be lost in esc(p.default, code=True).
- YAML param `mode` meaning says "'reactions' returns the four-reaction network", but the network has six reactions. The test name test_reaction_set_is_the_book_s_four_reactions has the same miscount. I did not edit it because params are off-limits.
- YAML param `steps` meaning says 'each sweep does six Metropolis exchange passes', but the number of passes is `relaxation`, which defaults to 12. I did not edit it because params are off-limits.
- YAML `origin` is 'Ono & Ikegami, 1999-2003', but one of the entry's refs, [637] (Ono, BioSystems 2005), is from 2005 and has Ono as sole author. `origin` is not a field I may change.
- YAML top-level `notes` claims this is 'why the book uses it to close the protocell discussion rather than the origin-of-information one'. The book gives no such reason, so this is an unsourced motivation. It is not printed on a page that has an explainer; I left it alone.
- The energies['anisotropy_field'] string in the generator says M_a puts its repulsion 'onto the two faces normal to its orientation o'. The docstring and YAML say repulsion is strongest along the ±o axis. The wording is at best ambiguous.
- The 2000 J. theor. Biol. paper (book ref [640]), which is the implementation's main source, is in the book's bibliography but is not cited in the text of 6.3.2. It is also missing from the entry's refs list, although it does appear in sources.

### oregonator

- No test file tests/chemistries/test_oregonator.py exists, but the generated page's footer links to it. The only test for this entry is test_oregonator_real_f_keeps_rate_equations in tests/chemistries/test_w1_structure.py, and it checks only how a fractional f is written as reactions. No test integrates the dynamics or checks either catalogued phenomenon.
- The default network's initial_state holds only A = B = 1, so X, Y and Z start at 0. Every reaction rate is then zero and a simulation started from the defaults never moves; the user has to seed Y or X. The page explains this, but the default could seed something.
- The YAML has no `sources` field, so the page has no 'Sources used for the implementation' list.
- The book's scheme (B + X -> 2X + Z; Z -> f Y) differs from Field's standard form (A + X -> 2X + 2Z; B + Z -> (f/2) Y, with A = bromate and B = malonic acid). So in the book, and in Chemart, B feeds the autocatalytic step, where in the chemistry bromate does. With A and B buffered the two forms have the same dynamics up to scaling. The page explains this, but the YAML's S text still calls A and B simply 'buffered reactants'.
- The book's rate constants (k1 = k5 = 1, k2 = k3 = 10, k4 = 2.5) are far from the chemical ones: in Field's scaling they give eps = 0.1 and q = 0.05, against about 0.0099 and 7.6e-5 for the real chemistry. My stability analysis puts the oscillation window at only 0.881 < f < 1.458. With these constants a stable steady state answers a push in proportion to its size, with no threshold (tested at f = 2), so the default constants may not give a clearly excitable medium.
- Scholarpedia (Field 2007) dates the FKN mechanism to 1974 in its text, but its own reference list gives 1972 (JACS 94, 8649). The page uses 1972.

### p-systems

- The wrong Paun 2006 section number (sec. 11 where it should be sec. 12, for strong/weak priorities) also appears in places this task was not allowed to edit: the 'priority' param meaning in p-systems.yaml, the Priorities decision ('weak is available for Paun 2006 sec. 11'), the docstring of select() in chemart/chemistries/p_systems.py, and the docstring of test_paun_strong_and_weak_priorities.
- The range text of the 'n' param says the divisibility computation halts 'after about n/k + 3 steps'. Runs show ceil(n/k) + 1 steps when k divides n and ceil(n/k) + 2 otherwise: n=7,k=3 takes 5; n=9,k=3 takes 4; n=2,k=5 takes 3; n=1000000,k=7 takes 142,860. It is not a param I could change.
- The Paun 2006 source entry puts 'the n^2 system of fig. 3' in sec. 9. Fig. 3 itself is in sec. 8 (line 564 of the preprint); the computation that walks through it is in sec. 9. This is a minor imprecision and I left it unchanged.
- The earlier job's jcss.pdf, guide.pdf and gentle.pdf in /home/marco/.claude/jobs/a14e4701/tmp/agents/p-systems/ are HTML pages (Elsevier landing pages), not the papers. This matches the decision that those originals were not reached.
- The earlier draft said the n=1000000, k=7 run takes about 17 s. My run took about 10 s. The page says 'about ten seconds'.

### prime-number-chemistry

- Book §17.1.1 (p. 349) says the prime number chemistry has a deterministic outcome: 'The computation produces the same result whatever the order we choose to pick the molecules.' With the catalytic rule of eq. 2.43 / NumberChem.py that is false for the exact final multiset. In 200 Chemart soup runs from {2, 3, 12}, 157 ended as {2,2,3} and 43 as {2,3,3}. In small soups, whether the run stalls also depends on the order. Only the kind of result (all primes when the run does not stall) is order-independent. The explainer says so.
- The source pack says 'Section 1.4 / 1.5 (not found)'. The YAML's 'book: 1 (eqs. 1.4-1.5)' is correct: these are equations 1.4 (8:4→2) and 1.5 (8:5→) in chapter 1, book.txt lines ~590-620. The pack generator read them as section numbers.
- The YAML intuition (not printed on the page) calls the chemistry 'the sieve of Eratosthenes as a reactor'. No source says this, and the mechanism differs: composites are divided, not crossed out. I did not change it because it is an analogy, not a factual claim, but it may be worth softening.
- The M param range says 'paper fig. 6 scans 10-200'. The figure's axis runs 0-200 and the lowest plotted soup size cannot be read from the text extraction, so I could not verify the 10. Left unchanged.
- At M = 100, maxn = 10000 and 700 generations, Chemart ends all-prime in 4 of 10 seeds (mean 0.91). The paper's fig. 5 shows a single all-prime run at M = 100, and its fig. 6 says M > 100 is nearly always 100%. That is consistent at the boundary but not a close match. No test checks M = 100.

### proof-ac

- YAML decision on the group-theory problem is wrong about the closure. It says 'level saturation with 2000 species only covers levels 1-3 (about 14 s)'. My run of generate_network('proof-ac', problem='group-right-inverse') at the default max_species=2000 took 3.8 s and hit the budget partway through level 2: 6 clauses at level 0, 183 at level 1, 1811 at level 2, none at level 3. I did not edit it because decisions are outside what I may change.
- The default closure's reported proof (analysis.proof) is 8 steps and goes through level 4. It is not the 7-step published proof from thesis table 6.1, which is only checked as a set of valid steps (test_published_proof_of_table_6_1). This is consistent with the decisions ('derivation of lowest level'), and the page explains it.
- Thesis fig. 6.1: the YAML param range says multiplicity '1-7000', but the extracted axis runs 0 to 7000. Minor, and I left it alone.
- ecal.pdf in the earlier job's folder is a 260-byte stub, so the Busch & Banzhaf ECAL 2003 paper was not available. This matches the YAML decision.
- The thesis also has two Chemart-absent experiments, Burnside (6.7) and Schubert's Steamroller (6.8), and neither is in the YAML phenomena. They are now covered in Results as not reproduced.

### raf

- Minor wording gap, no YAML change: the YAML decision calls inhibition (u-RAFs) 'NP-hard in general'. Hordijk & Steel 2012 (extended2012.txt l.333) say Mossel & Steel 2005 proved it NP-complete, and the book says NP-hard. Both are consistent (NP-complete implies NP-hard); the page gives both.
- Not a bug, but a reader could misread it: the 2011 definition lets reactions whose reactants and catalyst are all food count as RAFs. At n=8, f=1.0, 3 of 20 seeds have a one-reaction maxRAF '1+0<->10' (in seed 2 it is catalysed by the food molecule '11'). So `raf_exists` is True for RAFs the papers call trivial. The generator does not flag non-trivial RAFs (the 2011 paper's check: does any maxRAF reaction have a catalyst outside F). The page explains this.
- The book's reference [599] (Mossel & Steel 2005) is not in the YAML `sources`, although the book cites it for the inhibition NP-hardness result and it is the analytical confirmation of linear growth. I added it to the explainer's Further reading instead of editing sources, since it was not used for the implementation.

### random-catalytic-networks

- Chemart's random ensemble differs from the paper's. Stadler et al. give each reacting pair exactly one product, chosen at random (eq. 25), and control sparseness with p_el and copying with p_self. Chemart keeps each (pair, product) reaction independently with probability `density`, so a pair can make several products or none. `allow_direct_replication` only adds copy reactions alongside the others and does not match p_self. As a result, the paper's Figure 1 cannot be reproduced quantitatively, and no setting gives the pure replicator special case. Neither the YAML decisions nor the book mention this.
- The book (§7.2.9) says the random mix converges to a fixpoint where 'only a subset of species has survived'. In the paper, dense random networks converge to an interior fixed point where all species survive; a subset survives only in sparse networks. The book oversimplifies here.
- Many sparse Chemart networks (density ≤ 0.05, n=10) have no self-sustaining set. The dilution flux phi then decays towards 0 (about 1/t) and the system never settles; the paper calls such a system not 'active'. The single existing test does not cover this, and nothing warns users about it.
- There is no dedicated test file. The only test is test_random_catalytic_network_equation_and_no_direct_replication in tests/chemistries/test_w2_book.py. It checks eq. 7.32 and that no reaction makes one of its own reactants. No test integrates the network or checks convergence. The generated page footer still links a non-existent tests/chemistries/test_random_catalytic_networks.py: the generator's footer does not check that the file exists.
- The YAML has no `sources` entry, so the page's 'Sources used for the implementation' block is empty. The paper (Physica D 63:378-392) appears only as book ref [792].
- Paper preprint internal inconsistency: the life-cycle theorem (§6.2.2) says the fixed point is asymptotically stable 'if n ≤ 4', but the paper's own inequality (33), and my eigenvalue check, give stable for n ≤ 5, nonhyperbolic at n = 6 and unstable for n ≥ 7. The page mentions this.

### rbn

- catalog/chemistries/rbn.yaml R.arity is 2 (rendered as 'molecules per reaction: 2'), but the classic-RBN reactions have k+1 reactants (3 at the default K=2, e.g. x0_0 + x8_0 + x1_0 -> ...). Only RBN World reactions are binary. I did not change it because arity is not one of the text fields I may edit.
- In RBN World, a state left behind by a failed bond (e.g. B.2) goes back to its empty-site attractor (B.1) the next time it collides with anything, because every collision re-runs the atom's dynamics. So reactions like 'B.2 + E.1 -> B.1 + E.1' (10 of the reactions in the default seed-1 run) are recorded as catalysed by whatever the partner was. This follows the thesis procedure (6.2.1 and 6.2.6), but it inflates the apparent 'catalysis' count. I checked it by re-settling B.2 through World.settle.
- The ALife XII paper (Faulconbridge et al. 2010) says 'Less than 5% of alternative chemistries pass all the tests', but its table 5 lists 19 of the 200 tested chemistries (9.5%). The page quotes the 19 and the 200 and leaves the percentage out.
- With the default proportion-sum-one rule and random atoms, the RBN World networks are very small (seed 1: 7 species, 9 reactions; only element B can bond, and only with itself). By contrast, cycle-length-equal pulls almost all atoms into one molecule (seed 1: 98 of 100 atoms). This is not a bug, but a user may find the default run thin.

### reflexive-ac

- Reaction objects store the two reactants as a multiset, so the network loses which reactant was the sender. With machines=['M45','M61'], 89 reactant pairs each show up as two reactions with the same reactants and different products, and nothing in the network tells them apart. The page says this and points readers to compose(). The design decision 'Ordered pairs: s1 o s2 and s2 o s1 differ' cannot be read back from the Network.
- The M45 sequence in the phenomena, '2, 4, 16, 32, 64 states', mixes two different constructions. Self-composition gives 4. The chains M45∘M45∘M45... give 16, 32, 64, 128 (a chain of n distinct species has 2·2^n states because the product tracks which machine is sending). This is correct, but the list reads like one repeated operation.

### repressilator

- The generated page footer links tests/chemistries/test_repressilator.py, which does not exist. The repressilator tests are test_repressilator_is_cyclic in tests/chemistries/test_w1_structure.py and test_repressilator_needs_cooperativity_to_oscillate in tests/chemistries/test_w1_dynamics.py. The source pack also reported '(no test file)' for this reason.
- The YAML has no `sources` field, so the References part holds only the book citation. Elowitz & Leibler (2000), Loinger & Biham (2007) and Wolkenhauer et al. (2004) are now cited in the explainer's Further reading instead.
- No test covers the phenomenon 'oscillations survive stochastic noise at low gene copy number'. The tests check only the ODE with n=2 (oscillates) and n=1 (settles), plus the ring wiring.
- The book (19.3.2) says the oscillating parameter region 'can be determined by a stability analysis of the system [256]'. That analysis in Elowitz & Leibler Box 1 is for their own model (Hill-function repression with leakiness, time rescaled to the mRNA lifetime), not for the book's explicit gene-binding mass-action scheme. The book's rate constants are also not the paper's measured values, although the book's figures label time in seconds.
- The first eight reactions printed by the generator are all binding/unbinding plus gene 1 transcription/translation, so the default page shows none of the decay reactions. This is cosmetic only.

### rna-folding-ac

- YAML decision 6 says the defaults are '30 nt, 4 seed sequences, 60 species', but the `pool` param default is 8, and the default run (and test_default_closure_shape) uses 8 seed sequences. I did not fix it because the `decisions` field is outside what I may edit.
- YAML intuition calls this 'the most physically real chemistry in the catalog'. That is an opinion no source supports. It is not printed on a page that has an explainer, so I left it alone.
- A.reactor lists `graph-rewrite`, but Chemart's generator does no graph rewriting (only the original ToyChem toolchain does). The page could suggest otherwise.
- The generator picks the longest loop with max() over (size, stems, unpaired) tuples, so a tie goes to the loop with more stems (hairpin vs interior → ligation). No decision records this, and the default network's only ligase (CGCACGCUCGUUCAGGUCCACGUUAGUCCU) is exactly such a tie. The explainer now says so.
- The default closure contains a ligase but no ligation (all 8 reactions are cleavages), so the default network never shows the second reaction type.

### sac

- Model (iii) soup does not show the published behaviour. Copier string MODEL_III_COPIER[2], which decodes to [\L\0\0]$[\R\0\0]$[00*0011 → 00*0011.''M"\'"\'MMMM], releases the by-product ''M"\'"\'MMMM from every gene or gene intermediate it meets (54 of 69 observed species with seed 1). In Chemart's decoding this by-product, and ''\""\'"\'M\\\E\\\M, react with nothing. Result: with seed=1, organisation=spindle-membrane, method=soup, the cell has 79 strings after 20,000 steps, 60 of them this by-product, and max_membrane_M stays 0. After 100,000 steps it has 124 strings, 105 of them the by-product, and no L/R-tagged copy has appeared. The test test_model_iii_spindle_and_membrane leaves copier[2] out of its drive, so no test catches this. It is probably the doubled-prefix reconstruction flagged in the decisions.
- Suzuki's ECAL 2003 slides (p. 9) say 'Seventeen elementary characters are prepared', while book table 11.4, the YAML and the code use 19 symbols. The difference is not explained; possibly the slide excludes L and R, but that is unverified.
- The refs entry doi:10.1007/978-3-540-39432-7_9 could not be checked as Suzuki (2003): the earlier job's download of it is a 'Client Challenge' HTML page, not the paper.
- The generated 'at a glance' row says fidelity 'reconstructed — built from the original papers listed under References', but for sac no paper text was available: everything comes from the author's slides. This is generic generator wording, not a YAML fact.

### self-propelled-droplets

- A YAML decision and the `source` field of extras.interaction_law.chemotaxis (in the code) both say that Hanczyc & Ikegami (2010) section 3.1 states the droplet moves toward the HIGHEST pH. Section 3.1 of the recovered text says only 'pH-directed chemotaxis'. The 'toward the highest pH' wording comes from the book alone (§19.2.5). Horibe et al. (2011) say 'climbing chemical gradients'.
- A YAML decision says Hanczyc & Ikegami's appendix writes '100 mL' and '800 mL', which it calls physically impossible. The same pdftotext extraction renders µ as 'A' or 'm' elsewhere in that paper ('100 Am', '1, 5, 10, and 30 AL'), so the 'mL' may be an extraction artefact rather than an error in the paper. I could not check this against the PDF glyphs.
- Hanczyc (2014) writes: 'The maze solving droplets represent self-moving droplets of type 1 (with no reactive chemistry)'. The code labels maze-chemotaxis 'no onboard reaction fuel' but gives it two acid/base reactions, following Lagzi's abstract ('interplay between acid/base chemistry and surface tension'). The two accounts are not reconciled in the decisions, though the sentence in the review is ambiguous about which maze droplets it means.
- For system=fuel-surfactant and maze-chemotaxis, the initial state reuses precursor_M (default 0.5 M, the oleic-anhydride loading) as the concentration of 4-octylaniline or 2-hexyldecanoic acid in the oil (500 mM). It also reuses surfactant_mM for the fuel precursor or the carboxylate, and the default pH 11 carries over to the maze system. None of these values is published for those systems, so their initial states are placeholders rather than sourced values.
- The YAML lists `catalysts` in `provides`, but the default (oleic-anhydride) network does not provide it: net.summary() shows no catalysts, and only the fuel-surfactant system has a catalyst. The at-a-glance table therefore shows `catalysts` for the default page.
- Horibe et al. (2011) mix oleic anhydride and nitrobenzene 1:1 v/v, not at 0.5 M. The size_series fuel budgets use the 0.5 M default, so they do not match Horibe's droplets. The explainer says so.

### smn

- Book §18.3.1 has a typo in its text: it says 'Figure 18.3.1 compares the distribution of community sizes', but the figure is 18.6.
- The papers folder from the earlier job (/home/marco/.claude/jobs/a14e4701/tmp/agents/smn/) has no usable paper. springer_content.pdf, springer.html and flamm2010.pdf are all the same 3 KB HTML login page. strings_biorxiv (Moyer et al. 2020) never cites Ono et al. The ECAL 2005 paper is still closed access: Semantic Scholar marks it CLOSED and hides the abstract.
- Crossref gives the authors' full names as Naoaki Ono, Yoshi Fujiwara and Kikuo Yuta. The YAML and the book give only initials. This is not an error, just extra detail I used in the page.
- Chemart has no public ODE integrator. The only one is tests/chemistries/odes.py, a test helper, which is why the page cannot show a recipe that runs a metabolism and reports its letter mass.

### soas

- The YAML's `provides` lists sequence-structure-function, and the page's 'at a glance' table repeats it, but the generated network reports provides = [catalysts, initial-state, stoichiometry, topology], and test_network_shape asserts exactly that list. I did not touch it: decision 10 of the YAML explicitly argues that sequence-structure-function belongs in the list.
- Book §20.1 says the assembly line is 'optimized for a given manufacturing order'. The 2012 paper says the model 'only searches for a viable solution, without giving any importance to performance characteristics or efficiency', and the 2010 paper says it does not optimise for performance or for using fewer modules 'yet'. The page points this out; the book is imprecise here.
- In its reference list the 2010 JAmI paper cites Wermelinger (1998) in 'IEEE Proc Softw'. The journal's actual name is IEE Proceedings - Software, and Further reading uses that name.
- In the test file, test_2010_trace_of_the_two_part_tape_roller cites 'figs. 22-26' and the YAML sources cite 'figs. 15-26' for sec. 6.5. Both are consistent with the paper; noted only because the figure numbers were not otherwise checked.

### squirm3

- The param range for flood_period says '2007: 50000 and 200000'. Those values come from the 2004 ALife IX paper (alife9.txt lines 190-222) and the C++ default FLOOD_PERIOD = 50000. Hutton (2007) uses 8,000 iterations (experiment 1) and 30,000 (experiment 2) (cells2007.txt lines 924 and 988). I did not change it because params are off limits.
- The replicator does not reach the 2002 experiment 1 result. The paper reports 11 copies after 3,544 steps: 8 finished and 3 stuck halfway, with the e0 atoms used up. I ran seed_molecule='e8-a1-b1-c1-f1' with steps=3544 on seeds 1-5. The molecule split 6, 1, 1, 2 and 5 times, but only 0 or 1 finished copies were left at the end, and e0 was still available. Most copies end stuck mid-copy or tangled, and 10,000-step runs look the same. The test only checks that at least one split happens. The YAML notes' claim that the seed 'replicates within a few thousand steps' holds only in this weak sense.
- The 2007 cell never finishes dividing and never reads out an enzyme. I ran it on a 20x20 grid with 150 food atoms: seed 5 for 2,500 steps, seeds 1 and 2 for 3,000 steps. R1-R34 all fire, R34 included, but each run ends as a single bonded structure of 38-46 atoms, not two cells. R36 fires, but R37 (the base read-out) fires 0 times in all three runs, so no enzyme is ever made. This may be linked to the next problem.
- One decision says bond crossing 'matters only for membranes, not for the replicator'. But Hutton (2007), in the caption of Fig. 4 (cells2007.txt around line 326), says the no-crossed-bonds constraint is what stops c6 reacting with the wrong b7 during the 2007 base duplication. So the constraint also matters for gene copying under the membrane rules.
- Flooding differs slightly from the 2005 and 2007 papers. Both papers also dissolve the atoms bonded to flooded atoms: aca05 lines 207-209 ('the dissolving state is assigned to all the atoms that were bonded to that atom'), and cells2007 ('setting the states of the atoms and those connected to them to 0'). Chemart follows the C++ DoFlood, which only breaks the bonds of atoms inside the sector. The decision cites the 2005 'dissolving state' as if the two matched. Also, the quarter rotation order differs from the C++: Chemart goes TL, TR, BL, BR; the C++ goes TL, TR, BR, BL.
- method='closure' with the default max_species=120 did not finish within 5 minutes (it was killed by the timeout). The parameter default is not practical.
- Hutton (2002) p. 3 says the type variables stand for 'any of the other states a - f'. 'states' should read 'types'. This is a slip in the source paper; I did not quote it.

### sr-loops

- Code bug, chemart/chemistries/sr_loops.py Colony.observe: when a new loop appears detached (orphan), the first loop over `newborn` sets a local `parent` to the nearest living loop, but that never updates the `newborn` tuples. The second loop then still sees parent=None and records `∅ -> L`. So the YAML R notes ('or the nearest loop') and the decision ('the nearest living loop') do not match what the code does. Verified by instrumenting record(): with 9 to 29 loops alive, evoloop 200x200 for 8000 steps gave 56 `∅ -> L` events; the SDSR recipe gave 33 `∅ -> L086aaa`; the 30,000-step evoloop run gave 119 `∅ -> L149aaa`. The page warns the reader about this.
- Loops that touch merge into one 8-connected blob, and the loop that loses the merge is recorded as `L -> ∅`. So Langton runs show deaths even though Langton loops never dissolve: 49 on 100x100 over 4000 steps. Chou-Reggia runs also register fake variants, e.g. reggia-2 at 40x40 over 80 steps: `L005aaa -> L005aaa + L004aab` x15 and `L005aaa -> ∅` x24. The YAML R note says a death 'in the SDSR loop and the evoloop is Sayama's structural dissolution' but does not say that deaths are also merge artefacts.
- extras.analysis.ancestor_copy_steps is [] for rule=byl, and the replication records have copy=None. At step 25 the Byl mother and daughter are still 8-connected, so the blob never equals the ancestor configuration. The published 25-step period is only visible through the test's S.copies().
- Chemart's evoloop ancestor (Golly's Evoloop.rle, 149 cells) has 13 state-7 straight-growth signals and an 11x11 empty middle, so it is a size-13 loop. The YAML phenomenon ('a size-8 ancestor is replaced...') and the test comment ('far smaller than the size-8 ancestor') describe Salzberg et al.'s Figure 7 run, not Chemart's ancestor. The page now says this.
- In macro mode a species name records the cell count when the species was first registered, which is usually a half-built loop just after separation (e.g. L033atf), not the adult loop's size. The evoloop test's assertion `int(dominant[1:4]) < 60` therefore measures birth size, not loop size.
- Salzberg et al. (2004) contradict themselves on the C letter: the text says C is a single core state '1', while the Fig. 1(c) caption maps '0' states to C. The YAML source line follows the caption ('C = 0 filler'). Not changed.

### srsim

- Published scaffold geometry (ex/001.scaffold/*/spass.geo): all four S sites have theta=0 with phi 0/72/144/216. Because theta=0 is the pole, all four point in the same direction. The paper says so itself ('all the component vectors of S face towards one pole') and cancels the effect with AngularDeviation 180 and fAngle 0. Chemart copies this correctly, but anyone who runs model='scaffold' with k_angle>0 or a smaller angle_tolerance will pile every ligand onto one pole.
- The bond_angle param range text ('intermediate angles rings and helices') and R.notes ('rods, rings or helices depending on the angle ... (paper fig. 3)') go beyond the paper. Figure 3 makes helices from out-of-plane inclinations, and the paper says telling squares and helices apart from a worm-like chain needs dihedral angles, which neither SRSim nor Chemart implements. I did not edit these (params cannot be changed; the notes wording is loose rather than clearly wrong).
- With bond_angle=180 and angle_tolerance=15, seed 1 and the default 1500 steps give only dimers (23 dimers, largest complex 2). Even 6000 steps reach only trimers. The paper's 'rods' are not visible at default run lengths. The page says so.
- extras['rules'] prints the scaffold rule rates without rate_scale applied (e.g. phosphorylate @ 0.00047), while extras['kinetics'] k_micro includes it (x100 in the test). This could confuse readers of the BNGL text.
- I ran one read-only `git diff --stat` on catalog/chemistries/srsim.yaml to confirm my edit, although the task said to run no git commands. Nothing was changed by it.

### stringmol

- YAML A.notes gives the bind propensity as 1-(1-a/c)^n. The code (Container.__init__ area ratio), the cell_radius param text and the decisions all use (agent_radius/cell_radius)^2. Reading a/c as a ratio of radii gives the wrong value, so the notes could say (agent_radius/cell_radius)^2 explicitly. I left it unchanged because it is ambiguous rather than clearly wrong.
- YAML R.notes says '> move (I ... to F)', but the source and port move I to F+1 (exec_step; the decisions state this). The formal specification on the page is slightly inaccurate here.
- ALife XII paper (Hickinbotham et al. 2010) has two slips in its text on the origin of species 31. It says 'the new molecule (species 31) is created from most of a molecule of species 29 with a copy of species 9 pasted', which should be species 30. It says 'When species 9 binds to species 31, the bind site is shifted', which should be 30. Its Fig. 8, the YAML and Chemart's tests all use the correct species numbers.
- ALife XII paper says the copy of 9 is pasted 'over the penultimate symbol' of 29. The test (SP30 == SP29[:64] + SP9, where len(SP29) is 65) shows the copy is actually written over 29's last symbol. The page follows the test.
- ALife XII paper's hypercycle subclass counts (emergent 8, spontaneous 15, multispecies 14) add up to 37, more than the 30 trials it says had hypercycles. The page reports the numbers as printed.
- The YAML decay param range says 'spec v0.2 and ALife XII: 1/65^2'. Only the spec states 1/65^2; the ALife XII paper text gives no decay rate.
- The default container run (cell_radius 2500, 100 seeds, 3000 steps) is limited by binding, not by energy. Unspent energy piles up to 49,174 units and the population reaches only 123, far from the ~350-molecule steady state of ALife XII. This is not a bug, but the default does not show the energy-decay balance the entry describes.

### synthon

- The YAML declares `provides: [..., sequence-structure-function]`, but the default net.summary() does not list it, and it does list `catalysts`. The catalysts come from G5 reactions where a species appears on both sides, e.g. `HO:*H + HO:**H -> HO::H + HO:*H`. This is a possible mismatch between the declared and the computed capabilities.
- The 2005 report credits the synthon model to 'Kvasnicka and Koca' (Koca 1988a,b; Hladka et al.). The YAML sources cite only Koca. This is not an error, but it is incomplete.
- In the YAML's Results-facing text, table 1 G7 (A+ + e- -> A) is filed under 'dissociative recombination'. That follows the table layout, which prints no class name for that row, but chemically it is radiative recombination. I checked the row in the PDF (G7 A+ + e- -> A, 1e-11 cm3 s-1); code and tests match it.

### tominaga-stacked-strings

- Timing drift: the YAML max_species range and decisions say the 5000-species Adleman closure takes about 34-35 s. On this machine it took 62.3 s (2000 species: 6.3 s against the stated 5 s). This may be machine load, but the numbers in the catalog could be softened to 'under a minute or so'.
- Not a bug, but a reader might be surprised: the Benenson closure contains side reactions the published rules allow, such as rule (10) ligating a spent Fok I head (0#GGATGGCGCAGCTG/0#CCTACCGCGTCGACAGCG/) to the terminator sticky end. The paper says the head 'will not be used'. The analysis (accepted words) is unaffected. The explainer mentions this in one sentence.
- Figure 3 of the 2009 paper: the start molecule 0#HHHHHO/0#CCCCCC/0#HHHHHSCoa/ appears in the pdftotext too, so it is not only a rendering loss. The YAML decision that treats it as a slip matches the tests, and the explainer states it.

### toychem

- Energy model has no conjugation. In eht(), pi overlaps are placed only across double and triple bonds, never across the single bond between two sp2 atoms. Ethene and butadiene therefore get identical frontier levels (HOMO -13.75 eV, LUMO -6.16 eV, from my run), and so do the Diels-Alder barriers for butadiene+ethene and butadiene+butadiene (both 175.1 kcal/mol). This undercuts the frontier-orbital reactivity model for dienes. The 2003 paper mentions extra rules for resonance (lone pairs interacting with adjacent pi systems), and those are not implemented either.
- Reaction energies are not proportional to the published ones. Butadiene+ethene gives -112.1 and two butadienes -92.4 kcal/mol, against -16.76 and -8.15 in the 2005 Table 2. Scaling by the TAE factor (about 0.38) would give roughly -43 and -35, and the ratio between the two is 1.2 against 2.1 published. The decision saying 'relative energetics - orderings, reaction energies and barriers - carry over' overstates this: only the sign and the ordering carry over. The test only checks the ratio with rel=0.6.
- The YAML phenomenon says formose is not small-world (ECAL Table 1). Chemart's 11-species formose network comes out small_world=True (C 0.47 vs 0.29 random, L 2.07 vs 2.25), so this is not reproduced, and no test covers it. Chemart's formose network has 11 species against the published 48.
- barrier_cutoff is all-or-nothing on the default network. At a cutoff of 131.5 or more you get the full 11 species; at 131 or less you get only the 2 seeds, because the first keto-enol step has Ea 131.1. It does not reproduce the gradual growth of the ECAL figure 3 series. On the Diels-Alder network, once max_species truncates it, the cutoff is not even monotone. With max_species=120: cutoff 133 gives 181 reactions, while cutoff 0 (no gate) gives 139.
- The YAML R.notes says 'Which of the possible channels survives is decided energetically, not combinatorially'. In the default rewrite_mode='all' with barrier_cutoff=0, every channel survives, and energy only gates channels when a cutoff is set. I did not edit this because it is wording and judgement rather than a verifiable number.
- rewrite_mode='priority' keeps the most exothermic channel. The 2003 paper only says the rule with the highest 'priority value' is chosen, and for regioselectivity it says 'the rewrite with the smallest ∆E value is chosen' (∆E being the reactivity index, not the reaction energy). So 'most exothermic' is a Chemart interpretation, but the YAML presents it as the paper's mode.
- The keto-enol barriers are each molecule's own HOMO-LUMO gap. The forward and reverse barriers (131.1 and 138.0) therefore do not differ by the reaction energy (+25.9), so the barriers are thermodynamically inconsistent. The decisions acknowledge this in general terms only.

### typogenetics

- Book bibliography [785] calls Snare's thesis a Master's thesis. Its title page says 'Thesis for Bachelor of Science (Computer Science) Honours'. The YAML already has this right (BSc (Hons)).
- Morris [597] page numbers: the book gives pp. 341-368, but Snare's bibliography gives pp. 369-395, and Snare cites Morris pp. 379-382, which fits his range. Morris's PhD thesis [598] is dated 1989 in the book and 1988 in Snare. The explainer's Further reading notes both discrepancies.
- Snare fig. 3.2(b) labels CGATTCGAATCG's enzyme 'cop-swi-rpu-inc-swi-cop'. Under Hofstadter's table (Snare's own table 2.1a), GA codes ina, which gives cop-swi-rpu-ina-swi-cop, as Chemart and the YAML have it. This looks like a typo in the thesis; it does not change the outcome because that amino acid never runs.
- The earlier job's files wb_Typogenetics_paper_ECAL_old.pdf and wb_Typogenetics_paper_THEOCHEM_final.pdf (Kvasnicka et al.) and morris.pdf are only about 5 KB each: failed downloads, not papers. This matches the YAML decision that these works were not accessible.
- Book table 10.1 prints Hofstadter's table on both sides even though the caption says the right-hand one is Varetto's. The YAML already notes this erratum.

## Archive

### aevol

- Crash with default parameters: a genome shorter than about 11 bases makes decoding raise 'ValueError: window shape cannot be larger than input array shape' (from _windows / sliding_window_view in terminators or rbs_sites). min_genome_length defaults to 1, so a large deletion can produce such a genome. Reproduce with generate_network('aevol', seed=2, grid_width=4, grid_height=4, generations=60, max_genome_length=10000); direct check: A.transcribe/translate on np.zeros(10) fails, while np.zeros(12) works. The page tells users to set min_genome_length=100 as a workaround. The slow test only uses seeds 1 and 3, so it never hits this.
- The params' range fields contradict Knibbe et al. 2007 (params were not edited, per the rules): point_mutation_rate says 'Knibbe et al. 2007 scan 10^-6 to 10^-3', but the paper uses 5e-6 to 2e-4. genome_length says 'Knibbe et al. also evolve from 1000 to 20000'; the paper only says equilibria did not depend on initial size ('data not shown'), so this cannot be verified. generations says 'the published runs go to 10^5-10^7 generations'; the 2007 runs were 20,000 generations.
- The 2007 model differs from the ported 9.4.0 source: it uses a 28-bp promoter consensus instead of 22, and exponential-ranking selection over N=1000 instead of local 3x3 fitness-proportionate selection. The YAML and decisions do not mention this, and the decision saying the papers could not be fetched is now out of date.
- The decision says the default run takes 'about a second'. It took 3.0 s here.
- The generated 'fidelity' line says 'built from the original papers listed under References', but this entry was built from the aevol source code. That boilerplate is the generator's, not this entry's.
- Individual gives fitness 0.0 when the phenotype area is zero, not exp(-k*error). This may or may not match upstream compute_fitness; I did not check it against the source.
- The book (18.1.2) describes Parsons et al. 2011 as being about the evolvability of 'genetic regulatory networks'. The abstract summary I could find talks about genome structure and evolvability, not regulatory networks. The paper itself could not be read.

### arn

- tools/gen_catalog_pages.py esc(): `str(text or "")` turns a default of 0 (also 0.0/False) into an empty string, so the Parameters table shows link_threshold's default as `` instead of `0`. Likely affects other entries with zero defaults.
- catalog/chemistries/cpm-grn-evodevo.yaml decisions say Chavoya & Duthen [172] is 'covered by the arn entry', but arn does not implement it (arn.yaml's own decisions say the [172] morphogenesis extension is not generated).
- Banzhaf (2003) Table 1 gene counts (37 genes at 10,000 bits, 409 at 100,000) are not what Chemart's decoder gives for random genomes (27 and 288 on average over 20 seeds). Table 1 is close to the raw 0.39% promoter rate. Chemart follows the later papers' rule, where overlapping promoters count once and a gene needs room for both sites (it matches Kuo & Banzhaf's 340-440 at 131,072 bits). This is not a bug, but the YAML sources cite Table 1 without noting the mismatch.
- tests/chemistries/test_arn.py test_equal_start_settles_with_a_few_dominant_proteins checks 'point attractor reached' only as a change of <1e-2 between t=1.8e7 and 2e7, on 2,048-bit genomes. At the default length, seed 9 (13 genes) oscillates without settling (period about 3.7e7, still steady at t=4e8), so the test does not show that dynamics generally settle. It also does not check the damped oscillations named in phenomena.
- Phenomena 'small-world and scale-free topologies', 'network motifs matching natural GRNs', 'heterochrony' and 'evolvable to target dynamics' have no tests. They are listed as literature phenomena, and the page says plainly that they are not reproduced.
- Chemart has no ODE integrator in the library (only tests/chemistries/odes.py scaffolding), so the page's dynamics example integrates with SciPy by hand.

### avida

- Phenomenon 3 in the YAML and the comment in test_default_run_contains_replication_and_overwrite_events both say births overwrite organisms once the lattice is full. In fact overwrites start as soon as a local 8-cell neighbourhood is full. In the default seed-1 run the first overwrite comes at update 81, with about 28 of 144 cells occupied, and 122 of 279 births are overwrites while the world never fills (135/144 at the end). The test only checks that the run ends with more than 100 organisms.
- Book bibliography [629] gives pages 3-34; the YAML source says 3-35. Minor.
- Cooper & Ofria (the ecosystems paper) is dated 2002 in book ref [198] and 2003 in Ofria et al. 2009's bibliography. The explainer follows the book.
- Ofria et al. 2009 (sec. 4.1) says the EQU runs lasted about 17,000 generations and that at least 19 coordinated instructions are needed for EQU. Lenski et al. 2003 say 15,873 ancestral generations, and only that their shortest hand-written EQU program had 19 instructions, 'not proven' to be the shortest. The explainer uses Lenski's figures.
- Lenski et al. 2003 (Avida 1.6) give SIPs in proportion to genome length x merit. Chemart follows Avida 2's BASE_MERIT_METHOD 4, which uses the smallest of genome length, copied size and executed size. This is a real modelling difference from the paper, but the YAML decisions describe it only indirectly. The explainer states it.
- Book table 10.4 misprints (if-n-eau, nopC). This is already recorded in the YAML decisions.

### bnc-cell

- The YAML decision blames only the book for the import example ('12321000000 transports 1-2-3=2-1'). The paper's Methods has the same inconsistency: it writes '123321000000 specifies that molecule 1-2-3=2-1 is transported'. The code implies 1-2-3=3-2-1, but the molecule the paper names breaks the bond rule. The book copied the paper's molecule and dropped one digit of the code. Chemart's reading (1-2-3=3-2-1) is reasonable, but the decision should say the error is already in the paper.
- The paper disagrees with itself on the gene layout. Table S1 and the Results text give start, expression, type, specificity. The Figure 11 caption and labels give start, type, expression, specificity. Chemart does not decode genomes, so this has no effect on it.
- Figure 11's second import protein (1-2-2-3=3-2-2-2-2-2-3=2, 12 atoms) targets a metabolite, not a precursor. Random cell mode only lets importers target precursors, while an explicit proteins list accepts any target. Both are defensible, but the decisions do not mention it.
- The paper says the initial genome is 2,000 bp. Its Methods describe two 1,000-bp chromosomes, the first holding 880 bp of genes padded to 1,000. The two accounts agree, but the YAML only mentions the 1,000 bp chromosome.
- With the defaults, random cell mode builds cells whose parts usually do not connect: with seed=1 the only imported precursor is used by none of the three enzymes. This follows from the uniform draws, a Chemart choice, and is documented on the page. The default cell is not a working metabolism.

### corewar

- Decision 5 (can't edit) says 'the book says as much for assembler automata'. The book says it only of Coreworld (10.6.3).
- Book bibliography [432] dates Core War Guidelines to 1994 and links to guide2red.txt. The Guidelines are Jones & Dewdney, March 1984 (cwg.txt reproduction). The guide2red.txt file downloaded by the earlier job is a 404 page.
- Book bibliography [894] calls the venue 'Ninth International Conference ... (Alife XI)'. The paper itself says ALife IX, 2004.
- The book's 'trace of MOV 1 0 instructions' is a typo for MOV 0 1. The YAML decisions already record this.
- The imp-avalanche test only checks that the task series is monotone and reaches the limit. It does not check the growth rate the phenomenon used to claim; that claim was wrong and is now fixed.
- Dewdney's 1984 figure for Dwarf v. Dwarf (30% / 30% / 40% draws) is not tested, and the rules it assumes (zero counts as DAT) are not implemented.

### coreworld

- YAML decision says the defaults make a run take 'about half a second', but the default run (seed=1) took about 2.1 s here. The desert setting took 0.5 s. Decisions are outside what I may edit.
- In the jungle Table 2 sets resource_influx (0.5) equal to resource_max (0.5), so every cell is refilled to full after each update. extras.analysis.series.mean_resource is sampled after the refill, so it is always exactly 0.5 in the default run and carries no information. This is not a bug, but the series only tells you something in a desert.
- Book bibliography [894] (Vowk, Wait & Schmidt) contradicts itself: 'Ninth International Conference on the Simulation and Synthesis of Living Systems (Alife XI), 2004'. The Ninth conference was Alife IX (2004). The page cites it without the numeral.
- No test covers the paper's composition floods (phase transitions), SPL-fields, SPL-MOV organisms or the four epochs. The desert/jungle test only compares execution counts and split counts (jungle > 1.5x desert), not the structures the phenomena describe.
- The YAML sources credit the 1989 preprint to 'Rasmussen, Feldberg, Hindsholm & Knudsen' while origin uses the book's order (Rasmussen, Knudsen, Feldberg & Hindsholm). The two differ because the preprint and the Physica D paper list the authors in different orders. Not an error, but it can confuse readers.

### cpm-grn-evodevo

- The YAML decision 'Small defaults: a 32 x 32 lattice ...' contradicts the params, where grid defaults to 38. Decisions are outside what I was allowed to edit, so it is not fixed.
- The YAML scope decision says Chavoya & Duthen [172] is 'covered by the arn entry'. The arn entry implements only Banzhaf's base ARN. Nothing in arn.yaml mentions Chavoya or pattern formation, so the ARN pattern-formation extension is not implemented anywhere in the catalog.
- In the default run, cells are far below their target volume. After each scheduled cleavage both daughters get the reference target V = 50, but 10 steps per stage is too little time for them to regrow. The 4 survivors hold 12 to 18 sites each, and 4 of the 8 cells die. With steps_per_stage=40 on a 60 x 60 grid, the cells hold 43 to 48 sites. The papers do not say how cleavage divides the zygote's volume, so it is unclear whether resetting the target to V after a scheduled cleavage matches Hogeweg.
- Growth-triggered division, one of the model's central mechanisms, never happened in any run I made (the default, 8 seeds with final_steps=200, 3 seeds with population=6 and generations=4, and a 30 s run on a 60 x 60 grid). The largest target volume reached was 53, against the 100 needed. No test checks it: test_morphogenesis_is_a_side_effect_of_differentiation asserts seen['divisions'] > 0, which the scheduled cleavages alone always satisfy.
- The docstring of test_morphogenesis_is_a_side_effect_of_differentiation claims to cover 'the phenomena of Hogeweg (2000) sec. 3.2' (the morphogenetic mechanisms). It only checks that divisions, deaths and differentiations occur and that more than two expression patterns appear. It detects no engulfing, budding, intercalation or meristem.
- `analysis['differentiations']` counts every step of a period-2 cycle as a differentiation event (T3<->T4 x70/x69 in the default run). A reader could take it for real cell-fate changes.

### dimerization

- tests/chemistries/test_dimerization.py does not exist, but the generated page footer links to it (gen_catalog_pages.py line ~412 always writes tests/chemistries/test_<id>.py). The only chemistry-level test is test_dimerization_reaches_equilibrium_constant in tests/chemistries/test_w1_dynamics.py.
- The YAML phenomenon 'SSA fluctuates about the ODE equilibrium' is not tested anywhere, and Chemart has no stochastic simulator; only k_to_c is unit-tested (tests/test_core.py::test_k_to_c).
- The YAML has no `sources` list, so the References section shows only the book section. It also has no `refs` links (e.g. PyCellChemistry) even though reference_impl names Dimer.py.

### energy-gated-collision

- The YAML sets `constructive: true`, so the page's at-a-glance table says 'yes — the species set grows at run time'. That is false for this entry: the species set is fixed by the system passed in (S is `explicit`), and no new species ever appear. I did not change it because `constructive` is not one of the fields I was allowed to edit.
- The YAML `provides` list includes `thermodynamic-consistency`, but the generated network's `net.provides` does not (the summary prints energies, initial-state, rate-constants, rate-law, stoichiometry, topology). The test allows only that smaller set.
- Phenomenon 6 has two halves: cost per collision does not depend on the number of reactions, and the acceptance rate is exp(-Ea/RT). No test measures the first half. The tests only assert that the default run is more than 50% elastic.
- The parameter `A` never enters `run()`. It only labels the rate dicts and scales `k_observed`, so every reaction shares one prefactor and the simulation cannot give two reactions different values of A. The YAML's meaning text ('It scales every rate and cancels out of the acceptance measurement') is consistent with this but does not state the limitation.
- When several rules share a reactant pair, each collision draws one rule uniformly. Each rule's effective rate per pair collision is therefore (1/n_rules)·exp(-Ea/RT), which amounts to dividing its prefactor by the number of rules. This is documented as a decision but is not flagged as a change to the kinetics.
- The YAML `notes` call the algorithm 'Worth shipping as a Chemart reactor backend'. It is not exposed anywhere outside this entry: the grep found no other use in the code.

### evolve-series

- tests/chemistries/test_evolve_series.py::test_resistance_to_phenotypic_change_evolves only checks that mean sensitivity falls from the founders' 3. A mutated sensitivity locus is drawn uniformly from 1..4 (mean 2.5), so the fall happens without selection, and starting from 1 or 2 the sensitivity rises. The test does not check the EVOLVE II claim it cites.
- YAML decision on initial_mineral says '(about 216 matter units for 100 cells)'. The real total is 100 x 3 + 8 founders x 2 = 316 (analysis.matter_total = 316). Not fixed: decisions are outside what I may edit.
- The same decision says that at initial_mineral 3 'scavengers appear in every seed' of six. With the current code, seed 4 at defaults produces no scavenger in 120 steps (scavenger_episodes = []).
- The book's headline phenomenon (scavengers go extinct and reappear) never happens in the default 120-step run for seeds 1-6. It needs steps=600 and shows up in 3 of seeds 0-6 (1, 4, 5); in seeds 1 and 5 the gap is only 3 and 6 steps. Only one slow test checks it, on seed 4.
- Over 600 steps enzyme matches decay heavily (final fix efficiency 0.22-0.65 in seeds 0-6), and the 'none' guild (no fix and no respire enzyme; lives on light and builds offspring from others' organic matter) grows to 16-69 of the survivors. This emergent free-rider class is not described in the YAML phenomena.
- Most of the 'scavengers' counted in history are mixotrophs (fix and respire). In the default seed-1 run the final guilds contain no pure scavenger.
- net.summary() lists provides without 'sequence-structure-function', which the YAML's provides includes.
- Book bibliography [195] errors (author 'M. M. Pattee', pages 293-409) are confirmed by OpenAlex: it is H. H. Pattee, pp. 393-409. They are already recorded in the decisions.
- OpenAlex dates Brewster & Conrad's CEC paper (doi:10.1109/CEC.1999.781957) to 2003. The book and the YAML say 1999 (CEC 1999), which is the conference year.

### french-flag

- No tests/chemistries/test_french_flag.py exists, yet the generated page footer links it. The entry's only test is test_french_flag_three_segments in tests/chemistries/test_w1_structure.py, which checks that 9 cells give 3 blue, 3 white, 3 red.
- The second YAML decision mentions 'n_segments is implied by the number of thresholds', but there is no n_segments parameter. The wording is stale.
- The YAML's A lists reactor [lattice-2d, continuous-space], but the generator runs no reactor. It only emits one rate-less U_x -> fate_x reaction per cell, and the default is 1D. Neither the book nor the code supports continuous-space.
- The default threshold 0.667 is not exactly 2/3. With n_cells=4 the second cell (c=0.6667) comes out white, giving BWWR instead of a symmetric split. This is harmless but surprising.
- Fate names switch from blue/white/red to type0..typeK whenever the number of thresholds is not 2. This is undocumented except in the code (now explained on the page).
- Size invariance is built in: the gradient is normalised to x in [0,1] for any n_cells, so the stripe proportions never change. This is not the regulative size invariance of Wolpert's problem, and the page says so.
- The YAML has no phenomena list, so there is nothing for the tests to cover beyond the three-stripe check.

### hbcb-psd

- The ordering 'BM-BM links are the weakest bonds, deeper bonds stronger' appears in the YAML decisions, phenomena and parameter text as 'the model', but neither the abstract nor the book says which way the hierarchy is ordered. The abstract only says biomolecules are arranged 'according to the energy strength of their covalent bonds'. The direction is a Chemart inference. The page labels it as Chemart's reading, but the YAML still presents it as published.
- The flagship result (BM is the cheapest class to rebuild from) is true by construction: the route back from any deeper class includes every BM step plus more, for any positive bond energies. So test_decomposition_to_bm_is_the_cheapest_route and test_bm_is_optimal_for_every_hierarchy check arithmetic, not an emergent property. The bond-energy ordering constraint has no effect on which class wins. The page says this plainly.
- A lists reactors lattice-2d and well-stirred-multiset, but the network has no rates (all None), and Chemart supplies nothing that could run it on either reactor. The listed reactors describe the paper's ecosystem, not anything Chemart can simulate.
- The YAML says Sayama's 1998 SIVA-2r belongs to 'the same group'. Sayama was at the University of Tokyo and thanks Prof. Oohashi in the acknowledgements, and a related SIVA paper (Oohashi et al. 2001) is by the Oohashi group, so this is loosely right. Note that SIVA compares self-decomposition against none, not decomposition depths, so it does not bear on the 2009 HBCB result.
- The previous job's siva2r.txt is 8 bytes because the Sayama PDF is a scan (pdftotext gives nothing). I checked the SIVA facts by reading the rendered page images (pp. 1-3, 6-8).

### isologous-diversification

- The same misattribution of the 1998 paper to Kaneko & Yomo (it is by Furusawa & Kaneko) is still in the `enzyme` param meaning ('Kaneko & Yomo 1998 eq. 2'). Params were off limits, so it was not changed. It is also in the generator docstrings (composition(), classify() in chemart/chemistries/isologous_diversification.py) and in the test module docstring.
- The network screen tests the single cell from different random starting concentrations than the actual run uses. An accepted network can therefore still never divide: seeds 4 and 6 at defaults, and seed 0 with network_attempts=16. For example, seed 1's network started from default_rng(1) concentrations divides once, then settles on X8, which has no path to DF, and stops.
- With default settings the screen often fails. Over seeds 0-9 it found an accepted network for only 1, 4, 6 and 9; the rest fall back to the best rejected candidate (accepted=False). Only 6 of 10 seeds reach 8 cells.
- No run I tried went past stage 1. Default settings with max_cells=32 (seeds 0 and 3 reach 32 cells, seed 2 reaches 16) all stayed in stage 1 at t=300. So did the 1997 Fig. 18 parameter set (V=1000, nutrient=10, connections=2, max_cells=32, t_max=1000, seed 1), where the paper already sees differentiation by 32 cells. The decision's claim that 'raising max_cells, t_max and division_threshold to the paper's values runs them' is not demonstrated anywhere.
- The default composition of seed 1 is about 89% X7 and 11% X1, with the activity sum near 727, far above x_M=10. The medium's source is exhausted (X0 about 0). This looks more like the paper's winner-takes-all or tumour-like regime than a rich oscillation. The default V=100 comes from the different 1998 model and is combined with 1997 kinetics.
- enzyme='quadratic' only swaps the catalytic term inside the 1997 model. The rest of the 1998 model is not implemented: diffusion-only exchange, penetrable vs impenetrable chemicals, volume growth with normalisation (sum x = 1), division at volume doubling, and (1±ε) splitting with ε about 1e-6. The param text could suggest otherwise.
- In stage 1, analysis['recursivity'] and 'inherited_fraction' are non-trivial (e.g. 0.19 and 0.43 for seed 1) although the cells are identical: the composition still drifts between early generations. The numbers are easy to misread.
- In the Fig. 18 run, average_spread came back exactly 0.0 while snapshot_spread was 0.00035. The likely cause is that only one or a few cells count as 'mature' just after the last division. Not investigated.
- Book §18.5 credits isologous diversification to Furusawa and Kaneko [301], the spatial model. The originators are Kaneko & Yomo (1994/1997), and Chemart implements their non-spatial 1997 model, not the one the book describes and pictures (Figs. 18.10, 18.11, 18.13).
- The 1997 paper's text says 'When m = ℓ, the reaction is regarded as autocatalytic'. That is presumably a typo for j = ℓ (catalyst equals product), which is what Chemart implements.

### logistic-chemistry

- YAML decision says 'r and x0 are not given (defaults 1 and 0.1)'. That is true of the book, but the reference_impl PyCellChemistry Logistic.py (the script that made Figure 7.2) does set them: r = 1, K = 1, x0 = 0.01, NAV = 100 (one starting molecule), run to t = 50. So Chemart's default x0 = 0.1 does not match Figure 7.2's source. I did not edit it because decisions and params are outside what I may change. The page says so.
- The generated page footer lists tests/chemistries/test_logistic_chemistry.py, but that file does not exist. The only test is test_logistic_plateaus_at_carrying_capacity in tests/chemistries/test_w1_dynamics.py. The footer is hard-coded by gen_catalog_pages.py.
- The second phenomenon (stochastic runs fluctuate around K and overshoot it) has no test. Chemart has no SSA; only k_to_c is tested, in tests/test_core.py.
- The YAML has no `sources` list, even though Logistic.py is the reference implementation and holds the Figure 7.2 settings.
- The catalog lists 'catalysts' under provides for a network whose only 'catalysis' is X's own autocatalytic replication. Possibly a mislabel. I did not change it.

### lotka-volterra

- The generated page footer points to tests/chemistries/test_lotka_volterra.py, which does not exist. The only test is test_lotka_volterra_matches_book_reactions_and_conserves_invariant in tests/chemistries/test_w1_dynamics.py. The generator (tools/gen_catalog_pages.py line 412) builds that path by assumption, so other W1 chemistries (brusselator, logistic, and others) probably have the same dead pointer.
- Nothing tests the stochastic (SSA) phenomenon, 'erratic oscillations and extinction under SSA in small volumes'. Chemart has no built-in SSA, so the phenomenon was only checked with an ad-hoc Gillespie script. With the book's figure 7.5 settings (omega=20), all 10 seeds went extinct within 6 time units, before the first 9.3-unit cycle was complete. The book's figure presumably shows a run that lasts longer; I could not see the figure.
- The YAML's provides list includes 'catalysts', but the default network has no catalytic reaction. Catalysts appear only when an off-diagonal b_ij>0 is not paired with b_ji=-b_ij.
- The YAML's refs list omits bibliography entries the section also relies on: [5], [292], [712] and [834] (cyclic chains and mutation) and [329] and [926] (the SSA conversion in the figure 7.5 caption). This is an omission, not an error, so I did not change it.
- The book states the fixed point as 'x = c/b and y = a/b' (shorthand for kc/kb, ka/kb). The YAML phenomenon (kc/kb, ka g0/kb) is correct and more complete.

### n-economy

- YAML `origin` reads "Straatman, Banzhaf et al.; after Hanel's natural number economy". The book (§20.3, ref [377]) credits the natural number economy web page to J. Herriott and B. Sawmill (redfish.com); no Hanel appears. `origin` is not one of the fields I may edit, so the page header still shows it. The explainer gives the book's attribution.
- YAML `constructive: true` makes the at-a-glance table say "the species set grows at run time". Chemart's generator returns a fixed list: the product_set goods and the technology reactions, 7 species and 5 reactions by default, with no growth and no randomness. The book's product set P is fixed as well.
- No tests/chemistries/test_n_economy.py exists, although the page footer links to it. The only test is test_n_economy_default_is_book_eq_20_6 in tests/chemistries/test_w1_structure.py. It checks the five reactions and the factorisation of 260 only; the parameter validation has no test.
- The YAML has no `sources` entry, so the References section lists only the book. I added full citations for [816] and [377] under Further reading instead.
- The seed is ignored: the generator never uses rng, and net.extras is empty.
- Good 3 (money) is in the default product_set but in none of the reactions, so it is an isolated species.

### naming-game-ac

- Default behaviour does not match the claimed phenomenon. With kappa1 = kappa2 (the defaults, all rates 1), every name's total S_j + C_j is exactly conserved: the first three schemes only turn C_j into S_j or back, and the two mismatch reactions share reactants and move material in opposite directions, so d(S_j+C_j)/dt = (kappa1 - kappa2)·M·(S_j·C_tot - C_j·S_tot). No shared lexicon forms. The YAML's intuition and phenomena ('lateral inhibition ... converging to a shared lexicon') therefore hold only for kappa1 > kappa2 with M held constant. Checked by runs (4 s).
- M is consumed by eqs. 16.2-16.4 (as the book prints them) and never produced. From M=1 it is gone by t≈5 and all dynamics stop before any competition. The generator has no option to hold M fixed (buffered). The YAML gives no initial state and no extras.
- The kappa labels look swapped relative to the book. The book's prose says 16.5 and 16.6 'eliminate mismatching words (resp. adaptors)', which pairs κ1 (16.5) with replacing the word. Chemart's kappa1 replaces the adaptor (M + S_j + C_k -> M + S_j + C_j), and so do the YAML param meanings. The dynamics are symmetric under the swap, but the labels differ. Params were not changed, as instructed.
- There is no tests/chemistries/test_naming_game_ac.py, yet the generated page's footer links to it. Only the generic contract tests in tests/test_contract.py cover this entry, so nothing checks the lexicon-convergence phenomenon.
- The original paper (De Beule, Hovig & Benson 2011, Biosemiotics) is closed access (Unpaywall: closed), so the reconstruction of eqs. 16.5-16.6 and the conditions under which the authors saw convergence could not be checked.

### nk-landscape

- There is no tests/chemistries/test_nk_landscape.py. The NK tests live in tests/chemistries/test_w2_papers.py, but the generated page footer links to the missing test_nk_landscape.py.
- test_nk_fitness_follows_book_index_convention does not check the index convention. It only checks that there are 8 fitness values in [0,1) and that mu=0 leaves only copy reactions. No hand-computed M(i,j) lookup is compared.
- The epistatic table M(i,j) is not stored in extras.analysis, so users cannot inspect per-gene contributions. The page recovers it by redrawing default_rng(seed), which only works for the adjacent topology.
- The YAML origin is 'Kauffman, 1993', but Altenberg (1997) dates the model to Kauffman & Levin 1987 and Kauffman 1989. origin is not one of the fields I may edit.
- The YAML has no phenomena and no sources field.
- Book §18.4.1 gives 0 ≤ j ≤ 2^(K+1). It should be < 2^(K+1), since there are 2^(K+1) entries indexed from 0. The text also refers to 'Figure 18.4.1' where it means Figure 18.9.
- Book §14.2 (around line 13682 of the extract) describes 'NK fitness landscapes formed by random networks of N Boolean logic elements' and cites [445]. This conflates the NK landscape with random Boolean networks (NK networks).
- In a deterministic ODE every mutant is present at once, so the replicator-mutator network always converges to the global optimum's quasispecies. It cannot show populations stuck on local peaks, which is the main subject of the NK literature. The page says this.

### organization-computing

- No file tests/chemistries/test_organization_computing.py exists, but the generated page footer names it. The only tests are test_xor_gate_follows_truth_table and test_maximal_independent_set_rules in tests/chemistries/test_w1_structure.py. They check reaction structure only, not organization counts.
- Chemart puts the XOR inputs into initial_state as one molecule each. The papers (Matsumaru et al. 2005/2007, 2011) supply inputs as constant inflow reactions (∅ -> a1). The generated network therefore has 15 organizations (the gate with no input specified), not the single one the paper gets with both inputs. The page explains this, and the organization-collapse result is only reproduced once the user adds inflows by hand.
- The YAML declares A: reactor [ode, ssa], but every reaction has rate None, so neither reactor can run as generated. The page says Chemart provides structure only.
- The YAML has no `sources` field, so the page has no 'Sources used for the implementation' list. The papers used were the IJUC manuscript (MCSD2005), BIONETICS 2006 and the 2011 chapter preprint.
- The YAML intuition says 'read off a lattice of self-maintaining sets', and the notes say 'ORGANISATION LATTICE'. Matsumaru & Dittrich (2006) state explicitly that the MIS organizations 'do not form a lattice, because there is not a unique largest organization'. This was not edited because it is terminology rather than a clear error, and neither field is printed on the page.
- Book 17.3.3 cites the MIS algorithm first to [551, 553], then to [554], and calls it the 'maximum independent set problem' once (it is maximal). It also says a 'high concentration' of s1_j marks membership, while the section's framing is presence/absence.
- Citation discrepancy for book ref [552]: the book gives Matsumaru, Kreyssig & Dittrich, pp. 207–220. The authors' preprint header says 'Published as: P. Dittrich, P. Kreyssig ... Volume 1, Part 1, 67-78'. Not resolved.
- The YAML phenomena list 'flip-flops, oscillators, NAND chains', but the generator implements only the XOR gate and the MIS program.
