# One consolidated Universal Model: first migration outline

This is an editorial plan for one self-contained manuscript superseding the current TAUM/LG pair. It is not a new
companion paper and supplies no premises needed by either current manuscript. Neither manuscript is edited by this
proposal.

## Source and review scope

The [working manuscript](../manuscripts/consolidated/Universal_Model.tex) opens with the problem of perception and a table of explanatory tasks. Its first five sections cover preparation and its dynamical reach, polar supports and matched observables, neutral graphs and relative phase, the ground correspondence, and the determination algebra. Sections 6–8 treat asymmetric nonduality across two moments, cohomology and the index, and effective dynamics in the pair and orbital models. Section 9 develops the passive response of a self-adjoint system observed through the same port that prepares it. Section 10 follows a prepared plane as the apparatus turns and the system evolves, relating the geometry of that family to the measured response of §8. Section 11 asks which vectors a supplied compact symmetry fixes and when a port measures exactly the dynamics on that subspace. Section 12 compares four maps with the same algebraic relations and different spectra and responses. Section 13 identifies the currents that prepare and read the pair and relates them to angular momentum and the two quadrupoles on the bare orbital shell. The numbering can change as other result families are integrated.

- Both manuscript sources: Theory commit `8e9fde9`.
- Initial statement inventory and six proof annotations: `58b8354`.
- The [migration ledger](inventory/migration.csv) compares source clauses with the consolidated introduction and §§1–13. It distinguishes a statement’s source kind from its role in the draft and records passage locators, evidence and outstanding clauses. The all-channel strengthening of R5 is attributed to the consolidated manuscript (`UM`). The [ledger guide](inventory/README.md#migration-into-the-consolidated-manuscript) explains the fields.
- Supplements cover the pair and orbital equations used in §§1–3, 8 and 13. The [source reading queue](inventory/heading_queue.csv) identifies further passages without a direct ledger entry. Read the surrounding derivations when migrating a partially tracked section.
- The inventory now includes all used, declared theorem-like environments: 106 TAUM entries and 175 LG entries,
  including LG's 17 model choices and 29 doctrinal identifications.
- This outline maps all 50 top-level sections/appendices (16 TAUM, 34 LG), plus three major embedded subsections, in
  [SOURCE_MAP.csv](SOURCE_MAP.csv). Split destinations are deliberate; a source section need not migrate intact.
- This is a section/result-family placement review, **not a completed dependency audit of all 281 statements**.
  The six existing annotations were checked against their source passages and are retained unchanged. The other
  dependency descriptions below identify assumptions that must be checked during migration; they are not newly
  verified minimal-hypothesis theorems.

## Decisions already supported by the discussion

1. One canonical working manuscript will contain the common mathematics, the substantive cognitive/Buddhist
   development, and the physical research program. It will not be maintained in parallel Buddhist-on/Buddhist-off
   editions.
2. The final TeX must stand alone. Required proofs and technical appendices remain in that document. These planning
   files and the verification code are not substitutes for proofs in the manuscript.
3. Use Ψ and Φ for the two sector labels and A for the typed channel from Ψ to Φ throughout the consolidated manuscript. These labels specify mathematical roles. Part V gives them their physical and phenomenal reading. The inventory retains the identifiers of the source manuscripts.
4. State the long-term physical-unification ambition openly. Distinguish that aim from the architecture, conditional
   theorems, and realizations presently established. No Standard Model, physical gauge-ghost, or microscopic GTD
   derivation is added by consolidation.
5. Preserve the full cognitive and Buddhist argument. Generalizing a mathematical lemma must not delete its original
   philosophical use, sources, or motivating question.
6. Reuse a result through internal references once it has a canonical home. Retain a second version only when its
   hypotheses, domain, or conclusion genuinely differ.
7. Keep the positive Hilbert probability structure distinct from the Krein form throughout. Keep typed channels,
   graded extensions, polar maps, and the different projection constructions distinct.
8. Use `\noindent` for the abstract in the consolidated TeX.

## Proposed organization

The part and subsection codes below are planning identifiers, not final TeX labels or a proposed page limit.

### I. Problem, origin, and program

Begin with the problem that grew out of Dharmakīrti's account of perception. In the cognitive reading, A relates physical states to their presentation. The grasper–grasped distinction arises within the presentation. Explain why presentation, the readout of selected content, and the physical conditions of that readout are separate questions. State the broader ambition toward physical unification and identify the missing selection principles.

Introduce the five optional realization hypotheses in plain terms. End the opening with a table locating the tasks of appearance, physical conditions, determination, reflexivity, error and transformation, with the status of each proposed answer. Part V develops the cognitive argument in its own sequence and gives the reasons for its choices.

### II. Common mathematical architecture

**II.1 Carrier, channel, adjoints, and returns.** One treatment of the two complex Hilbert sectors, positive and
indefinite products, the closed typed channel, natural domains, graded orientations, odd-adjoint conversion,
canonical Hilbert/Krein branches, and positive partner returns. Include LG's correction to an untyped nilpotent
isometry here. Move the distinguished-reference subsection out of the bare-channel assumptions.

**II.2 Polar supports, correspondence, and observables.** Polar decomposition and rotations; active-support
isometric graphs, neutral geometry, transition-probability and effect transport; paired observable algebra; sector
factorization and phase blindness. Separate the maximal constructed algebra from a physically selected observable
algebra. Explain that transport of Born probabilities is not a derivation of the Born rule. Ground-inclusive
specializations belong in III.1.

Sections 2 and 3 transport spectral measurements between the active supports, characterize the block-diagonal operators that preserve a phase graph, and distinguish the sector compressions of a graph state from its partial trace. The source material from LG §§16–19 and TAUM's polar section is allocated as follows.

| Source material | Destination |
| --- | --- |
| Polar reflection and rotation, the neutral graph and its projectivization, sector factorization, phase blindness | II.2, drafted in `sec:neutral-graphs-phase`. The odd-branch factorizations and normalized transitions remain with their operator domains in the common architecture. |
| `choice:active-physical-observable-algebra`, `choice:correspondence-versus-generation` | II.2 states the two selections. III.3 treats generators and V applies the choices to cognition. |
| `cor:ground-inclusive-BRLC` | III.1, drafted in `sec:ground-correspondence` with the reference and bridge hypotheses |
| `id:BRLC-doctrinal`, `id:one-event-two-compressions` | V, after the mathematical graph and paired representation |
| `q:operational-structural-phase` | Stated at the end of the sector-factorization discussion in II.2; carry to VI as the open question of operational access to the structural phase operators. |

**II.3 Kernels, cohomology, and index.** Closure-aware Hodge decomposition, partner spectral pairing, reduced
cohomology, bounded transforms, Fredholm index, and conditional heat-supertrace statements. State the roles of
closed range, Fredholmness, and trace-class assumptions where used. Introduce LG's two independent map hypotheses here: R4 is injectivity (`ker A = 0`) and R5 is strict final defect (`ker A* ≠ 0`). Each can be stated before choosing a reference vector or a Hamiltonian. V motivates and interprets their joint adoption in the cognitive realization.

The working section `sec:cohomology-index` now carries TAUM’s two-term complex, reduced harmonic representatives, spectral pairing, closed-range criterion, graded Fredholm and heat identities, and LG’s single-channel index interpretation and shift. A rank-one variation of the shift keeps the ground bridge and moves the missing line. The extension count distinguishes the minimum new input dimension for surjectivity from the input dimension that cancels the index.

| Source result | Location in the working draft |
| --- | --- |
| TAUM `def:cohomology`, `thm:Hodge-decomposition`, `cor:harmonic-representatives` | `eq:index-cohomology` and `prop:index-harmonic-pairing` |
| TAUM `thm:partner-unitary-equivalence`, `cor:spectral-projection-pairing`, `cor:eigenvector-pairing`, `thm:closed-range-criterion` | The proposition, its proof and the closed-range criterion immediately following |
| TAUM `def:Fredholm-channel`, `thm:index-formulas`, `cor:odd-Fredholm`, `eq:supertrace`, `thm:heat-index`; LG `def:fredholm-awareness`, `thm:index-as-phenomenal-excess`, `id:index-unpaired-count` | `eq:index-graded`, `eq:index-heat` and their single-channel interpretation |
| LG `eq:bounded-transform-awareness`, `eq:unilateral-shift-example`, `eq:shift-index-minus-one` | §7.3, including the ground-preserving family `eq:index-shift-family`; the extension calculation follows in `eq:index-extension-counts` |
| LG `rem:annihilation-not-destruction`, `eq:shifted-physical-return-model` | §7.1 explains why the forward differential vanishes on the output sector and tests the shifted input return on the active support. The full output return has the final defect as its zero eigenspace. |

### III. Additional structures and their consequences

These are optional extensions of the core, not a cumulative list every realization must satisfy.

**III.1 Reference vectors, grounds, and the ground bridge.** Separate a selected unit vector, a cyclic representation,
a Hamiltonian ground, uniqueness, and a spectral gap. Develop reference fidelity and energy estimates under the
assumptions they actually use. Introduce selected-line bridges, restricted polar orientations, weak versus strong
compatibility, the reference-plane Pauli/phase geometry, and ground inclusion in the active graph. LG's R1 and R2 select a unique cyclic ground in each supplied representation. R3 sends the input ground to a nonzero multiple of the output ground. Keep strong polar compatibility under its additional hypothesis. V adopts these three conditions alongside R4 and R5 from II.3.

The working ground section (`sec:ground-correspondence`) integrates LG’s restricted ground channel, five strong-compatibility conditions and ground-inclusive correspondence with TAUM’s conditional reference-fidelity energy identity and gap bound. It explains the cognitive interpretation of the weak bridge and tests whether it places those grounds on the aligned graph. Section 4.4 gives the neutral rays of the ground plane and the weights of its balanced vector on the full polar graphs. Section 4.5 states the energy–fidelity identity and the gap bound on the full form domain.

Section 5 constructs the cyclic representation of the ground state, relates variance to its centered orbit, and gives examples separating uniqueness, cyclicity and a gap. Section 5.3 asks whether the luminous state is pure or faithful and proves the two properties can coexist only on the scalar algebra. The five-site example distinguishes generation from the ground from the ability to discriminate states by direct measurement. The interpretation of appearances as determinations of luminosity is introduced here. The disposition-dependent generators, conceptual judgement and grasper–grasped frame remain with Part V.

Section 5.1 uses the variance to decide when a bounded operator preserves the ground line. Section 6.5 gives the corresponding criterion for an entire unitary group, including what can be inferred from two return times. Section 6.6 states when a moment’s chosen ground projector is the whole bottom spectral projection.

**III.2 Supplied symmetries and moving subspaces.** Compact symmetry averaging, invariant compression, moving
projection fields, frame covariance, Grassmann curvature, and characteristic classes. The group/action and, for
bundle geometry, the base and projection field are supplied data. Distinguish the mathematical connection from a
physical Yang-Mills identification. A global symmetry may constrain an interface without a spacetime-local theory.

Working §10 (`sec:moving-subspace`) carries TAUM’s moving decompositions and characteristic classes into the manuscript. It derives the curvature from the tangent maps, constructs parallel transport and proves integer first Chern flux for finite-rank families. The pair example (`eq:moving-pair-flux`, `eq:moving-pair-holonomy`) gives the Chern number of the sphere of distinct planes and relates the excitation’s geometric phase to its spectral weight in §8. Working §11 (`sec:symmetry-average`) carries the compact-group average and its compatibility with the sector grading, preparations and dynamics. In the pair, quarter turns suffice for the same invariant compression at every coupling. The character projections distinguish fixed vectors from invariant rays and explain how a symmetry can preserve a unique vacuum up to phase. Convergence of projections and compressed responses across lattice refinements remains open.

**III.3 Selected operations and dynamics.** Adjoint/polar return, trace-nonincreasing CP return maps, reference-line
stabilizers, no-unitary-attraction result, the specific deconditioning semigroup, and momentary reference-state
systems with covariance/contraction. The semigroup remains a model choice. Retain Schur-Feshbach elimination as a
separate subsection with its bounded coupling, self-adjoint diagonal generators, and resolvent assumptions.
Do not present reference-preserving CPTP maps or a particular contraction as consequences of the bare channel.

Working §6 proves the trace-norm bounds and generator of canonical deconditioning. It characterizes channels that carry a ground between moments, gives sufficient conditions for stabilization and distinguishes the probability of a successful return from the state prepared on success.

The working section `sec:preparation-effective-dynamics` uses the Hilbert Schur complement to compute the measured resolvent on the prepared plane. The TAUM source section `sec:effective-dynamics` remains assigned to III.3 for its general signed factorization, including the Krein branch, the oriented return and the opposite dissipative signs. III.4 can then reuse the positivity of the self-energy. Sections 8 and 10 now give both returns in TAUM’s `thm:two-return-forms`. The Schur self-energy in `eq:prepared-schur` inserts the complementary resolvent between the two maps. The comparison following `eq:moving-mixed-return` identifies the identity as the intermediate operation for curvature.

**III.4 Positive response and operational access.** Working §9 (`sec:passive-response`), together with §1, carries the 26 source entries in TAUM’s `sec:boundary-response`. It gives the Herglotz and Schur kernels of a bounded reciprocal port, the spectral and temporal representations, and the unitary uniqueness of the minimal realization. It also treats finite-rank line weights, opposite symmetry sectors and heat-smoothed singular inputs. The two-level example at `eq:passive-signed-witness` gives a negative spectral weight with a real spectrum. The signed-response discussion in III.3 can refer to that example when it distinguishes fixed reciprocal ports from independently calibrated preparation and readout. Section 9 also supplies the Hilbert-space response clause of TAUM’s final architecture theorem (`thm:architecture`). Its general signed reduction remains due in III.3.

TAUM already gives the argument for counting relations among the weight ranges and explains why antiunitary symmetry preserves the matched spectral measures. Working §9 adds the pair example where enlarging the port lowers the indirect dimension and the overlapping-branch example where the chosen analysis direction lies partly outside the reached space. The realization result carried over here is uniqueness on the reached space. A characterization of all realizable functions would require a separate existence theorem.

### IV. Checked mathematical and physical realizations

Begin with the scalar and abstract three-state examples. Then place the finite pair construction, its self-channel
completion, controls, and measurement support conditions beside the continuum orbital parent, selected carrier,
bias window, and filtered readout. They test consistency and clarify what is independently supplied; they do not
derive the cognitive interpretation or microscopic GTD.

Preserve the consequential limitations: Singh's minimal symmetric completion removes the pair's reachable defect;
the orbital parent is a continuum model with a selected finite carrier; the readout needs its filter; and full-matrix
cyclicity does not by itself select a ground. Keep the finite preparation-invariance criterion and the distinction
between an input Hamiltonian and a compression of the output Hamiltonian.

Section 8 compares the angle dependence of the pair and orbital spectra. Rotating the preparation on LG’s exact filtered carrier gives distinct tests with fixed and co-rotated detectors. For the bare orbital shell, parity gives a two-quantum gap to configurations that can mix with it. The section bounds errors in energies, leakage and amplitudes, then illustrates the singlet-gap correction at second order.

Working §12 (`sec:four-maps`) compares a scalar map, the translation generator, oscillator creation and the derivative on a twisted interval. Its three-state preparation distinguishes the orbit span from the largest missed-direction weight in a finite window. The line and oscillator show how spectral multiplicity constrains a scalar port. The oscillator heat trace yields the completed zeta factor on its half-plane of convergence. The smoothed interval response distinguishes opposite twists that have the same mass.

Working §13 (`sec:pair-currents`) shows which one-body mixers preserve the ground and when the zero-imbalance outcome measures the weight on the missed line. The trap comparison distinguishes the two chiral states reached by the excitation at nonzero field from the three states reached when the preparation also supplies the ground. The ledger locates all 57 tracked clauses from LG’s pair-transfer and orbital-source subsections.

VI can discuss the exploratory colour encoding and the remaining task of establishing its readout. Placement of the separate time project awaits a review of its connections to the architecture.

### V. Cognitive realization and Buddhist interpretation

Adopt the five realization hypotheses explicitly for the cognitive model and explain their motivation and costs.
Retain the doctrinal map and its sources; the ground/emptiness/luminosity argument; the appearance and determination
algebras; dispositions and generator classes; the operator-stage model; derived grasper/grasped frames; sems and
rig pa; reflexive terminality; non-alteration, deconditioning, momentariness, and the Third/Fourth Yoga account;
Volition, actualization, and the proposed interpretation of excess/index.

The arguments in working §6 have three planned homes. Preserve the labels and the links from the two-flash discussion during reorganization.

| Planned home | Material now in §6 |
| --- | --- |
| III.1, with the determination algebra | Define the joint presentation space in `eq:joint-presentation-space`. Use `eq:pair-two-setting-ranges` to show how a common defect disappears. Prove `prop:cyclic-common-excess` here and keep the five-site illustration. Its first equivalence needs a cyclic vector. Locating the vector inside the joint span gives the stronger conclusion about an operation leaving that span. |
| III.3, with dynamics | State `eq:two-moment-correspondence` for specified transitions and pairings. Prove `prop:two-moment-deconditioning` for a fixed pairing, including the unique physical relaxation. Keep the proof of all nine semigroup clauses together with stabilization across moments and the restricted and global return after relaxation. The [deconditioning ledger entry](inventory/README.md#worked-entry-for-lgs-canonical-deconditioning-theorem) locates each clause. |
| V, with the cognitive argument | Keep `conj:asymmetric-nonduality`, the cost of enlarging the physical sector, the grounds for counting a missing direction as an appearance, and the two flashes together. Cite the algebra and dynamics results when asking which operation first creates an excess component and whether it affects the judgement. In the move, clarify “include the observer” as including the observer’s bodily state and physical history in Ψ. |

Refer internally to general theorems instead of repeating them in new notation. Keep elementary mathematics local
when it belongs to a specifically chosen cognitive model: extracting every lemma into Part II would obscure the
reason for its assumptions. Operational, causal, and doctrinal identifications remain distinguishable, and none
becomes an empirical result merely because its mathematical objects have been constructed.

#### Grasper–grasped frame and recollection (planned)

Take up the question left by §6’s two flashes. On the reading adopted here, the immediately following mental perception takes the earlier seeing as its object. Both the flash as it appeared and the experiencing of it belong together to the objective aspect of the new cognition. The new act is itself manifest and can in turn be remembered. Locate the apprehending aspect (*grāhakākāra*) within Φ. Self-manifestness belongs to the whole act. For the distinction between the two occurrences and the possibility of further recollection, return to Jinendrabuddhi on *Pramāṇasamuccaya* I.11ab and the dissertation, p. 550 nn. 96–97.

Use the Monday–Tuesday example to pose the open question. Monday’s and Tuesday’s seeings can have the same ray. What in the history of Wednesday’s recollection makes it a recollection of Tuesday’s seeing? Distinguish this question from what the later brightness judgment can discriminate.

- **Relation between rays and matched effects.** For an illustrative isometry $W$ between phenomenal state spaces, let $B=WW^*$ and $q=\|B\eta\|^2>0$ for the unit later state $\eta$. Require $[B\eta]=[W\phi]$ for the unit earlier state $\phi$. By §2, the conditional probability of $WEW^*$ given $B$ is $\langle\phi,E\phi\rangle$ for each earlier effect $E$. This leaves the rest of the later act unspecified.
- **Whose vagueness?** If the current act’s vagueness effect $G'$ satisfies $G'W=WG$ for the earlier act’s vagueness effect $G$, then $p_{\rm current}=q\,p_{\rm earlier}+(1-q)r$, with $0\le r\le1$. The range of $W$ reduces $G'$, so the two contributions add. A vivid recollection of a vague seeing can test that identification. A vivid recollection of a vivid seeing still involves two acts.
- **Supporting note on a sharp aspect.** For a paired unit state $\phi\in D$, an aspect projection $P$ and $w=\|P\phi\|^2$, the squared norm of $(I-\Pi_D)P\phi$ is at most $w(1-w)$. The outside components of $P\phi$ and $(I-P)\phi$ cancel in the whole state. The complement $I-P$ may include more than the other differentiated aspect. A component outside this particular $D$ leaves the common-excess conjecture open.

### VI. What would complete the physical account

One synthesis and open-problem section, merging the existing summaries. Identify which observables, operations,
interfaces, dynamics, symmetries, and couplings are supplied and which have been derived in each realization.
Connect the long-term unification aim to actual missing restrictions rather than to the breadth of examples alone.

The current gauge question belongs here as an investigation: can an independently justified redundancy or constraint
select currently free structure while preserving the intended observable content? Keep distinct frame covariance,
phase-insensitive selected observables, invariant operators, covariant operations, and a genuine BRST/ghost
realization. Neither the negative sector sign nor nilpotency establishes that identification.

The seven formal open questions from LG have the following destinations. The source labels preserve the distinction between a question about a chosen cognitive realization and a condition on a physical protocol.

| Source question | Destination and issue to preserve |
| --- | --- |
| `q:matter-algebra` | III.1 introduces the represented algebra needed for cyclicity. VI asks which physical operations select a nontrivial choice in a realization. |
| `q:luminous-purity` | Section 5.3 leaves both properties open under R2. If the chosen algebra also admits a nonzero noncyclic vector, its cyclic ground state cannot be pure. V considers the choice of algebra. |
| `q:phenomenal-algebra` | III.1 states what cyclicity requires. V asks which global algebra of phenomenal determinations to choose. |
| `q:operational-structural-phase` | II.2 gives the criterion for detecting relative phase. VI asks which effects or controls realize it. |
| `q:content-affect-algebra-choice` | V asks which subalgebras represent object content and affect in a given realization. |
| `q:projection-or-effect-duality-frame` | V asks which frame represents the grasper–grasped distinction inside Φ, including the choice between sharp and unsharp effects. |
| `oq:status-phenomenal-excess` | V asks which phenomenal determinations, if any, densely generate the final defect and what those directions represent. |

### Appendices

- **A. Domains and technical proofs:** retain TAUM's domain appendix and integrate any additional LG requirements.
  Long physical-realization proofs may live here if that improves the main argument; they remain inside the TeX.
- **B. Notation and operator types:** one charter for Ψ, Φ and A. Distinguish all projection types.
- **C. Attribution and verification scope:** preserve scholarly provenance and precise primary citations. State what
  computations check and where they do not replace proof. Consolidate the bibliography and remove obsolete
  references to the superseded manuscript as though it were an external premise.

## Result families needing deliberate coalescence

| Family | Source anchors | Migration decision |
| --- | --- | --- |
| Oriented algebra | TAUM `sec:algebra`; LG `prop:imported-algebra` | One domain-aware statement, not two versions distinguished only by the channel letter. |
| Polar/reflection geometry | TAUM `def:polar-phases`; LG active graph and ground-polar sections | General polar geometry once; reference-line and ground-compatible specializations separately. |
| Fidelity/excitation | TAUM `sec:reference-state`; LG `sec:ground-fidelity` | Separate reference geometry, eigenvector decomposition, ground positivity, and gap estimates. |
| Reference neutrality | LG `thm:ground-graphs-neutral` | Isometric-graph theorem first; two-dimensional maximality and ground interpretation separately. |
| Born correspondence and phase blindness | LG `sec:born-rule-locus-correspondence`, `sec:paired-observable-algebra`, `sec:ssp-factorization-phase` | Keep general active-support results without requiring global injectivity or strict excess. Preserve selected-algebra qualifications. |
| Momentary stabilization | LG `def:momentary-phenomenal-system`, `prop:exact-realization-absorbing`, `thm:momentary-ground-stabilization` | Generalize presentation to supplied pure references and stipulated CPTP covariance/contraction; retain cognitive interpretation in V. |
| Response-generated defect | TAUM `cor:finite-response-defect`, `prop:permanent-port-invisibility`; LG pair/orbital sections | Prove general response statements once and apply them in IV. Ground data are extra; the codomain choice matters. |
| Index and excess | TAUM `sec:cohomology`; LG `sec:fredholm-index-excess` | Coalesce mathematics; retain injective/Fredholm specialization and distinguish infinite defects. |
| Return and actualization | LG `sec:volition-return-orientation`, `sec:mental-causation-actualization` | General polar/weighted return and CP map first; availability and causal interpretation remain additional. |

## Hypothesis review during migration

For each retained or coalesced formal statement, record its minimal mathematical requirements, how a particular
realization supplies them, and its interpretation. Use the existing six reviewed rows as examples, not as evidence
that the remaining rows have been audited. Explicit local standing hypotheses are acceptable; avoid both hidden
assumptions and repetition of the whole axiom system in every theorem.

The working distinctions are: bare closed channel; selected algebra; selected reference; cyclic representation;
Hamiltonian eigenvector/ground/unique ground/gap; weak bridge; strong polar compatibility; injectivity; final defect;
closed range/Fredholm/trace conditions; supplied dynamics/ports; and smooth projection-bundle data. They do not
form one implication chain. A result mentioning a ground may only use its normalized vector; a result with no ground
terminology may still require a supplied reference or dynamics.

Keep mathematical status separate from provenance: a theorem may be standard, imported, or newly established.
Model choices and doctrinal identifications remain explicit. Proof independence from an interpretation is the rule;
a forward reference to an application is not a logical dependency.

## Safe migration sequence

1. Use the agreed Ψ, Φ and A convention in the consolidated draft and retain source identifiers in the crosswalk.
2. Before editing the manuscripts, preserve their reviewed state with a clearly named snapshot/tag.
3. Establish the consolidated preamble and the canonical homes for the shared channel/polar/reference results.
4. Migrate one result family at a time, preserving domains, support conditions, weak/strong distinctions, citations,
   and the associated checks. Do not silently substitute stronger global assumptions for local ones.
5. Split mixed sections as recorded in the source map. Give each model choice, doctrinal identification and formal open question an explicit destination. Record any decision to retire one.
6. Reconcile the four colliding labels (`app:notation`, `eq:conditional-excitation-vector`, `eq:sector-projections`,
   `eq:two-term-complex`). Retain aliases for genuinely coalesced statements where safe; otherwise use an explicit
   old-to-new mapping. Assign stable labels to entries that currently lack them if checks need to reference them.
7. Update `checks/test_labels.py` for the actual destination and any crosswalk, rather than disabling label checks.
8. Build and inspect the consolidated PDF; run the full mathematical checks after the actual migration. Validate
   references and attribution separately from compilation, which cannot establish semantic correctness.
9. Supersede the current pair only after review. Their existing sources and compiled PDFs remain unchanged in this
   outline branch. No new working copy of either mathematical manuscript is introduced here.
