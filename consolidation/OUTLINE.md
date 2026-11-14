# One consolidated Universal Model: first migration outline

This is an editorial plan for one self-contained manuscript superseding the current TAUM/LG pair. It is not a new
companion paper and supplies no premises needed by either current manuscript. Neither manuscript is edited by this
proposal.

## Source and review scope

- Both manuscript sources: Theory commit `8e9fde9`.
- Initial statement inventory and six proof annotations: `58b8354`.
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
3. Notation remains open. Nothing here chooses `B` over `A`, or sector subscripts `+/-` over `Psi/Phi`. Existing labels
   are retained as source identifiers, not as a commitment to their final spelling. A uniform `Psi/Phi/A` convention is
   compatible with stating the general theory independently of its phenomenal interpretation.
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

Explain the original problem of relating physical states and their presentation to an observer, its origin in
Dharmakirti, and the broader ambition toward physical unification. Explain the logical roles of the subsequent
parts before introducing the operator machinery.

Coalesce TAUM's introduction and LG's charter. Keep the distinction between a general architecture and a selected
physical theory; replace the two-document hierarchy with an internal hierarchy of data, hypotheses, consequences,
and interpretation. The opening should identify where additional physical principles would have to make the
architecture selective. It should not imply that all subsequent results depend on all five LG axioms.

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

**II.3 Kernels, cohomology, and index.** Closure-aware Hodge decomposition, partner spectral pairing, reduced
cohomology, bounded transforms, Fredholm index, and conditional heat-supertrace statements. State the roles of
closed range, Fredholmness, and trace-class assumptions where used. The interpretation of a final defect as
phenomenal excess belongs in V, not in the definition of the kernel.

### III. Additional structures and their consequences

These are optional extensions of the core, not a cumulative list every realization must satisfy.

**III.1 References, grounds, and channel subclasses.** Separate a selected unit vector, a cyclic representation,
a Hamiltonian ground, uniqueness, and a spectral gap. Develop reference fidelity and energy estimates under the
assumptions they actually use. Introduce selected-line bridges, restricted polar orientations, weak versus strong
compatibility, the reference-plane Pauli/phase geometry, and ground inclusion in the active graph. State injectivity
and nonzero final defect as independent subclass conditions. LG's R1-R5 receive neutral formulations here; their
joint adoption and cognitive meanings are made explicit in V.

**III.2 Supplied symmetries and moving subspaces.** Compact symmetry averaging, invariant compression, moving
projection fields, frame covariance, Grassmann curvature, and characteristic classes. The group/action and, for
bundle geometry, the base and projection field are supplied data. Distinguish the mathematical connection from a
physical Yang-Mills identification. A global symmetry may constrain an interface without a spacetime-local theory.

**III.3 Selected operations and dynamics.** Adjoint/polar return, trace-nonincreasing CP return maps, reference-line
stabilizers, no-unitary-attraction result, the specific deconditioning semigroup, and momentary reference-state
systems with covariance/contraction. The semigroup remains a model choice. Retain Schur-Feshbach elimination as a
separate subsection with its bounded coupling, self-adjoint diagonal generators, and resolvent assumptions.
Do not present reference-preserving CPTP maps or a particular contraction as consequences of the bare channel.

**III.4 Positive response and operational access.** Hilbert-self-adjoint realizations with bounded ports; Herglotz
and Schur kernels; spectral measures; minimal realizations; finite response-defect count; permanent invisibility;
symmetry-resolved departure; and heat-smoothed singular ports. Keep the ambient invisible complement distinct
from the response-generated final defect. Krein self-adjointness alone does not establish positive response.

### IV. Checked mathematical and physical realizations

Begin with the scalar and abstract three-state examples. Then place the finite pair construction, its self-channel
completion, controls, and measurement support conditions beside the continuum orbital parent, selected carrier,
bias window, and filtered readout. They test consistency and clarify what is independently supplied; they do not
derive the cognitive interpretation or microscopic GTD.

Preserve the consequential limitations: Singh's minimal symmetric completion removes the pair's reachable defect;
the orbital parent is a continuum model with a selected finite carrier; the readout needs its filter; and full-matrix
cyclicity does not by itself select a ground. Keep the finite preparation-invariance criterion and the distinction
between an input Hamiltonian and a compression of the output Hamiltonian.

Retain the translation, oscillator, and quasiperiodic examples for their distinct spectral/domain lessons. Arithmetic
material remains an explicitly limited example/open direction unless a later editorial decision removes it.
An expanded introduction does not license silently dropping it.

The separate colour and time projects can be situated here or in VI with exact statements of their relationship to
the architecture. Their full texts are not automatically absorbed by this TAUM/LG merge, and no unverified bridge
to either project is asserted.

### V. Cognitive realization and Buddhist interpretation

Adopt the five realization hypotheses explicitly for the cognitive model and explain their motivation and costs.
Retain the doctrinal map and its sources; the ground/emptiness/luminosity argument; the appearance and determination
algebras; dispositions and generator classes; the operator-stage model; derived grasper/grasped frames; sems and
rig pa; reflexive terminality; non-alteration, deconditioning, momentariness, and the Third/Fourth Yoga account;
Volition, actualization, and the proposed interpretation of excess/index.

Refer internally to general theorems instead of repeating them in new notation. Keep elementary mathematics local
when it belongs to a specifically chosen cognitive model: extracting every lemma into Part II would obscure the
reason for its assumptions. Operational, causal, and doctrinal identifications remain distinguishable, and none
becomes an empirical result merely because its mathematical objects have been constructed.

### VI. What would complete the physical account

One synthesis and open-problem section, merging the existing summaries. Identify which observables, operations,
interfaces, dynamics, symmetries, and couplings are supplied and which have been derived in each realization.
Connect the long-term unification aim to actual missing restrictions rather than to the breadth of examples alone.

The current gauge question belongs here as an investigation: can an independently justified redundancy or constraint
select currently free structure while preserving the intended observable content? Keep distinct frame covariance,
phase-insensitive selected observables, invariant operators, covariant operations, and a genuine BRST/ghost
realization. Neither the negative sector sign nor nilpotency establishes that identification.

### Appendices

- **A. Domains and technical proofs:** retain TAUM's domain appendix and integrate any additional LG requirements.
  Long physical-realization proofs may live here if that improves the main argument; they remain inside the TeX.
- **B. Notation and operator types:** one charter, after the notation decision. Distinguish all projection types.
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

1. Review this outline and the unresolved notation; no notation migration is made by this branch.
2. Before editing the manuscripts, preserve their reviewed state with a clearly named snapshot/tag.
3. Establish the consolidated preamble and the canonical homes for the shared channel/polar/reference results.
4. Migrate one result family at a time, preserving domains, support conditions, weak/strong distinctions, citations,
   and the associated checks. Do not silently substitute stronger global assumptions for local ones.
5. Split mixed sections as recorded in the source map. Carry all model choices and identifications forward to an
   explicit destination; reconsideration must be a deliberate editorial decision.
6. Reconcile the four colliding labels (`app:notation`, `eq:conditional-excitation-vector`, `eq:sector-projections`,
   `eq:two-term-complex`). Retain aliases for genuinely coalesced statements where safe; otherwise use an explicit
   old-to-new mapping. Assign stable labels to entries that currently lack them if checks need to reference them.
7. Update `checks/test_labels.py` for the actual destination and any crosswalk, rather than disabling label checks.
8. Build and inspect the consolidated PDF; run the full mathematical checks after the actual migration. Validate
   references and attribution separately from compilation, which cannot establish semantic correctness.
9. Supersede the current pair only after review. Their existing sources and compiled PDFs remain unchanged in this
   outline branch. No new working copy of either mathematical manuscript is introduced here.
