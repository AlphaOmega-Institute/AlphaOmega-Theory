# Statement inventory for the TAUM–LG consolidation

Working input for a merged outline. It lists every formal statement in the two manuscripts so that each result can be
given one home, with what it needs and what it contributes. The manuscripts are not modified by anything here.

The first proposed organization and section-level migration map are in [the consolidation outline](../OUTLINE.md).
That map is not a completed proof-by-proof dependency audit.

```sh
python consolidation/inventory/inventory.py     # from the repository root; rewrites statements.csv
```

Standard library only.

| File | What it is |
| --- | --- |
| `inventory.py` | Reads both manuscript sources and writes `statements.csv` |
| `statements.csv` | Generated: one row per theorem-like environment (281 rows, including model choices and doctrinal identifications) |
| `annotations.csv` | Written by hand: the reviewed columns, keyed by manuscript and label |

## Columns

**Extracted from the source.** `manuscript`, `line`, `section`, `kind`, `title`, `label`, `status`, `has_proof`.
`label` is the first label inside the environment; a few definitions carry only an equation label. `status` is LG's
explicit status label where one follows the statement.

**Automated prompts.** These say where to look. They are not dependency results.

- `auto_refers_to`: statements whose labels this statement or its proof cites.
- `auto_cites_axiom` (LG): whether one of LG's five axioms is cited directly or through other statements.
- `auto_ground_objects`: whether ground vectors, ground projections or ground-plane objects appear in the statement, in
  its proof, or in a statement it refers to. Matched by macro, by raw notation and by the words ground, vacuum, luminous.
- `auto_taum_cited_nearby` (LG): TAUM is cited within a few lines.

**Reviewed.** Filled only from `annotations.csv`, by someone who has read the statement and its proof.

- `proof_minimally_needs`: the hypotheses the proof actually uses.
- `realization_supplies`: how the realization provides those objects.
- `interpretation`: the reading the realization gives the result.
- `basis`: what was read, and by whom the point was made if not the annotator.

## What the automated columns can and cannot tell you

Two examples show both directions of error.

- *A proof can need less than the surrounding axioms.* `thm:ground-graphs-neutral` sits among statements that use both
  grounds, and its `auto_ground_objects` is "statement". Its proof uses only that the selected map between the two lines
  is unitary; uniqueness, cyclicity and energy of the ground vectors play no part. The flag records where the objects
  come from in LG, not what the theorem requires.
- *A pattern can miss objects written another way.* An earlier version of the ground-object pattern matched macros only
  and reported the momentary-continua section as free of ground objects. That section writes its luminous ground
  projections in raw notation. The pattern now includes raw notation and the three words above, and can still miss
  others.

Citations are a weaker guide still. Of LG's 170 statements other than its axioms, 16 cite an axiom directly or through
other statements, while 90 show ground objects. A missing citation is not a missing hypothesis: standing assumptions
and earlier definitions legitimately carry it. The point for the consolidated text is to make each result's hypotheses
identifiable without repeating the axiom system in every theorem.

Neither flag measures dependence on injectivity or on strict excess.

## Current counts

| | Statements | With a proof | With an explicit status label |
| --- | --- | --- | --- |
| TAUM | 106 | 52 | 0 (TAUM uses a summary table) |
| LG | 175 | 71 | 14 |

LG's total includes 29 doctrinal identifications and 17 model choices. Both kinds are part of the consolidation;
excluding them would omit formal statements of the interpretive commitments and selected models. The extractor
refuses any used, declared theorem-like environment that is not in its list of supported kinds.

Six statements are reviewed so far. Add a row to `annotations.csv` after reading a proof; the script refuses annotations
whose label does not exist.

Thirteen TAUM remarks have no label. Their temporary annotation keys are `@` followed by the source line
number (for example, `@100`); the CSV retains a blank label and a separate line column. Such keys are not stable across
source edits. Assign or reconcile stable labels when those entries are migrated, rather than treating line numbers as
permanent identifiers.

## Migration into the consolidated manuscript

[The clause ledger](migration.csv) compares the TAUM and LG sources at `8e9fde9` with the introduction, §§1–17 and Appendix A of the consolidated TeX. It records where a reader can find a claim and what the manuscript still owes. The six annotations above concern a different question, the minimal hypotheses of selected proofs.

Section 9 brought all 26 entries from TAUM’s passive-response section into the draft. Together with §1, it states 40 previously tracked clauses and five from the surrounding discussion. Splitting the Hilbert response in TAUM’s architecture theorem from its signed reduction brought the ledger to 219 explicit, 109 recoverable and 157 owed.

Section 10 carries TAUM’s moving decompositions, projector curvature and characteristic classes. Of 21 previously tracked clauses in those two sections, 19 moved from owed to explicit and two were already explicit. Six new rows cover source passages that lacked anchors. Three clauses elsewhere in TAUM now have explicit counterparts. The polar-geometry remark distinguishes its algebraic sign from holonomy, the two-return-forms theorem compares tangent returns with the Schur self-energy, and the architecture theorem collects the curvature and Chern-character results.

Section 11 brings five owed TAUM clauses into the manuscript, including the compact-symmetry clause of its architecture theorem. One new row tracks the distinction between the Haar projection and the other projections in the construction. The discussion of character projections also states LG’s result that a symmetry commuting with a Hamiltonian acts on its simple ground by a character.

Section 12 gives TAUM’s scalar, translation, oscillator and interval maps, the three-state preparation and the closing comparison table. It makes 16 owed and three recoverable clauses explicit and adds two rows for previously untracked passages. A separate LG row reserves for Part V the qualification that the phase of a unique ground does not explain individual differences.

Section 6 now gives all nine clauses of LG’s canonical deconditioning theorem. The theorem and the surrounding arguments on pure-ground distance, stabilization across moments and conditional return bring 11 owed and ten recoverable LG clauses into the text.

Section 13 makes 11 owed and 17 recoverable clauses explicit. All 30 tracked clauses in LG’s finite pair-transfer subsection and all 27 in its orbital-source subsection now have explicit counterparts.

Sections 2–5 now state 45 previously recoverable and seven previously owed clauses on polar transport, graph geometry, ground fidelity and cyclic representations. A new TAUM row covers which Hamiltonians the fidelity result applies to. LG’s luminous-purity question is explicit as an open question in §5.3.

Sections 4–7 make ten recoverable LG clauses and one owed clause explicit. They give the criteria for an operator or a unitary group to preserve the ground, state the spectral requirement for a unique ground at each moment, and test the return operators as candidate Hamiltonians. Section 7.1 explains why the forward differential vanishes on the output summand before any dynamics are specified. Section 15 gives LG’s interpretation of this property as the no-regress aspect of reflexive awareness.

Section 5.4 makes two recoverable and five owed LG clauses explicit. It defines the space generated by a family, compares the two orders of a pair of generators, and states when two density operators agree on a chosen content algebra. Under LG’s assumption on the generators of its exhausted profile, the ground state restricts to a character. The choice of content and affective algebras is recorded as an open question.

Section 5.5 makes four owed LG clauses explicit. It defines the sharp duality frame, proves that both aspects vanish on the ground, gives a spectral construction and states the open question about the choice of frame. The section identifies which differentiated directions a fixed frame can miss. Section 16 explains what the conditional aspect ray contributes to the question of an independently existing subject.

Appendix A and the accompanying passages in §§2, 3, 7, 10 and 12 make 16 recoverable and 21 owed clauses explicit. One new TAUM row tracks the fundamental decomposition and its sector projections. The two odd branches now have proofs on their full operator domains, and the polar factors have bounded extensions on the active carrier. The oscillator compares a single polar rotation with the mode-dependent rates of the weighted generator. This left the general signed Schur reduction for §8.

The signed Schur reduction in §8 makes five recoverable and ten owed TAUM clauses explicit. It gives both adjoint signs, the full inverse on the common resolvent set and the positivity properties of the return. All 192 tracked TAUM clauses are explicit. This left 46 owed LG clauses with a planned home in Part V. The heading queue below is a separate record of passages outside the formal entries.

Section 14 makes sixteen owed LG clauses explicit. They cover the claim statuses, the names of the two sectors, the vacuum and luminosity, the relational ground and its companion, the exclusion of a lower energy level for the Fourth Yoga, the reading of the active graph as a cognitive event, and the open status of a single map’s defect. One new row records the conjunction at `eq:unconditioned-ground-two-conditions` under the exhausted profile. LG supplies the two endpoint conditions and leaves the evolution of the family unspecified.

Section 15 makes eighteen owed LG clauses explicit. It locates exact realization at the ground, distinguishes a later judgment or record, and states what canonical deconditioning can attain in finite time. It relates LG’s readings of continuity across moments and of the Third and Fourth Yogas to the conditions on the state, the channels and the active algebra. The Fourth Yoga row records both an identification and an open question. The draft retains the shared ground and the exhausted profile at the Fourth, and leaves open how to represent the cessation of the undifferentiated appearance of luminosity. The ordinary attainability of exact ground orthogonality is recorded as an open question.

Section 16 makes nine owed LG clauses explicit. It models how a disposition profile permits a sensory image and its conceptual identification, constructs a judgement through exclusion and superimposition, and tests its conditional state against specified content effects. The grasper and grasped are aspects of the same event.

Section 17 makes the final three tracked LG clauses explicit as interpretations. It distinguishes the Hilbert adjoint from the polar return, interprets return through the final support as constrained actualization, and asks which physical preparations count as presentations. The manuscript’s conjecture about common phenomenal excess keeps its status as a conjecture.

| Source | Explicit | Recoverable | Owed |
| --- | ---: | ---: | ---: |
| TAUM | 192 | 0 | 0 |
| LG | 306 | 0 | 0 |
| UM | 1 | 0 | 0 |
| **All 499 clauses** | **499** | **0** | **0** |

The UM row records the conjecture introduced in the consolidated manuscript.

Filter by `source` and `source_key` to read all clauses of an entry. For TAUM and LG, the key is the source label or the inventory’s `@line` key for an unlabelled environment. The source line locates the environment or supplementary label. `UM` identifies a claim introduced in the consolidated draft. For example, LG’s R5 assumes a missing direction for one specified map. The draft conjectures that one genuine appearance direction is missed by **every** admissible map. The strengthening has its own UM entry.

Three fields answer different review questions.

- `source_kind` records how the source presents a statement, or marks supplementary source material or a draft extension.
- `draft_status` records the role the reviewer finds in the draft passage. A clause can serve more than one role, separated by `|`. For example, R1–R5 are axioms in LG and hypotheses for a realization in the draft. An owed clause has status `absent`.
- `draft_labels` locates the draft passage. An explicit clause with only part, section or subsection labels also carries a short `draft_quote`. Its TeX must occur exactly once after whitespace is normalized, so a reviewer can relocate the passage after a move. `evidence_or_owed` explains what the passage establishes or what work remains. `planned_homes` gives the intended position in the final outline.

| Coverage | How to read the evidence |
| --- | --- |
| `explicit` | The cited passage states or proves the indicated clause. A narrower case gets its own row when the general source result is still due. |
| `recoverable` | The reviewer gives a specialization or short inference from named passages. State the consequence or supply an explicit cross-reference before closing the source entry. |
| `owed` | The draft still needs the indicated statement, proof, construction or interpretation. The destination label is blank until that material is present. |

The ledger includes every formal entry in the source inventory. Supplements track mathematics outside those environments, including the pair Hamiltonian, preparation and readout, and the orbital parent, carrier and response. One source passage can contribute a theorem, a modelling choice and an interpretation with different destinations.

The [heading queue](heading_queue.csv) lists nonempty source heading bodies with no ledger row anchored in them. Each body runs up to the next heading. The screen skips empty or label-only bodies and the bibliography, with no minimum length. Read queued passages to decide whether to add a supplement or point to coverage elsewhere. Recheck partially tracked bodies when migrating the family. A single entry removes a heading from the queue even if nearby derivations or prose still need a home.

All 499 tracked clauses now have manuscript locators. The separate queue still has 107 source heading bodies without a ledger anchor in that body. For example, the LG headings “Raw ground return” and “Raw Volition and spectral reweighting” remain in the queue. Section 17 treats their domain and weighting distinctions under rows anchored at LG’s definitions and identifications. The new LG row is anchored in “Recognition without fabrication” at line 6463. TAUM’s domain summary for the Schur–Feshbach generators (lines 4421–4443) repeats the product-domain and bounded-coupling requirements proved in §8. Its ledger rows are anchored at the original propositions, so that summary remains in the generated queue. A new row covers TAUM’s “Fundamental decomposition”. Its unlabelled “Equal-carrier realization” (lines 505–523) is also covered in Appendix A by the Pauli formulas at `eq:balanced-channel-branches`. The row for those formulas is anchored later in TAUM, so this heading still appears in the queue. The fidelity row for “Relation to the transition channel” now covers that TAUM heading. TAUM’s unlabelled closing table (lines 4025–4048) has been reviewed. The `sec:examples` ledger row locates its counterpart at `subsec:four-maps-table` in §12. That row is anchored at the opening of the source section, so the generated queue still lists the closing heading.

### Worked entry for LG's canonical deconditioning theorem

| Source clause | Proof or statement in §6 |
| --- | --- |
| Channel on trace-class operators | `eq:moment-deconditioning` and its expansion `eq:relaxation-blocks`. |
| Complete positivity and trace preservation | Reset and trace calculation after `eq:moment-deconditioning`. |
| Semigroup law | Part (a) of `prop:relaxation-semigroup` composes the block coefficients. |
| Strong continuity on the trace class | Part (a) bounds the change in the induced operator norm on the trace class. |
| (i) Target projection fixed | Part (c) identifies the fixed operators at positive time. |
| (ii) Exact decay of ground-deviation weight | `eq:relaxation-decay` gives the exponential factor for every state. |
| (iii) Unique fixed normal state at any positive time | Part (c) shows every fixed trace-class operator is a multiple of the target projection. |
| (iv) Trace-norm estimate with constant 2 | `eq:relaxation-distance` gives the estimate and `eq:relaxation-sharp-distance` proves sharpness of the constant and of the decay rate. |
| Generator | Part (b) gives the bounded trace-class generator in `eq:relaxation-generator`. |

From the repository root, validate the ledger and check that the reading queue matches it. After adding source entries, regenerate the queue by running `heading_queue.py` without `--check`.

```sh
python consolidation/inventory/check_migration.py
python consolidation/inventory/heading_queue.py --check
```

The checker reports editorial counts. The ledger weights neither the length nor the difficulty of the clauses and makes no percentage-complete estimate.
