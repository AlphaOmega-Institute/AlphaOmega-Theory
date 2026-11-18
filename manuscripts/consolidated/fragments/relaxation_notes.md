# `relaxation.tex`: LG's deconditioning theorem and the results that depend on it

Source: LG `thm:canonical-deconditioning-semigroup` (line 6033) with the definitions and
results around it in LG's sections on deconditioning, momentary stabilization and return
after deconditioning. Base: `0cd012d`. Planned home: Section 6, as four subsections
between "Deconditioning and the partner of the luminous ground" and "Remembering the
first flash". The preview numbers them 11.1 to 11.4 because the builder continues from
the last section. In place they are 6.4 to 6.7.

    python manuscripts/consolidated/fragments/build_fragment.py relaxation
    python checks/run_all.py -k relaxation
    python consolidation/inventory/proposed/build_relaxation_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger.

## Source results a reader can now find

| LG | Fragment |
|---|---|
| `thm:exact-rigpa-density` | last sentence of `prop:relaxation-pure-distance` and the sentence after its proof |
| `thm:canonical-deconditioning-semigroup`: continuity | `prop:relaxation-semigroup` (a) |
| the same: generator | `prop:relaxation-semigroup` (b), `eq:relaxation-generator`, `eq:relaxation-jumps` |
| the same: (i) ground fixed, (iii) unique fixed state | `prop:relaxation-semigroup` (c) |
| the same: (ii) exact decay | `eq:relaxation-decay` |
| the same: (iv) trace-norm estimate | `eq:relaxation-distance`, from `eq:relaxation-pure-distance`; sharpened in `eq:relaxation-sharp-distance` |
| `eq:ground-deviation-lyapunov` and the sentence after it | paragraph after `eq:relaxation-jumps` |
| `def:ground-nonaltering-channel`, `def:ground-attractive-family` | first paragraph of `subsec:relaxation-unitary`; `prop:relaxation-carried` for what fixing the ground implies |
| `thm:no-unitary-deconditioning` and the paragraph after it | `eq:relaxation-unitary` and the text after it |
| `rem:reduced-irreversibility` | `eq:relaxation-dilation`, `eq:relaxation-dilation-trace` and the text after them |
| `def:momentary-ground-covariance` | `eq:relaxation-ground-carried`, `eq:relaxation-contraction` |
| `prop:exact-realization-absorbing` | `prop:relaxation-moments` (a) |
| `thm:momentary-ground-stabilization` | `prop:relaxation-moments` (b) and (c), `eq:relaxation-product` |
| `def:realization-event-stabilization` | last paragraph of `subsec:relaxation-moments` |
| `def:conditional-projective-actualization`, `prop:return-operation-trace` | `eq:relaxation-return` and the text after it |
| `thm:deconditioning-contracts-excess`, second display | `eq:relaxation-paired-weight` |
| `eq:restricted-return-after-deconditioning` | `eq:relaxation-restricted-return` |
| `eq:global-return-after-deconditioning` | `eq:relaxation-return-bound` |

LG's four-term formula `eq:canonical-deconditioning-channel` is displayed as
`eq:relaxation-blocks`. Its row is already explicit and is left as it is.

## Revision after review of `acc05ae`

The review found one mathematical error and asked for three tightenings. All four are
made, and the rest of the fragment is unchanged.

1. **Error corrected.** I wrote that a step contracting every state carries the ground
   and that the converse fails. The converse holds with the factor one. If a channel
   carries the ground, the effect that measures the weight outside the next ground has
   zero expectation in the ground, so by positivity it annihilates the ground and lies
   below the complementary projection. The weight outside the ground therefore cannot
   increase for any state. This is the reviewer's argument. It is now
   `prop:relaxation-carried`, stated as an equivalence for a channel between two spaces,
   and `prop:relaxation-moments` (c) applies it. What the unitary example shows is that
   carrying the ground need not give a strict decrease. The product bound along a given
   sequence is as before.
2. **Sharpness.** The sentence "both bounds are attained" was wrong for the upper bound.
   The lower bound is attained by states commuting with the ground projection. The upper
   bound is an equality for pure states at time zero and is strict at every positive time
   off the ground. The fragment now proves this with a sharper bound,
   `eq:relaxation-sharp-distance`, which is an equality for every pure input and contains
   the earlier pure-input formula. The constant 2 and the rate γ/2 are optimal because
   the ratio tends to one for pure inputs approaching the ground.
3. **Attraction.** A family is now called ground-attractive on a basin when it fixes the
   ground and the weight outside it tends to zero there, as in LG's definition.
4. **Return.** The probability belongs to the effect and the state afterwards to the
   instrument. The fragment says so and gives a second instrument with the same effect
   that leaves a different state.

## Change in the clause count

The block I proposed was the theorem's nine clauses. Before writing rows I searched the
whole ledger for clauses that name the theorem, its formula, the trace-norm bound it
rests on, or the results that apply it. That search added fifteen clauses outside the
theorem. They are in the block because each is a statement about this relaxation or a
step of its proof.

`consolidation/inventory/proposed/relaxation_rows.csv` has 21 rows, all proposed as
explicit. Ten were recoverable and eleven were owed. The theorem goes from 3 explicit, 3
recoverable and 3 owed to 9 explicit. If the rows are accepted as they stand and nothing
else changes, the ledger goes from 493 rows (254 explicit, 108 recoverable, 131 owed) to
493 (275, 98, 120), and LG from 132, 88, 84 to 153, 78, 73. No row is added. The
four-maps block, on its own branch, is counted separately, and the two sets of rows do
not overlap.

Clauses near these that the fragment does not carry, all still owed:

- the cognitive identifications (`id:rigpa-density`, `id:nonalteration-doctrinal`,
  `id:continuity-without-substance`, `id:mental-causation`), the reading of the model as
  exhaustion of dispositions (`choice:canonical-deconditioning-model`, second clause) and
  `choice:no-recognition-operator`. These belong to Part V.
- `prop:terminality-nonalteration-distinct`, which concerns the typed differential of
  Section 7 and is recoverable there.
- `thm:exhausted-profile-character` and the two Yoga identifications after it.

## What the proofs need

From the consolidated manuscript: the formula `eq:moment-deconditioning` with its
statement that the map is a channel, `eq:moment-excess-decay` and
`prop:two-moment-deconditioning` (Section 6); the outcome `eq:matched-defect-outcome`
(Section 2); the reference lines `eq:restricted-ground-orientation` and the partner
`eq:actual-ground-partner` with strong compatibility (Section 4); the unequal-cell path
(Sections 4 and 6).

Prerequisites the fragment supplies because the manuscript did not have them:

- The bound of the trace distance to a pure state by twice the square root of the
  weight outside it. LG cites Fuchs and van de Graaf for it. The fragment proves it in
  four lines from the pure case, the triangle inequality and the concavity of the square
  root, and cites the paper. Parts (e) of the first proposition, (b) of the third, and
  the return bound all rest on it.
- The equivalence of zero weight outside the ground with being the ground state. LG
  proves it from positivity and uses it for the unique fixed state. Here it follows from
  the same bound.
- The trace norm and its notation, the returned state, and the definitions of fixing the
  ground and attracting a basin.

Standard facts used without proof: a density operator is a convex combination of the
projections onto its eigenvectors; a bounded generator gives a semigroup as its
exponential; a partial trace; an isometry extends to a unitary when the two orthogonal
complements have the same dimension; a finite sum of oscillating exponentials returns
arbitrarily close to its initial value (Dirichlet's approximation theorem).

Three bibliography entries are added in `relaxation_bibliography.tex`, copied from LG
(Fuchs and van de Graaf; Gorini, Kossakowski and Sudarshan; Lindblad). I could not open
the publishers' pages from this environment, so the entries are as LG gives them.

## Mathematics added to the source

These go beyond LG and are offered for review as additions.

1. Continuity in the operator norm on the trace class, with the bound
   `4(1 - exp(-γt))`. LG states strong continuity. The stronger statement is what makes
   the generator bounded and the semigroup its exponential.
2. The fixed trace-class operators at a positive time are the multiples of the ground
   projection. LG states uniqueness among states.
3. The lower bound in `eq:relaxation-distance`, the sharper upper bound
   `eq:relaxation-sharp-distance` with equality for pure inputs, and the conclusion that
   the constant 2 and the rate γ/2 are optimal while LG's bound is strict at positive
   times. A state that commutes with the ground projection approaches at rate γ, and one
   with coherence approaches at rate γ/2.
4. An explicit one-step dilation, `eq:relaxation-dilation`. LG says that a dilation may
   exist and that the model does not identify its environment. The fragment constructs
   one for each time, keeps LG's caveat, and adds that one unitary group cannot reproduce
   all times on a finite-dimensional larger space.
5. `prop:relaxation-carried`: a channel carries the ground exactly when it increases the
   weight outside the ground for no state. The direction from carrying the ground to the
   inequality is the reviewer's. LG assumes both conditions in its stabilization theorem.
   The first already gives the factors one, the bound uses only the factors, and
   stabilization needs products that tend to zero. The fragment says each of these.
   It also gives a conditional step that fixes the ground and increases the weight
   outside it, to mark that the proposition is about trace-preserving maps.
6. The condition `Σ γ_n τ_n = ∞` for deconditioning steps, with the two examples.
7. A bound for the normalized return, `eq:relaxation-return-bound`, valid from every
   initial appearance state. LG states the limit.
8. The numbers on the unequal-cell path at half relaxation.

## Sentences that interpret

My computations have held up in review better than my interpretive sentences. These are
the ones to read as proposals:

- "A unitary can carry one chosen vector to y. It cannot do so while leaving y in place."
  This restates LG's paragraph after the no-go theorem.
- "The construction shows that an environment of this kind is possible. It does not
  identify one in an observer."
- "In the cognitive model the first is an event of exact realization and the second is
  its asymptotic stabilization." LG's names are "realization event" and "asymptotically
  stabilized in realization". The fragment's own terms are "at the ground" and
  "stabilizes".
- "The return is therefore a choice of operation, and what a later step can register
  depends on that choice as well as on the effect." LG says the weight becomes a success
  probability only when the return is a branch of a separately supplied instrument.
- "A conditional step, normalized after a selected outcome, can fix the ground and
  increase the weight outside it." The example is mine. The determinations of Section 5
  act on amplitudes in this conditional way, which is why I thought the limit of the
  proposition worth marking. I have not drawn that connection in the text.

## Names

`δ_y(σ)` for the weight outside the ground (LG's ground deviation; the subscript keeps
it apart from the δ of the pair model). `𝒢_y` for the generator and `J_j` for the jump
operators. `𝒯_t` for a general family of channels, following the manuscript's `𝒯_n`.
`𝖡` for a basin. `𝓗_e` and `ι_t` for the environment and the isometry; I avoided `W`,
which is the silent space of the pair model. `𝓗_{Φ,n}`, `y_n`, `σ_n`, `δ_n`, `q_n`,
`γ_n`, `τ_n` for the moments. `σ_S` for the returned state.

## Left for integration

- The proof of `prop:two-moment-deconditioning` says that the output of any state
  converges in trace norm to `P_y`. It could cite `eq:relaxation-distance`.
- `checks/fragments/test_relaxation.py` has eight tests and replaces
  `checks/pending/test_deconditioning_owed_clauses.py`, whose three tests it contains.
  On this branch the suite gives 126 passed and 2 skipped, and 128 with `--full`. With
  the four-maps block as well the totals would be 134 and 136.
- The check cites existing `UM` labels only. The fragment's labels can be added to its
  `LABELS` once they are in the manuscript.
- This branch does not add `consolidation/inventory/proposed/README.md`, which the
  four-maps branch adds, so that the two branches do not touch the same file. Both
  branches change the same two lines of totals in `checks/README.md`.
