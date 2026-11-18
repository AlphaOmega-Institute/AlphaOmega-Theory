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
| the same: (iv) trace-norm estimate | `eq:relaxation-distance`, from `eq:relaxation-pure-distance` |
| `eq:ground-deviation-lyapunov` and the sentence after it | paragraph after `eq:relaxation-jumps` |
| `def:ground-nonaltering-channel`, `def:ground-attractive-family` | first paragraph of `subsec:relaxation-unitary` |
| `thm:no-unitary-deconditioning` and the paragraph after it | `eq:relaxation-unitary` and the text after it |
| `rem:reduced-irreversibility` | `eq:relaxation-dilation`, `eq:relaxation-dilation-trace` and the text after them |
| `def:momentary-ground-covariance` | `eq:relaxation-ground-carried`, `eq:relaxation-contraction` |
| `prop:exact-realization-absorbing` | `prop:relaxation-moments` (a) |
| `thm:momentary-ground-stabilization` | `prop:relaxation-moments` (b), `eq:relaxation-product` |
| `def:realization-event-stabilization` | last paragraph of `subsec:relaxation-moments` |
| `def:conditional-projective-actualization`, `prop:return-operation-trace` | `eq:relaxation-return` and the text after it |
| `thm:deconditioning-contracts-excess`, second display | `eq:relaxation-paired-weight` |
| `eq:restricted-return-after-deconditioning` | `eq:relaxation-restricted-return` |
| `eq:global-return-after-deconditioning` | `eq:relaxation-return-bound` |

LG's four-term formula `eq:canonical-deconditioning-channel` is displayed as
`eq:relaxation-blocks`. Its row is already explicit and is left as it is.

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
3. The lower bound in `eq:relaxation-distance`, the exact distance for a pure input
   `eq:relaxation-pure-input`, and the conclusion that neither the constant 2 nor the
   rate γ/2 can be improved. A state that commutes with the ground projection approaches
   at rate γ, and one with coherence approaches at rate γ/2.
4. An explicit one-step dilation, `eq:relaxation-dilation`. LG says that a dilation may
   exist and that the model does not identify its environment. The fragment constructs
   one for each time, keeps LG's caveat, and adds that one unitary group cannot reproduce
   all times on a finite-dimensional larger space.
5. `prop:relaxation-moments` (c): a step that contracts every state carries the ground.
   LG assumes both conditions in its stabilization theorem. Its bound uses only the
   contraction, and the fragment says so.
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
- "Its trace is the probability of a branch only once a measurement with the effect Π_D
  has been realized." LG says the weight becomes a success probability only when the
  return is a branch of a separately supplied instrument.

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
- `checks/fragments/test_relaxation.py` has seven tests and replaces
  `checks/pending/test_deconditioning_owed_clauses.py`, whose three tests it contains.
  On this branch the suite gives 125 passed and 2 skipped, and 127 with `--full`. With
  the four-maps block as well the totals would be 133 and 135.
- The check cites existing `UM` labels only. The fragment's labels can be added to its
  `LABELS` once they are in the manuscript.
- This branch does not add `consolidation/inventory/proposed/README.md`, which the
  four-maps branch adds, so that the two branches do not touch the same file. Both
  branches change the same two lines of totals in `checks/README.md`.
