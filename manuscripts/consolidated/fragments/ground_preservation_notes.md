# `ground_preservation.tex`: returns, the differential, and operators that keep the ground

Sources: LG lines 2240 to 2400 ("The strong polar ground bridge", "The return-kernel
obstruction", "Return operators are not automatically Hamiltonians"), lines 5600 to 5960
("The odd terminality law", "Ground non-alteration and stabilizer structure", "Logical
independence from terminality"), and `def:momentary-phenomenal-system` (line 6157). Base:
`99bf2c3`. The fragment has five parts, each headed by a marker that names its place:

| Part | Place | Length |
|---|---|---|
| A | Section 4.1, after the paragraph that ends "Neither R4 nor R5 is required for the proposition." | one short paragraph |
| B | Section 5.1, after the paragraph that ends "no normalized vector Xy/‖Xy‖." | one paragraph |
| C | Section 6.5, after the paragraph that ends "reduced description of a unitary evolution on a larger space." | one paragraph |
| D | Section 6.6, after the paragraph that ends "The second concerns the particular states σ_n." | three sentences |
| E | Section 7.1, after the paragraph that ends "were specified independently." | three paragraphs |

    python manuscripts/consolidated/fragments/build_fragment.py ground_preservation
    python checks/run_all.py -k ground_preservation
    python consolidation/inventory/proposed/build_ground_preservation_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger. The preview prints the
five parts one after another with their markers. No new subsection, no new bibliography
entry.

## Revision after review of `0c11cb5`

The review kept the block and asked for one correction and four smaller changes. All are
made.

1. **The shifted return acts on the active input.** LG proposes `A*A − a²I` as an
   active-support model, and I had checked it on the whole input space, where a kernel
   lies below `x`. On `S` the answer can differ. For `A = diag(0, 1, 2)` and
   `x = y = e_1`, with `a = 1`, the grounds are strongly compatible. The whole-space
   operator has eigenvalues −1, 0, 3, while on `S` they are 0 and 3 and `x` is the unique
   ground. The text now states the condition on `S`, that `a²` is the bottom of the
   spectrum of `A*A` there, and gives the example. The check computes both. The row for
   `eq:shifted-physical-return-model` is revised. I also added that without strong
   compatibility `x` is not an eigenvector of `A*A` at all, since an eigenvalue would equal
   `a²`, which is condition (i).
2. **The output side needs only R3 and R5.** R5 makes `ker A*` the ground eigenspace of
   `AA*`, at the eigenvalue zero, and R3 puts `y` in the range of `A`, orthogonal to that
   whole eigenspace. Strong compatibility only adds the eigenvalue `a²`. The text and a new
   check with a weak bridge and a two-dimensional defect say so.
3. **Bounded operators in Section 5.1.** LG's hypothesis that `X` is bounded is restored.
4. **Return times in Section 6.5.** The sentence that condition (i) is needed at every
   time is replaced. The return times are those with `|⟨y, U_t y⟩| = 1`, a closed subgroup
   of the real line, so one return is not enough and two returns with irrational ratio
   are. The check tests the 2π example and a three-level case in which 2π is no longer a
   return.
5. **Rank one in Section 6.6.** A projection onto a chosen vector has rank one even in a
   degenerate ground space. What LG requires is that it be the whole spectral projection
   at the bottom, which is the uniqueness part of R2. The text, the check and the row now
   say this.

## Source results a reader can now find

| Source (LG) | Fragment |
|---|---|
| `prop:strong-polar-consequences`, the output return | A, `eq:ground-output-return` |
| `def:ground-line-preserving-operator` | B, `eq:ground-line-kept` and the sentence after it |
| `prop:ground-preservation-reduction` | B, from "The condition is one-sided." to the end |
| `thm:unitary-ground-preservation` | C, `eq:relaxation-unitary-group` and its proof |
| `def:momentary-phenomenal-system`, each moment's ground projection as the whole spectral projection | D |
| `prop:no-iterated-awareness-tower`, the mathematical clause | E, `eq:index-no-second-step` and the four sentences after it |
| `id:rangrig-terminal-law`, the consequence that no second act is needed | E, "The opening section takes ..." and the next sentence |
| `rem:annihilation-not-destruction` | E, second paragraph, first three sentences |
| `prop:terminality-nonalteration-distinct`, the mathematical clause | E, second paragraph |
| `prop:vacuum-not-return-zero-mode` | E, `eq:return-kernels` |
| `eq:shifted-physical-return-model` | E, third paragraph after `eq:return-kernels` |

## Change in the clause count

The block was proposed as ten recoverable clauses. I searched the whole ledger for other
clauses that are not explicit and name returns, nilpotency, terminality, ground
preservation, unitary groups or moments. One more is covered: the owed remark
`rem:annihilation-not-destruction`, whose evidence asks to attribute the zero second leg
to the typing of the map and to separate it from dynamical destruction of the ground.
Part E does that in two sentences. Its status is given as a proved result, because the
statement concerns the typing of the map and its independence from any evolution. LG's
word "ontological" is not carried.

`consolidation/inventory/proposed/ground_preservation_rows.csv` proposes the 11 clauses
as explicit. No row is added.

The ledger goes from 497 rows (377 explicit, 33 recoverable, 87 owed) to 497 (388, 23,
86). LG goes from 229, 18, 58 to 240, 8, 57, and TAUM is unchanged.

Clauses near these that the fragment does not carry:

- the interpretive halves of `prop:no-iterated-awareness-tower`, `id:rangrig-terminal-law`
  and `prop:terminality-nonalteration-distinct`, and the identifications
  `id:nilpotency-rangrig` and `id:nonalteration-doctrinal`, all for Part V. The sentence
  "Whether the vanishing of d_A on the output summand expresses this reflexivity is a
  further identification" leaves the first two open in the text;
- TAUM `prop:differential-basic` (the adjoint of the differential, density, closedness),
  `thm:architecture` (the nilpotent pair and the two odd sums) and `prop:Omega` (the signed
  reference image and its return). They concern the signed orientations and belong with
  the fundamental decomposition;
- LG `prop:exhausted-words-ground-preserve` and `choice:exhausted-profile-ground-triviality`,
  which extend the reducing condition of part B to the words of a profile. They belong to
  the profile apparatus.

After this block eight LG rows remain recoverable. Six are the odd-extension group (LG 523
to 830), which overlaps the fundamental decomposition. The other two,
`prop:word-order-sensitivity` and `def:content-equivalent-states`, belong to the profile
apparatus.

## What the proofs need

From the consolidated manuscript: Proposition 4.1 with conditions (i) and (v); the
variance of Section 5.1 (`eq:determination-variance`); `eq:relaxation-unitary`,
`eq:graph-evolution-phase` and the character in Section 11; the moments of Section 6.6;
the differential and the returns of Section 7.1; the shift example and the channels that
attain `eq:index-channel-coherence`; the unequal-cell path; the pair preparation; R2 and
R5.

Used without proof: the spectral theorem, in the form that the derivative of `e^{-itG}h`
at zero exists exactly for `h` in the domain of `G`; that a continuous homomorphism from
the real line into the unit circle is `t ↦ e^{-itg}`.

## Mathematics added to the sources

Offered for review as additions.

1. Part B: the ground line is kept exactly when the variance vanishes, and the ray exactly
   when in addition the mean is nonzero; the operator `|y⟩⟨w|` as a one-sided example; a
   unitary keeps the line when its adjoint does.
2. Part C: the example `G = diag(0, 1)` showing that a single return is not enough.
3. Part E: the returns on the unequal-cell path and in the pair preparation select other
   states than the grounds; under strong compatibility the shifted output return is `−a²`
   on the defect; the channels attaining `eq:index-channel-coherence` as the example of one
   preparation with ground-carrying and ground-moving channels.
4. Part E, from the review: the example `A = diag(0, 1, 2)` separating the shifted return on
   the active input from the one on the whole input, and the exclusion of `AA*` from R3 and
   R5 alone.
5. Part C, from the review: the return times as a closed subgroup of the real line.

## Sentences that interpret

- "In this sense the model generates no sequence x ↦ Ax ↦ A(Ax) ↦ ⋯. An output vector is
  not an input to the same map."
- "The opening section takes a cognition to be manifest in its own occurrence. On that
  reading no second act is needed to make it manifest." This is the recoverable clause of
  `id:rangrig-terminal-law`. The identification itself stays owed, and the next sentence
  says so.
- "It does not describe an evolution that removes y. It records that the map is not applied
  to output vectors."
- "In the examples the returns select other states."

## Names

`eq:ground-output-return`, `eq:ground-line-kept`, `eq:relaxation-unitary-group`,
`eq:index-no-second-step`, `eq:return-kernels`. "Keeps the ground line" and "keeps the
ground ray" for LG's "preserves the ground line" and "preserves the ground projectively".
`g` for the eigenvalue of the generator, as LG's `g_0`. Part B uses `w` for a unit vector
orthogonal to the ground; Section 4 does not use `w` nearby.

## Left for integration

- The five marker lines, each preceded by the comment `%%% PREVIEW MARKER`, are to be
  deleted.
- Part A refers forward to the unequal-cell path of Section 4.3, and part E to the shift
  example of Sections 7.3 to 7.5.
- For a trial I inserted the five parts at the places named above in a copy of the
  manuscript at `99bf2c3`. It builds to 74 pages without warnings.
- `checks/fragments/test_ground_preservation.py` has seven tests. On this branch the suite
  gives 172 passed and 2 skipped, and 174 with `--full`.
- The check cites existing `UM` labels only. The fragment's five labels can be added once
  they are in the manuscript.
