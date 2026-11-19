# `duality_frame.tex`: a sharp frame on the differentiated part of a state

Sources: LG lines 5236 to 5405 ("Derived grasper–grasped frames inside the phenomenal
sector": the frame, its grading, its vanishing at the ground, the duality generator, the
corollary on the subject ket and the open question on the form of the frame). Base:
`5bfb9db`. Planned home: one subsection at the end of Section 5, after "Families of
generators and the words they form". The preview numbers it 13.1. In place it is 5.5.

    python manuscripts/consolidated/fragments/build_fragment.py duality_frame
    python checks/run_all.py -k duality_frame
    python consolidation/inventory/proposed/build_duality_frame_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger. The scope is the one
agreed for this block: the mathematics of the sharp duality frame, with the support bound
made explicit as allowing the frame to miss part of a differentiated state, and the
normalization of a nonzero frame component. What would count as an independently
existing subject is left to the cognitive model.

## Source results a reader can now find

| Source (LG) | Fragment |
|---|---|
| `def:derived-duality-frame` and the grading `eq:derived-duality-grading` | `eq:duality-frame` and the two sentences after it |
| `thm:duality-frame-vanishes-ground` | `prop:duality-frame-ground`, first statement |
| the split `eq:internal-grasper-grasped-split` | `eq:duality-frame-norms` |
| `choice:duality-generator`, the spectral construction | the paragraph "One way to obtain a frame is spectral." |
| `cor:no-independent-subject-ket`, its mathematical part | the last paragraph, first three sentences |
| `q:projection-or-effect-duality-frame` | the last sentence |

## Change in the clause count

`consolidation/inventory/proposed/duality_frame_rows.csv` proposes four owed clauses of
LG as explicit: `def:derived-duality-frame`, `thm:duality-frame-vanishes-ground`,
`choice:duality-generator` and `q:projection-or-effect-duality-frame`. No row is added.

The ledger goes from 497 rows (395 explicit, 21 recoverable, 81 owed) to 497 (399, 21,
77). LG goes from 247, 6, 52 to 251, 6, 48, and TAUM is unchanged.

Clauses near these that the fragment does not carry:

- `cor:no-independent-subject-ket` stays owed. Its clause is that no independently
  existing subject follows from the aspect projection, which is a question for Part V.
  The text supplies its mathematical part: a nonzero component normalizes to a vector
  determined by the state and the frame, it is the state after a Lüders outcome, it
  changes with either, and there is none at the ground. The evidence of the row could say
  that only the question of an independent subject is owed.
- `choice:vasana-dependent-duality-frame` and `id:derived-duality-location`. The text
  attributes the dependence of the frame on profile and state, and the reading of the two
  projections as aspects, to the cognitive model.

## What the proofs need

From the consolidated manuscript: the ground projection `P_y` and the notation of Section
5; the sector symmetry `J` of Section 3, for the remark that the two gradings differ.

Used without proof: for projections, `P ≤ Q` is equivalent to `P = PQ`; the square root of
a positive operator; the spectral theorem for a bounded self-adjoint operator.

## Mathematics added to the sources

Offered for review as additions.

1. The splitting `eq:duality-frame-norms` of the differentiated part into the two frame
   components and a part outside the frame, with the condition for that part to vanish for
   every state and an example in which it does not. This is the reviewer's emphasis.
2. A zero mean of the grading does not characterize the ground: a balanced vector inside
   the frame has it too, with all its weight in the frame.
3. The vanishing at the ground holds for unsharp effects bounded by `I − P_y`.
4. Every sharp frame is the spectral frame of its own grading, so the spectral
   construction is general.
5. The normalized component as the state after a Lüders outcome, with its probability.

## Sentences that interpret

- "LG organizes the differentiated part of an appearance by two projections."
- "LG reads P_+ and P_- as the grasper and grasped aspects of differentiated content, and
  lets the frame depend on the disposition profile and on the state. That reading and that
  dependence belong to the cognitive model."
- "The ground is not a balanced combination of the two aspects." This is LG's sentence,
  followed here by the computation that supports it.
- "LG lets a realization assign such an operator to each profile."
- "What would count as an independently existing subject is a question for the cognitive
  model."

## Names

`P_+`, `P_-`, `P_d` and `Γ` for LG's grasper, grasped and dualized projections and its
grading, so that the text does not carry the cognitive names in its symbols. `η_φ` as in
LG and in Section 4.5. `E_±` for unsharp effects. `Υ` for the generating operator, because
`D` is the final support and LG's `D_vas` would collide with it.

## Left for integration

- The marker line, preceded by the comment `%%% PREVIEW MARKER`, is to be deleted.
- For a trial I inserted the subsection at the end of Section 5 in a copy of the
  manuscript at `5bfb9db`. It builds to 77 pages without warnings.
- `checks/fragments/test_duality_frame.py` has five tests. On this branch the suite gives
  181 passed and 2 skipped, and 183 with `--full`.
- The check cites existing `UM` labels only. The fragment's four labels can be added once
  they are in the manuscript.
