# `polar_graph_clauses.tex`: clauses of Sections 2 to 4 that a reader had to recover

Sources: LG lines 2160 to 4350, that is "Equivalent rank-one formulations", the chapters
on the restricted ground map and the ground plane, the graph maps, the transport of
observables and the paired algebra, and LG `prop:three-part-phenomenal-decomposition` at
line 7079. Base: `e0b1926`. The fragment has four parts, each headed by a marker that
names its place:

| Part | Place | Length |
|---|---|---|
| A | Section 2, first subsection, after the paragraph that ends "requires further data." | one paragraph |
| B1 | Section 3, end of the first subsection | one paragraph |
| B2 | Section 3, "Paired observables and the sector coordinate", after its second paragraph | two paragraphs |
| C | Section 4, a new subsection after "Grounds selected by two path Hamiltonians" | about one and a half pages, one proposition |

    python manuscripts/consolidated/fragments/build_fragment.py polar_graph_clauses
    python checks/run_all.py -k polar_graph
    python consolidation/inventory/proposed/build_polar_graph_clauses_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger. The preview prints the
four parts one after another with their markers, and numbers the subsection 13.1 because
the builder continues from the last section. In place it is 4.4, or 4.5 if it follows the
fidelity subsection that is waiting for integration.

## Source results a reader can now find

| Source (LG) | Fragment |
|---|---|
| `prop:pure-projector-correspondence` | A, `eq:matched-pure-projections` |
| `thm:alphaA-star-isomorphism`, spectra and functional calculus | A, the sentences from "Second" to "of the same kind" |
| `cor:paired-spectral-measures` | A, `eq:matched-spectral-measure` |
| `thm:WA-isometry`, closed range and adjoint | B1, `eq:graph-adjoint-blocks` and the three sentences after it |
| `prop:phase-graph-projection`, the two special projections; `cor:active-graph-decomposition` | B1, `eq:graph-aligned-opposite` |
| `def:paired-observable-embedding`, the unit; `thm:paired-algebra-representation` | B2, first paragraph up to the sentence on the two isomorphisms |
| `thm:paired-algebra-graph-intertwining`, the compression | B2, `eq:graph-paired-compression` |
| `prop:mixed-state-correspondence`, `eq:mixed-sector-compressions` | B2, `eq:graph-density-blocks` and the sentences after it |
| `prop:phase-graph-factorization`, ranges and projectors | B2, "In the sector coordinates ..." |
| `rem:direct-sum-not-partial-trace` | B2, from "The partial trace ..." to the end |
| `prop:ground-bridge-equivalences` | C, first paragraph, "The bridge can also be stated without representatives ..." |
| `lem:ground-channel-bounded`, `thm:ground-channel-polar`, `prop:u0-uniqueness` | C, first paragraph, `eq:ground-plane-restriction` |
| `def:ground-graph-lines`, `thm:ground-graphs-neutral`, `prop:normalized-ground-pair`, `prop:nonduality-without-identity` | C, second paragraph, `eq:ground-plane-vectors` |
| `prop:projectivized-neutral-locus` | C, third paragraph |
| `def:ground-polar-operators`, `thm:ground-polar-algebra` | C, `eq:ground-plane-operators` and the paragraph after it |
| `prop:ground-graph-eigenspaces`, `eq:Xi-Krein-cross-pairing` | C, `eq:ground-plane-eigenlines`, `eq:ground-plane-witt` |
| `id:companion-branch-status`, the clause on equal sector compressions | C, "The sector components of z_- ..." |
| `thm:ground-compression-full-polar` | C, `prop:ground-plane-reduction`, `eq:ground-plane-compression` |
| `prop:three-part-phenomenal-decomposition` | C, `eq:ground-plane-three-parts` |

## Change in the clause count

The block was chosen as the mathematical clauses of LG lines 2400 to 4700 that are not
explicit: 24 recoverable and one owed, after setting aside `cor:BRLC-FS-isometry`, which
the fidelity block already proposes. I then searched the whole ledger for other clauses
that are not explicit and name these results. Four more are covered by the text:

- two recoverable clauses with the status of a proved result whose rows are keyed to
  interpretive headings: the first clause of `id:companion-branch-status` and
  `rem:direct-sum-not-partial-trace`;
- two recoverable clauses outside the line range: `prop:ground-bridge-equivalences` (line
  2189) and `prop:three-part-phenomenal-decomposition` (line 7079).

`consolidation/inventory/proposed/polar_graph_clauses_rows.csv` proposes these 28
recoverable clauses as explicit. No row is added and no owed row changes.

On the base the ledger goes from 496 rows (324 explicit, 78 recoverable, 94 owed) to 496
(352, 50, 94). LG goes from 181, 61, 63 to 209, 33, 63, and TAUM is unchanged. Together
with the fidelity rows the totals would be 497 rows (363, 44, 90).

Clauses near these that the fragment does not carry:

- `eq:W0-SSP-summary` stays owed. Its clause is an identification: the ray of the balanced
  vector as the relational ground and its two compressions as the two aspects. The
  mathematics under it is in part C: the vector from the restricted map alone, its two
  sector components of squared norm one half, and their normalizations. The evidence text
  of the row could say that only the identification is owed.
- The interpretive rows of the same region: `id:relational-nondual-ground`,
  `rem:ground-measurability`, `id:emptiness-luminosity-neutral-relation`, the second
  clause of `id:companion-branch-status`, `id:BRLC-doctrinal`,
  `rem:rejected-sector-measurement-reading`, `id:one-event-two-compressions`, and the
  model choice `choice:determination-algebra`.
- LG `thm:unitary-ground-preservation`. Part C has one sentence on how sector
  Hamiltonians move the points of the circle. It does not state that theorem.
- LG `eq:shifted-physical-return-model`, and TAUM `eq:balanced-branches`,
  `prop:normalized-cycle` and `prop:Omega`, which belong with the fundamental
  decomposition.

## What the proofs need

From the consolidated manuscript: `eq:polar-active-supports` and
`prop:matched-born-probabilities`, including the word "normal" (Section 2); the graph
maps, the Krein form, the polar operators with their zero extensions,
`eq:polar-graph-projector`, `prop:phase-readout`, the sector coordinates, the map
`ι` and `eq:graph-evolution-phase` (Section 3); `eq:restricted-ground-orientation`,
`prop:ground-polar-compatibility`, and from `eq:actual-ground-partner` the partner `r`,
the fidelity `F_x` and the fact that the product of `x` with `r` is real and positive
(Section 4); the unequal-cell path; the pair preparation of Section 1.

Used without proof: uniqueness of the spectral measure of a self-adjoint operator and the
bounded Borel calculus of a normal operator; uniqueness of the polar decomposition of a
bounded operator once the partial isometry has the kernel of the operator; that an
isometry on a complete space has closed range; three values of the zeta function.

Three points where a prerequisite could be missed:

1. LG states the compression `G* ι(X) G = X` as equivalent to the intertwining
   `ι(X) G = G X`. For a single operator the compression is weaker. A block-diagonal
   operator with unrelated blocks `B` and `VLV*` has the compression `(B + L)/2` on every
   graph and maps no graph into itself unless `L = B`. The fragment derives the
   compression from the intertwining and states that formula.
2. LG's uniqueness of `u_0` is a division by `a`. What it excludes is a phase: the
   partial isometries between the two lines are the multiples `e^{iβ} u_0`. The fragment
   says so.
3. The proposition of part C uses `r`, `F_x` and the positivity of their product from
   the subsection "Which input state does the output ground test select?". If the new
   subsection is placed before that one, these become forward references.

## Mathematics added to the sources

Offered for review as additions.

1. Part A: the unbounded self-adjoint case; the spectral measure of the extension by
   zero, whose outcome zero collects the defect and a zero eigenvalue of the input
   observable; the remark that the occupation detector has no such eigenvalue.
2. Part B1: the product of two graph maps at different phases and the squared norm
   `cos²((φ − η)/2)`; the remark that the graph is closed also when the range of `A` is
   not.
3. Part B2: a block-diagonal operator maps a graph into itself exactly when it is paired;
   the state without its off-diagonal blocks is the equal mixture of the graph states at
   two opposite phases; the partial trace written as the sum of the two diagonal blocks
   after the lower one is carried back.
4. Part C, first paragraph: the bounded adjoint of the restriction extends `P_x A*`; the
   example on sequences with `a² = 15/π²` and `F_x = ζ(3)²/(ζ(2)ζ(4))`, where the output
   reference vector is outside the domain of `A*`.
5. Part C: every statement LG makes for the aligned and anti-aligned vectors is made for
   the whole circle of phases where it holds there.
6. Part C: all points of the circle have the same mean energy under sector Hamiltonians
   with the reference vectors as eigenvectors, and the evolution moves the phase at the
   rate of the energy difference.
7. Part C: under the weak bridge alone the compressions of the full polar operators to
   the plane are the plane operators times `√F_x` (`eq:ground-plane-compression`).
   Commutation of the plane projection with the full reflection is by itself equivalent
   to strong compatibility (`prop:ground-plane-reduction`, (b)). LG proves the direction
   from strong compatibility.
8. Part C: the weights of the balanced vector on the two graphs of the full polar map
   (`eq:ground-plane-graph-weights`), with the values 0.974 and 0.026 on the unequal-cell
   path.
9. Part C: the polar map pairs `D ⊖ ℂy` with `S ⊖ ℂr`.

## Sentences that interpret

- "The transport of a pure state is the pure state of the transported vector."
- "A measurement of X on ψ and a measurement of VXV* on Vψ therefore have the same
  distribution of outcomes, and not only the same mean."
- The two sentences on a detector that registers the defect separately and on the
  occupation detector.
- "The two operations are different. The carrier S ⊕ D is a direct sum of two sectors and
  is not the tensor product S ⊗ D of two systems". LG's remark says "presently a direct
  sum, not a bipartite tensor product of two independent systems".
- "The bridge thus fixes the two rays and the relative phase of their representatives".
- "Equal weights in the two sectors define the circle. The bridge selects one point of
  it".
- "Whether R_0 or −iI_0 can be measured is the question of Section 3." LG has a warning
  to the same effect.
- "Mean energy does not distinguish the points either." This is a mathematical statement placed
  where LG has the identification that the anti-aligned companion is not a state below
  the ground. The identification and its other exclusions are not written.
- "Only the component in D⊥ has no partner in S."
- The title speaks of the two grounds. The mathematics of the subsection uses two
  reference vectors and the weak bridge.

## Names

`z^0_φ`, `z_+`, `z_-` for LG's `Ξ_0(θ)`, `Ξ_0^±`, because `Ξ` is the frame of Section 10
and Section 3 writes `z` for vectors of the direct sum. `K_0` for the plane and `P_xy` for
its projection, because `Π_0` is the aligned graph projection in the manuscript. `R_0`,
`I_0` for the polar operators of `u_0`, and `−i I_0` for LG's `Q_0`, as Section 3 writes
`−i I`. `N^0_φ` for the graph lines. "Anti-aligned" is introduced in B1 as a second name
for the opposite graph of Section 3. `E_X` for a spectral measure and `1_Δ` for an
indicator. `E_Ψ`, `E_Φ` for the eigenvalues of the two reference vectors; the fidelity
subsection writes `E_0^Ψ`, `E_0^Φ` for ground energies, and the two could be unified.

## Left for integration

- The four marker lines, each preceded by the comment `%%% PREVIEW MARKER`, are to be
  deleted.
- Section 4 already says that `u_0` is the polar part of `AP_x` and that the balanced
  vector lies on the graph of the restricted line map. Both sentences could point to the
  new subsection, or the first paragraph of part C could move up to follow
  `eq:restricted-ground-orientation`.
- Parts C and the fidelity subsection can stand in either order.
- For a trial I inserted the four parts at the places named above in a copy of the
  manuscript at the base. It builds to 67 pages without warnings.
- `checks/fragments/test_polar_graph_clauses.py` has seven tests. On this branch the
  suite gives 152 passed and 2 skipped, and 154 with `--full`. With the fidelity checks
  it would be 159 and 161.
- The check cites existing `UM` labels only. The fragment's 16 labels can be added once
  they are in the manuscript.
