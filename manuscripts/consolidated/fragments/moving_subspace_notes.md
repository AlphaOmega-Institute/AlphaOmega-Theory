# `moving_subspace.tex`: TAUM's moving decompositions, projector curvature and characteristic classes

Source: TAUM lines 2008 to 2508, 16 formal entries and 3 supplementary equations in the
ledger. Base: `8a68b22`. Planned home: outline III.2.

    python manuscripts/consolidated/fragments/build_fragment.py moving_subspace
    python checks/run_all.py -k moving_subspace
    python consolidation/inventory/proposed/build_moving_subspace_rows.py

Nothing here changes `Universal_Model.tex`, the shared ledger or the existing checks.

## Source results a reader can now find

| TAUM | Fragment |
|---|---|
| setting of the section, `eq:projection-field`, `eq:moving-J` | `eq:moving-projection-field`, `eq:moving-symmetry`, `prop:moving-unitary` |
| `lem:projection-tangent`, `def:tangent-channel` | `prop:moving-tangent`, `eq:moving-tangent-map` |
| `thm:tangent-realization`, remark at 2155 | `eq:moving-block-forms`, `eq:moving-two-odd` and the paragraph after them |
| `eq:fiberwise-Laplacian` | `eq:moving-squares` and the same paragraph |
| `def:Grassmann-connection`, `thm:projector-curvature` | `eq:moving-connection`, `prop:moving-curvature` |
| `def:mixed-return`, `cor:curvature-AV`, remark at 2238 | paragraph after `prop:moving-curvature`, `eq:moving-mixed-return` |
| `def:local-gauge-data`, `thm:frame-covariance` | `eq:moving-connection-matrix`, `eq:moving-frame-curvature`, `prop:moving-frame` |
| `def:holonomy`, remark at 2326 | `eq:moving-transport`, `eq:moving-frame-transport` and the two paragraphs around them |
| `eq:curvature-form`, `thm:chern-character`, `eq:higher-chern-character`, `eq:ch2-return` | `eq:moving-curvature-form`, `prop:moving-chern` |
| `cor:line-flux`, `eq:physics-flux-integral`, remark at 2477 | `cor:moving-flux` and the three paragraphs after it |

## Change in the clause count

The ledger has 21 clauses for this section: 2 explicit and 19 owed. In
`consolidation/inventory/proposed/moving_subspace_rows.csv` the 19 owed clauses are
explicit and the 2 explicit ones are unchanged. Six further rows cover source passages
that had no row, two of which are entries of the reading queue (TAUM 2008 to 2015 and
2366 to 2389). If the rows are accepted as they stand, the ledger goes from 485 rows
(219 explicit, 109 recoverable, 157 owed) to 491 (244, 109, 138), and TAUM from 180
(87, 20, 73) to 186 (112, 20, 54). These are proposals.

One row needs a ruling. `thm:tangent-realization` asks for "both oriented sums". The
fragment gives them as `d + d*` and `d* - d` in the notation of Section 7. TAUM's own
letters for the oriented pair are not in the manuscript. I marked the clause explicit. If
the names are wanted, it is one cross-reference once the block that introduces them
lands.

## What the proofs need

From the consolidated manuscript:

- Section 3 for the symmetry, the indefinite form and the adjoint for that form. They
  are named and compared. No proof uses a result of that section.
- Section 7 for the two-term complex, its Laplacian and the two odd operators named after
  `eq:index-graded`. The fragment identifies its objects with these. The proofs are
  self-contained.
- Section 8 for the worked example: the vectors and weights of the pair model,
  `eq:pair-invariance-amplitude` and `eq:prepared-schur`.

Standard results used without proof:

- Calculus of operator-valued functions that are differentiable in norm, including the
  product rule and the symmetry of second derivatives.
- Existence and uniqueness of the propagator of a linear differential equation with a
  bounded, piecewise continuous, anti-Hermitian generator, and Liouville's formula for
  its determinant.
- Stokes' theorem and the existence of arbitrarily fine triangulations of a closed
  surface.
- The Chern-Weil theorem, for two statements only: the class of a characteristic form
  does not depend on the connection, and it is the Chern character. Both are cited and
  not reproved. Closedness and integrality of the flux are proved in the fragment.

From TAUM outside the section. Its text cites four sections (`sec:setup`,
`sec:classification`, `sec:algebra`, `sec:polar`) and no single result. The fragment
replaces them by Sections 2, 3 and 7. TAUM says the tangent operators "fall under the
classification already proved". The fragment proves that they are odd and does not use a
classification. No proof here needs the fundamental decomposition, the oriented pair or
the domain results that are being ported separately.

## Where the fragment departs from TAUM

- Theorems, lemmas, definitions and remarks appear as propositions, displayed
  definitions and prose, as the manuscript's preamble requires.
- The curvature is stated for arbitrary vector fields of class C1, with the bracket
  term. TAUM states it for commuting coordinate fields, which is the last display of the
  proposition.
- Frame covariance is proved from the action of the connection on coefficient vectors.
  TAUM quotes a differential-form identity.
- The flux corollary is stated for every finite rank, with the trace. TAUM states rank
  one.
- Holonomy is defined without a frame and for any dimension, and then in a frame. TAUM
  defines it in a frame, for a loop inside one frame domain.

## Added mathematics, to be reviewed as new text

- `prop:moving-unitary`: the direct rotation of two nearby projections, with proof. TAUM
  asserts that the dimensions are locally constant and that the ranges form a bundle.
- The formula for a local frame.
- Compatibility of the connection with the inner product.
- The transport equation `eq:moving-transport` and the intertwining
  `eq:moving-transport-intertwining`.
- The proof that the characteristic forms are closed, from the fact that an odd number
  of derivatives of the projection has no trace.
- The proof of integer flux by a triangulation and determinants of edge transports.
- The formula for the holonomy of a loop bounding a disc, in rank one.
- All of `subsec:moving-pair`: the two tangent maps in the pair model, the curvature, the
  flux of minus one over all settings and one period, the holonomy, and the statement
  that the tangent map of waiting is the coupling of the Schur formula times -i/hbar.
- The scope statement in the opening that the connection is not a dynamical field. The
  outline asks for this distinction.

Two sentences attribute known results and need a ruling on whether they belong: the
name "quantum geometric tensor" for the mixed return, and "geometric phase of a cyclic
evolution" for the holonomy in the example.

## References

Kobayashi and Nomizu, and Husemoller, are copied from TAUM's bibliography. Provost and
Vallee (1980) and Aharonov and Anandan (1987) are new. I could not check either online
from this environment, so journal, volume and pages are from memory and need
verification. The same holds for the section numbers given for Kato's book.

## Names to settle at integration

- `\mathcal L_x` is the moving subspace. TAUM's letter is the manuscript's port space.
- `\Xi` is a frame. TAUM's letter is the pair model's silent directions.
- `q` is a change of frame. TAUM's letter is the ground vector.
- `\mathcal F^\Pi` and `\mathcal A_\mu` are the curvature and the connection matrix, as
  in TAUM. `\mathcal F_2` is the two-particle space of the pair model.
- `B_Y` is the tangent map, as in TAUM. It is an instance of the map `A` of Section 7.
- `\mathsf Z` and `\mathsf S` are a closed surface and a disc. `t_\delta` is the period.
- `X` is the parameter manifold and `\mathcal K` the fixed Hilbert space.

## Not in this block

Outline III.2 also lists compact symmetry averaging and invariant compression (TAUM's
section on the Haar projection). That section is not ported here.

## Checks

`checks/fragments/test_moving_subspace.py` has nine tests. Families of projections are
built with exact first derivatives. Three tests go beyond the displayed identities:

- Transport is compared in three ways: the equation without a frame, the equation in a
  frame, and the limit of successive projections along the path.
- The flux is computed by quadrature and by the edge determinants used in the proof, for
  a bundle of rank two over a torus, with integers -2, 0, -1 and 0 for four choices of
  parameters. The two computations agree, which also tests the signs in the proof.
- A bundle of rank two over the four-sphere has vanishing first Chern form and second
  Chern number of modulus one. This tests the normalization of `eq:moving-second-chern`.
