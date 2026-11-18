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

Review of `5e937be` accepted all 27 rows, and a trial integration validated at these
counts. `thm:tangent-realization` is explicit: the fragment gives the two oriented sums
as `d + d*` and `d* - d` in the notation of Section 7, and substituting TAUM's own
definitions of its oriented pair gives the same two sums.

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
- `eq:moving-holonomy-flux`: the determinant of the holonomy of a loop bounding a disc,
  for every finite rank, and the paragraph that separates curvature at a point from
  holonomy of a loop.
- `eq:moving-scalar-form`: the scalar form of the mixed return in finite rank, where it
  vanishes, and the hypothesis it needs in infinite rank.
- All of `subsec:moving-pair`. The moving subspace there is the prepared plane inside the
  output space, and the splitting into the two sectors stays fixed. The example gives the
  two tangent maps and their imaginary ratio, the curvature, the Chern number of minus
  one for the sphere of distinct prepared planes (settings from 0 to pi/4, waiting times
  modulo one period, each end circle collapsed to a point), the vanishing integral over
  the plane's full period pi/2, the holonomy of one waiting period, and the statement
  that the tangent map of waiting is the coupling of the Schur formula times -i/hbar.
- The scope statement in the opening that the connection is not a dynamical field. The
  outline asks for this distinction.

Two attributions were kept in review, each with a restriction. The mixed return is
called the quantum geometric tensor of Provost and Vallee only for a family of rays. The
geometric phase of a cyclic evolution, after Aharonov and Anandan, is the phase on the
excitation in the example, not the holonomy of the plane as a whole.

## References

Kobayashi and Nomizu, and Husemoller, are copied from TAUM's bibliography. Provost and
Vallee (1980) and Aharonov and Anandan (1987) are new. They could not be checked online
from the environment in which the fragment was written. The review of `5e937be` checked
both, and the references to Kato's book for nearby projections and for transport.

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
built with exact first derivatives. Inverse square roots of positive matrices use
`lib.linalg.hermitian_power`, which stays in double precision. With SciPy 1.14.1,
`scipy.linalg.sqrtm` returned an extended-precision array that `numpy.linalg` refused, in
two tests here and in one of the response checks. The suite passes with SciPy 1.14.1 and
with the current release. Three tests go beyond the displayed identities:

- Transport is compared in three ways: the equation without a frame, the equation in a
  frame, and the limit of successive projections along the path.
- The flux is computed by quadrature and by the edge determinants used in the proof, for
  a bundle of rank two over a torus, with integers -2, 0, -1 and 0 for four choices of
  parameters. The two computations agree, which also tests the signs in the proof.
- A bundle of rank two over the four-sphere has vanishing first Chern form and second
  Chern number of modulus one. This tests the normalization of `eq:moving-second-chern`.
