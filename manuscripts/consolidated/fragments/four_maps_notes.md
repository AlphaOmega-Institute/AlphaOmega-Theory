# `four_maps.tex`: TAUM's canonical realizations

Source: TAUM lines 3628 to 4048, with 9 formal entries and 7 supplementary equations in
the ledger. Base: `0cd012d`. Planned home: outline IV, checked mathematical and physical
realizations.

    python manuscripts/consolidated/fragments/build_fragment.py four_maps
    python checks/run_all.py -k four_maps
    python consolidation/inventory/proposed/build_four_maps_rows.py

Nothing here changes `Universal_Model.tex`, the shared ledger or the existing checks.

## Source results a reader can now find

| TAUM | Fragment |
|---|---|
| opening of the section and the closing table | opening of the section, `subsec:four-maps-table` |
| `ex:scalar-channel` | `subsec:four-maps-scalar`, `eq:four-maps-scalar`, `eq:four-maps-scalar-response` |
| `ex:split-excitation-response` and the paragraph after it | `subsec:four-maps-three-state` |
| `ex:translation-channel`, `eq:no-translation-gap`, `eq:free-response` | `subsec:four-maps-translation` |
| `thm:oscillator-channel`, `cor:oscillator-heat-index` | `prop:four-maps-oscillator` |
| `thm:oscillator-zeta`, remark at 3901 | `prop:four-maps-mellin` and the paragraph after it |
| `eq:oscillator-response` | `eq:four-maps-oscillator-port`, `eq:four-maps-oscillator-response` |
| `ex:quasiperiodic-channel`, remark at 4015 | `eq:four-maps-interval-operator`, `eq:four-maps-interval-spectrum` and the text after them |
| `eq:unsmoothed-point-mass`, `eq:theta-mass`, `eq:quasiperiodic-smoothed-response`, `eq:normalized-smoothed-point` | `eq:four-maps-interval-point` to `eq:four-maps-interval-response` |

## Change in the clause count

Before proposing rows I searched the whole ledger for clauses that name these results.
All of them lie in this section: 20 clauses, of which 1 is explicit, 3 are recoverable
and 16 are owed. In `consolidation/inventory/proposed/four_maps_rows.csv` the 19 that
were not explicit are explicit and the other is unchanged. Two further rows cover source
passages that had none: the opening of the section, which is an entry of the reading
queue, and the paragraph after the three-state example. If the rows are accepted as they
stand, the ledger goes from 493 rows (254 explicit, 108 recoverable, 131 owed) to 495
(275, 105, 115), and TAUM from 188 (121, 20, 47) to 190 (142, 17, 31). Review of
`b301a28` accepted the block in substance.

The rows file also holds one owed row agreed in that review, for LG. After
`eq:vacuum-symmetry-phase` LG says that the phase freedom of a unique ground ray is not a
model of individual differences among sentient beings. Section 11 carries the
mathematics of that equation and not this sentence, so the sentence is recorded as owed,
planned for Part V. With it the ledger would stand at 496 rows (275, 105, 116).

TAUM's closing table (lines 4025 to 4048) is carried in `subsec:four-maps-table`. That
passage has no source label, so no row can be keyed to it. Its entry in the reading queue
can be closed by hand at integration.

## What the proofs need

From the consolidated manuscript:

- Section 7 for the differential, the two odd operators, the Laplacian, the criterion
  for a closed range stated after `prop:index-harmonic-pairing`, the heat-trace formula
  `eq:index-heat`, and the unilateral shift.
- Section 9 for the canonical response (`cor:passive-sources`), the reached space
  (`prop:passive-two-spans`), uniqueness (`prop:passive-uniqueness`), the sector
  proposition (`prop:passive-opposite-sectors`) and the smoothed ports
  (`prop:passive-heat-distributional`, `prop:passive-heat-discrete`).
- Section 1 for the count of indirect directions (`cor:finite-response-defect`) and the
  first form of the three-state example.

Standard results used without proof:

- The Fourier transform is unitary and turns the derivative into multiplication, and the
  derivative is self-adjoint on its Sobolev domain.
- The ladder operators of the oscillator are closed and adjoint to each other.
- The first derivative with a twisted condition at the ends is self-adjoint, and its
  exponentials form an orthonormal basis.
- Euler's integral for the gamma function, the series for the zeta function, and the
  exchange of a sum and an integral under absolute convergence.
- The Poisson summation formula.

From TAUM outside the section. Its text cites five results: the closed-range criterion,
the index formulas, minimal uniqueness, the finite response count and the heat-smoothed
port. All five are explicit in the manuscript.

## Where the fragment departs from TAUM

- TAUM's map `B` is the manuscript's `A`, and its two odd operators are written
  `d_A + d_A^*` and `d_A^* - d_A`, as in Section 7.
- Examples appear as subsections. The oscillator results are two propositions.
- Planck's constant is explicit in the transfer formula of the three-state example.
- The half-length of the interval is `\Lambda` and its variable is `x`. TAUM's letters
  are the manuscript's angular momentum and polar part.
- The port vectors of the translation and oscillator examples are `\eta` and `\eta_r`.
  TAUM's letter is the manuscript's ground vector.
- The closing table uses the manuscript's `tabular`, with a fourth feature added to the
  translation row: the map is not Fredholm.

## Added mathematics, to be reviewed as new text

- The scalar map at the number zero.
- For translations: the range is dense, the reduced cohomology vanishes and the map is
  not Fredholm; the second form of the response, as a measure on the spectrum of the
  return; and the statement that a one-dimensional port is never minimal there, with an
  explicit vector orthogonal to the reached space.
- For the oscillator: the polar decomposition, whose isometry is the shift of Section 7,
  and the remark that the heat-trace formula applies here although it did not apply to
  the shift; the identity for every summable function of the number operator, of which
  TAUM's two heat traces are cases; and the pole of the completed factor at one.
- For the interval: the exponent one half of the scale on which evaluation at a point is
  bounded; the second form of the mass by Poisson summation, with its leading term and
  exponentially small corrections; and minimality of the smoothed pair, so that the
  response determines the length and the twist.
- The sentence in the three-state example that the inequality on the bias makes the
  ground unique.

Three comparisons were drawn out in review of `b301a28` and are now in the text:

- Three-state example. The port and the direction it misses are the same at every bias.
  The bias decides whether the evolution carries the prepared excitation into that
  direction: it is the indirect space when the lines are split and permanently invisible
  when they coincide.
- Oscillator and line. Both use a port of rank one. The oscillator port reaches every
  level. On the line a port whose Fourier transform vanishes nowhere still misses part
  of the space, because each energy has two momentum directions.
- Interval. The smoothed mass is the same for opposite twists at every smoothing scale.
  The first-order response tells them apart, unless the twist is 0 or pi. Neither
  depends on the point.

Two corrections from that review, both in the Mellin passage:

- The calculation uses the positive integers as eigenvalues of multiplicity one. An
  earlier version said it used only the set of eigenvalues. Two copies of the oscillator
  have the same set and twice the trace.
- Integration term by term is justified by the sum of the integrals of the absolute
  values. An earlier version summed the absolute values of the integrals.

The opening sentence now says densely defined closed map, since it uses the adjoint.

One edit agreed in review belongs to the manuscript and is left for integration: in
Section 11 the exchange of the two orbitals "fixes the vacuum", which should name the
empty Fock state, because the same section uses vacuum for the ground of a Hamiltonian.

## Names to settle at integration

- `a`, `a^*` and `N` are the ladder operators and the number operator. `N` also counts
  points in Section 9.
- `\mathsf S` is the unilateral shift.
- `\Gamma(s/2)` in `prop:four-maps-mellin` is Euler's gamma function. Everywhere else
  `\Gamma` is a port. The proposition says so.
- `D_x` and `D_{\Lambda,\vartheta}` are the two first-order operators, `\lambda_n` and
  `e_n` their eigenvalues and eigenfunctions. `e_n` is also the oscillator basis.
- `\phi_n` are the low-frequency unit vectors, `\eta` and `\eta_r` the port vectors,
  `b_\tau` the smoothed point vector and `Z_\tau` its squared norm. `b` is also the
  number of the scalar example and the prepared excitation of the pair model.
- No new references.

## Checks

`checks/fragments/test_four_maps.py` has eight tests.

- The oscillator map is tested exactly, between the first M and the first M + 1 levels,
  where the two returns are diagonal with entries 1 to M and 0 to M.
- The translation example is tested through its Fourier multipliers: the response of a
  Gaussian against a closed form, the two forms of the response against each other, and
  the vector orthogonal to the reached space.
- The Mellin transform is integrated numerically for real and complex s and compared
  with the completed zeta factor. Its pole at one is checked, the doubled trace of two
  copies is checked, and the sum of the integrals of the absolute values is compared
  with its closed form.
- The twisted interval is tested through its eigenfunctions and through a Hermitian
  difference operator with the same condition at the ends. The two forms of the mass
  agree to eleven digits.
- One test applies the common relations to all the maps and shows the gap of the
  translation generator closing on circles of growing length.
- The run uses no numerical integration that raises a warning. The check runner fails if
  a warning is raised inside a test, because the warning adds a name to the module while
  the runner iterates over it.
