# `fidelity_energy.tex`: excitation energy and distance from a ground, on the form domain

Sources: TAUM lines 1197 to 1349 ("Vacuum fidelity", "Gap consequence", "Relation to the
transition channel") and LG lines 1926 to 2030, with LG `thm:projective-ground-recovery`,
`prop:no-state-below-ground` and the graph clause of `cor:BRLC-FS-isometry`. Base:
`47a7f56`. Planned home: one subsection at the end of Section 4, after "Grounds selected
by two path Hamiltonians". The preview numbers it 13.1 because the builder continues
from the last section. In place it is 4.4.

    python manuscripts/consolidated/fragments/build_fragment.py fidelity_energy
    python checks/run_all.py -k fidelity_energy
    python consolidation/inventory/proposed/build_fidelity_energy_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger.

## Source results a reader can now find

| Source | Fragment |
|---|---|
| TAUM `thm:reference-fidelity` for any unit vector; LG `prop:ground-fidelity-factorization` on the form domain | `prop:fidelity-energy` (a), `eq:fidelity-energy-factor` |
| TAUM `def:conditional-excitation-scale` and the remark after it | `eq:fidelity-energy-factor` and the paragraph after the proof |
| TAUM `cor:fidelity-gap` with `eq:Eexc-form-domain`; LG `cor:ground-gap-fidelity` | `eq:fidelity-energy-form`, `prop:fidelity-energy` (b), `eq:fidelity-energy-gap` |
| TAUM `eq:FS-gap-bound` | `eq:fidelity-energy-distance`, `eq:fidelity-energy-gap` |
| LG `eq:fidelity-state-decomposition` | `eq:fidelity-energy-decomposition` and the two sentences after it |
| LG `prop:no-state-below-ground`, the spectral and form bounds | sentences after `eq:fidelity-energy-form`; `prop:fidelity-energy` (b) for uniqueness |
| LG `thm:projective-ground-recovery` | `eq:fidelity-energy-ratio` and the text after it |
| LG `cor:BRLC-FS-isometry`, the graph maps | paragraph containing `eq:fidelity-energy-transport` |
| TAUM "Relation to the transition channel" | last paragraph |

## Change in the clause count

I searched the whole ledger for clauses that name fidelity, the gap, the form domain or
the Fubini-Study distance and are not yet explicit. Ten belong to this block: four owed
and six recoverable. `consolidation/inventory/proposed/fidelity_energy_rows.csv` proposes
all ten as explicit and adds one row for a TAUM passage that had none, the subsection
"Relation to the transition channel", keyed by its label `eq:cyclic-reference`. That
passage is an entry of the reading queue, which the new row closes.

If the rows are accepted as they stand, the ledger goes from 496 rows (324 explicit, 78
recoverable, 94 owed) to 497 (335, 72, 90). TAUM goes from 142, 17, 31 to 147, 15, 29,
and LG from 181, 61, 63 to 187, 57, 61.

Clauses near these that the fragment does not carry:

- LG `choice:no-fidelity-lower-bound`, second clause, and the Fourth Yoga clause of
  `prop:no-state-below-ground`. Both are interpretive and planned for Part V.
- TAUM `prop:cycle-positivity`, the two clauses on the form domains of the returns and
  the common two-sector energy. They belong with the fundamental decomposition. The last
  paragraph of the fragment uses the form of `A*A` on `Dom A` and does not replace them.

## What the proofs need

From the consolidated manuscript: `eq:ground-return-energy` and `eq:actual-ground-partner`
with the weak bridge and strong compatibility (Section 4); the Fubini-Study distance
(Section 2); the graph maps (Section 3); the two paths (Section 4) and the site state
(Section 5); the centered determinations (Section 5); the number operator (Section 12);
the closed-range criterion (Section 7); `eq:prepared-schur` (Section 8); R2.

Used without proof: the spectral theorem for a self-adjoint operator bounded below, in
the form that `H - E_0` has a positive square root whose domain contains `Dom H`; that
the closed form of `A*A` is the squared norm of `A` on `Dom A`; the values of the zeta
function.

## Mathematics added to the sources

Offered for review as additions.

1. Part (c): for bounded `H` the scale `Z` is at most the width of the spectrum, so the
   excitation energy is bounded above by the width times one minus the fidelity.
2. The three examples of `Z`: the two paths, and a state of the number operator in the
   form domain and outside the operator domain, with its values in terms of the zeta
   function.
3. The sequence showing that a fidelity tending to one does not bound the energy of an
   unbounded Hamiltonian, set against the converse that a gap gives.
4. The transport statement `eq:fidelity-energy-transport`, and the energy `1/16` of the
   image of the coarse ground on the fine path.
5. The reading of the gap hypothesis for `A*A` as the closed-range bound.

## Sentences that interpret

- "The number `Z(ψ)` is the mean excitation energy of `η_ψ`, the state that an ideal
  projective ground test leaves when `ψ` fails it." TAUM calls `Z` the conditional
  excitation scale. I named the instrument because Section 6 distinguishes an effect
  from the instrument that realizes it.
- "The polar map relates the fidelities of the two sides and not their energies."
- "A reference vector with nonzero image is in general not a ground of `A*A`."

## Names

`F_x(ψ)` for the fidelity, so that the `F_x` of Section 4 is `F_x(r)`. `d_FS` for the
Fubini-Study distance. `E_exc`, `η_ψ`, `Z(ψ)` as in TAUM. `Δ` for the gap, as in Section
4. `E_max` for the top of a bounded spectrum. `φ_k` for the vectors of the unbounded
example and `c`, `η` for the two components of a vector.

## Left for integration

- The sentence of Section 4 that follows `eq:ground-return-energy` gives the gap bound
  for the partner. It could point forward to this subsection.
- Placed at the end of Section 4, the subsection refers forward to Sections 5, 7, 8 and
  12 for its examples and for the operators named in its last paragraph. In place the
  manuscript builds to 66 pages without warnings.
- `checks/fragments/test_fidelity_energy.py` has seven tests. On this branch the suite
  gives 152 passed and 2 skipped, and 154 with `--full`.
- The check cites existing `UM` labels only. The fragment's labels can be added once
  they are in the manuscript.
