# Fragments for review

A fragment is a section written for the consolidated manuscript and kept in its own file
until it is integrated. It is body text only. It uses the preamble of
`../Universal_Model.tex` and refers to labels of that file. Nothing here changes
`Universal_Model.tex`, the shared ledger or the existing checks.

    python manuscripts/consolidated/fragments/build_fragment.py passive_response
    python checks/run_all.py            # includes checks/fragments/
    python consolidation/inventory/proposed/build_passive_response_rows.py

The first command typesets the fragment with the main file's preamble, real numbers for
references into Sections 1 to 8, and the main bibliography plus
`passive_response_bibliography.tex`. It needs `pdflatex`. The built PDF is committed
for reading.

## `passive_response.tex`: TAUM's passive-response family

Source: TAUM lines 2845 to 3628, 26 formal entries. Base: `1419df8`.

### Source results a reader can now find

| TAUM | Fragment |
|---|---|
| `ass:selfadjoint-response`, `def:boundary-response` | opening of the section, `eq:passive-response` |
| `thm:response-Herglotz` | `prop:passive-herglotz` |
| `def:Nevanlinna-kernel`, `thm:positive-response-kernel` | `eq:passive-nevanlinna-kernel`, `prop:passive-kernel` |
| `cor:canonical-transition-response` | `cor:passive-sources` (a) |
| `cor:Hilbert-block-response` | `cor:passive-sources` (b) |
| `rem:positivity-boundary` | paragraph after the corollary |
| `def:Schur-response`, `lem:Cayley-denominator`, `thm:Schur-contractivity` | `eq:passive-schur-response`, `prop:passive-cayley` |
| `def:dBR-kernel`, `thm:positive-Schur-kernel` | `eq:passive-dbr-kernel`, `prop:passive-schur-kernel` |
| `def:compressed-spectral-measure`, `thm:spectral-temporal-response` | `eq:passive-spectral-measure`, `prop:passive-spectral-temporal` |
| `def:response-generated-subspace` | `prop:passive-two-spans` and the definition of a minimal pair |
| `thm:canonical-minimal-reduction` | `prop:passive-two-spans`, `prop:passive-minimal-reduction` |
| `thm:minimal-uniqueness` | `prop:passive-uniqueness` |
| `cor:finite-response-defect` | Corollary 1.3 of the draft, unchanged, with three remarks in `subsec:passive-minimal` |
| `prop:permanent-port-invisibility` | Proposition 1.1 of the draft, unchanged, with the paragraph on finite windows |
| `prop:opposite-sector-response` | `prop:passive-opposite-sectors` and the two paragraphs after it |
| `ass:distributional-trace`, `thm:heat-smoothed-distributional-port` | `subsec:passive-heat`, `prop:passive-heat-distributional` |
| `ass:discrete-response-spectrum`, `thm:heat-smoothed-channel` | `prop:passive-heat-discrete` |
| `rem:heat-channel-normalization` | last paragraph |

### Change in the clause count

The ledger at `1419df8` has 40 clauses for these 26 entries: 9 explicit, 8 recoverable,
23 owed. Thirteen entries were entirely owed and three entirely explicit. In
`consolidation/inventory/proposed/passive_response_rows.csv` all 40 are explicit, and 5
further clauses are separated out from prose that follows TAUM's statements (scope of
the definition, enlarging the port, finite windows, the origin of matched measures, and
the cautions about crossing and return).

If the rows are accepted as they stand, the ledger goes from 479 clauses (182 explicit,
117 recoverable, 180 owed) to 484 (218, 109, 157), and TAUM from 174 (50, 28, 96) to
179 (86, 20, 73). These are proposals. The labels and quotes point into the fragment,
so the rows cannot enter `migration.csv` until the text has a permanent home.

One row needs a decision. `rem:positivity-boundary` is carried as a statement of scope
with status `open question`. The fragment asserts that a generator self-adjoint for the
indefinite form is in general not Hilbert-self-adjoint and gives no witness, as in TAUM.
`checks/pending/test_port_calibration.py` has one, and the signed-response section is
the natural place for it.

### What the proofs need

From the consolidated manuscript:

- Proposition 1.1 (`prop:preparation-split`) for the reached space, the indirect space
  and permanent invisibility. Corollary 1.3 (`cor:finite-response-defect`) and
  `prop:semibounded-waiting` are cited and not reproved.
- Section 7 for the fact that the return of a closed densely defined map is positive and
  self-adjoint. Part (a) of the corollary depends on that sentence.
- `eq:prepared-schur` in Section 8 for part (b), with self-adjoint diagonal blocks and
  bounded coupling as stated there.
- The pair model and the trap of Section 8, as examples only
  (`eq:pair-two-setting-ranges`, `subsec:orbital-apparatus`).
- Section 3 is named once, for the indefinite form in the statement of scope. No proof
  uses it.

Standard results used without proof. Kato, already in the draft's bibliography, covers
the first four, and the fragment cites it only at the perturbation step:

- The spectral theorem and functional calculus: norm analyticity and adjoint symmetry of
  the resolvent, the resolvent identity, the spectral integrals for the resolvent and
  the unitary group, the Laplace integral for the resolvent, and the strong limit of
  the scaled resolvent at large imaginary argument.
- A closed subspace whose projection commutes with every resolvent reduces the operator,
  and the restriction to a reducing subspace is self-adjoint with the restricted
  resolvent.
- Unitarily intertwined resolvents give intertwined operators with matching domains.
- A bounded symmetric perturbation of a self-adjoint operator is self-adjoint on the
  same domain (cited to Kato).
- A positive kernel is the reproducing kernel of a unique Hilbert space (cited to
  Aronszajn). The fragment states this consequence and does not reprove it.
- The scale of spaces of a self-adjoint operator and its duality with the carrier as
  pivot, for the distributional trace.

From TAUM, outside the passive-response section. Its text cites four outside labels.
`sec:setup` and `sec:effective-dynamics` are replaced by Sections 1, 7 and 8 of the
draft. `ass:block-generator` is needed only on its Hilbert branch, which Section 8
states. Its signed branch is still owed and no proof here uses it.
`eq:realization-datum` is a later list of the data of a realization. The fragment keeps
the one sentence that refers to it and does not port the list. None of the proofs uses
the objects of TAUM's fundamental decomposition beyond the two sectors and the map
between them, which Section 1 supplies.

### Where the fragment departs from TAUM

- Theorems, lemmas, definitions and remarks appear as propositions, displayed
  definitions and prose, because the draft's preamble defines only proposition,
  corollary and conjecture.
- Planck's constant is explicit, as in the draft.
- The Cayley-denominator lemma is the first statement of `prop:passive-cayley`.
- The block corollary is stated for the retained and complementary spaces of
  `eq:prepared-schur` and not for TAUM's two sectors.
- The canonical response uses the closed map of the draft. TAUM states it for its
  transition channel.
- The heat port is defined as the adjoint of a bounded operator on the carrier. This
  gives TAUM's formula and avoids a separate extension of the Gaussian to the negative
  scale.

Added mathematics, to be reviewed as new text:

- A full proof that the resolvent span and the time span agree for arbitrary spectra.
- In the uniqueness proof, the continuity step at equal points of the resolvent
  identity, which TAUM omits, next to the one at conjugate points, which it mentions.
- A proof that the count of Corollary 1.3 is the dimension of the space of relations
  among the weight ranges.
- A pair-model example in which enlarging the port lowers the count from one to zero.
  The check also has a case where it rises.
- An explicit overlapping example for the caution that the analysis line need not lie
  in the reached space.
- The reason an antiunitary symmetry still gives equal measures.

### Names to settle at integration

- `\Theta` is the unitary of `prop:passive-uniqueness`. TAUM's letter clashes with the
  draft's `W` for the silent directions of the pair model.
- `R(z)` and `R_j(z)` are resolvents, and `R_k` without an argument is a weight of
  Corollary 1.3. TAUM has the same pairing.
- `S_\Gamma` is the Schur response. Section 2 uses `S` for the initial support.
- `C_\Gamma(t)` and `C_0(t)` are time-domain functions. Section 8 uses `C` for the
  coupling, which also appears in part (b) of the corollary.
- `\mathcal X` is the carrier, since two of the uses take it to be the input sector or
  a general Hilbert space and not `\mathcal H_\Phi`.
- Three references are new: Donoghue, Aronszajn, and de Branges and Rovnyak, in
  `passive_response_bibliography.tex`, copied from TAUM's bibliography.

### Checks

`checks/fragments/test_passive_response.py` has nine tests. Each constructs its own
operators and tests the stated identities on them. Its `LABELS` name TAUM labels only,
because the fragment's labels are not yet in `Universal_Model.tex`.
