# `symmetry_average.tex`: TAUM's compact symmetry averaging and invariant compression

Source: TAUM lines 1078 to 1187, three formal entries. Base: `91c3b9b`. Planned home:
outline III.2. Review of `b08bad6` decided to keep it as Section 11, next to Section 10:
that section asks what happens when the prepared space moves, and this one asks when a
symmetry singles out a space on which the dynamics can be computed exactly.

    python manuscripts/consolidated/fragments/build_fragment.py symmetry_average
    python checks/run_all.py -k symmetry_average
    python consolidation/inventory/proposed/build_symmetry_average_rows.py

Nothing here changes `Universal_Model.tex`, the shared ledger or the existing checks.

## Source results a reader can now find

| TAUM | Fragment |
|---|---|
| opening of the section, `eq:unitary-representation`, `eq:Haar-projector` | `eq:symmetry-average`, `eq:symmetry-invariant-space` |
| `thm:Haar-projector`, first part | `prop:symmetry-projection`, `eq:symmetry-absorption` |
| `thm:Haar-projector`, compatibility with the grading | `prop:symmetry-compatibility` (a), `eq:symmetry-sector-sum` |
| `cor:finite-group-average` | `eq:symmetry-finite-average` and the sentence after the proof |
| remark at 1177 | last paragraph of `subsec:symmetry-compatibility`, and `subsec:symmetry-pair` |
| `thm:architecture`, compact-symmetry clause | `prop:symmetry-projection`, `prop:symmetry-compatibility` (a) |

## Change in the clause count

Before proposing rows I searched the whole ledger for clauses that name these results.
There are five, all owed: the four of this section and the compact-symmetry clause of
TAUM's closing theorem. In `consolidation/inventory/proposed/symmetry_average_rows.csv`
all five are explicit, and one further row covers the opening of the section, an entry
of the reading queue (TAUM 1078 to 1085). If the rows are accepted as they stand, the
ledger goes from 492 rows (247 explicit, 109 recoverable, 136 owed) to 493 (253, 109,
131). With this block every clause of TAUM's closing theorem that is not about the
signed reduction or the two odd branches is explicit. Review of `b08bad6` accepted all
six rows.

## What the proofs need

From the consolidated manuscript:

- Section 3 for the adjoint for the indefinite form, in one identity of part (a).
- Section 2 for the support projections, named in part (b).
- `prop:passive-two-spans` for part (c): the reached space is the closed span of the
  resolvent vectors, and a projection that commutes with every resolvent reduces the
  operator.
- `prop:preparation-split` for the remark on permanent invisibility.
- Section 8 and `subsec:moving-pair` for the worked example, and `prop:shell-accuracy`
  and `prop:passive-opposite-sectors` for the two sentences about the trap.

Standard results used without proof:

- Existence of the Haar probability measure of a compact group and its invariance under
  left translation, right translation and inversion.
- A bounded sesquilinear form defines a bounded operator.
- A decreasing sequence of projections converges strongly to the projection onto the
  intersection of their ranges.

From TAUM outside the section: nothing. Its text cites no other label.

## Where the fragment departs from TAUM

- The average is defined by its matrix elements. TAUM uses a vector-valued integral. One
  sentence identifies the two: strong continuity makes each orbit Bochner integrable, and
  the inner product of its integral with a test vector is the displayed formula.
- The representation is written `\kappa` and a group element `a`. TAUM's letters are the
  manuscript's density operator and ground vector.
- TAUM's theorem, corollary and remark appear as one proposition, a sentence and a
  paragraph with an example.

## Added mathematics, to be reviewed as new text

- In `prop:symmetry-projection`, the statement that the average maps every closed
  invariant subspace into itself. Parts (b) and (c) below rest on it.
- The paragraph on fixed vectors and invariant states. A pure state is invariant exactly
  when the group multiplies its vector by a character, and `eq:symmetry-character-projection`
  is the projection onto the vectors of a given character. The invariant projection is
  the case of the trivial character.
- The paragraph on subgroups: a subgroup fixes at least as much; the cube rotations on
  homogeneous quartic polynomials fix a plane where all rotations fix a line; an
  increasing sequence of closed subgroups with dense union gives projections that
  decrease strongly to the invariant projection of the whole group.
- Part (b) of `prop:symmetry-compatibility`: a closed covariant map takes invariant
  vectors to invariant vectors, its adjoint is covariant, and its support projections
  commute with the representations.
- Part (c): the invariant space reduces a Hamiltonian that commutes with the group, and
  a port into the invariant space has its reached space inside it and the response of
  the compressed pair. This is the "invariant compression" of TAUM's title, which its
  section names and does not state.
- The two sentences on the trap: parity as the rotation by pi, and the eigenspaces of
  the angular momentum as ranges of twisted averages.
- All of `subsec:symmetry-pair`, as revised in review.
  - The exchange of the two orbitals has an invariant projection equal to the prepared
    projection at the setting pi/4.
  - The orbital rotations are the preparation rotation. The whole circle commutes with
    the Hamiltonian only at s = j. The quarter turns commute with it at every coupling
    and have the same invariant projection, of rank four, so a port into that space has
    the compressed response although the continuous symmetry is broken. The preparation
    is not such a port, since the prepared excitation is orthogonal to that space.
  - For a cyclic subgroup, n is its order as a group of rotations of the orbitals. Its
    action on two-particle vectors can have a kernel.
  - The vectors fixed by the orthogonal group of the two orbitals form the ground line.
    Every ray of the silent space is also an invariant pure state, because the exchange
    only changes its sign. A fixed vector and an invariant state are different
    requirements. An earlier version called the ground the only fixed state, which was
    wrong for states.
- The scope statement in the opening. The outline asks that a global symmetry be kept
  apart from a spacetime-local theory, and TAUM's table of established layers says the
  average "yields an invariant projection, not a preferred dynamics".

One sentence carries an item from TAUM's list of open extensions (lines 4199 to 4247):
the convergence of the invariant projections of finite groups to that of a continuous
group, with the compressed Hamiltonians and responses. That passage has no source label
and stays in the reading queue for its other items.

## Names to settle at integration

- `\kappa` is the representation, `a` and `b` are group elements, `da` is the Haar
  probability measure.
- `\Pi_G` and `\mathcal K^G` are the invariant projection and the invariant space, as in
  TAUM. `\Pi_G^\Psi` and `\Pi_G^\Phi` are those of the two sectors.
- `H^G` is the restriction of a Hamiltonian to the invariant space.
- `\sigma` is the exchange of the two orbitals of the pair model. Sections 2 to 7 use the
  same letter for a density operator.
- `\Pi_\sigma` and `\Pi_{\rm rot}` are the two invariant projections of the example.
- No new references.

## Checks

`checks/fragments/test_symmetry_average.py` has seven tests.

- The average is compared with the projection onto the joint fixed space for
  permutation groups, a circle group and SU(2) on two, three and four spins, where the
  ranks are 1, 0 and 2.
- The cube rotations act on polynomials of degree 0 to 6 with 1, 0, 1, 0, 2, 0, 3
  invariants, against 1, 0, 1, 0, 1, 0, 1 for all rotations.
- The pair-model test builds the action of a change of the four modes on the
  two-particle space and tests every statement of the example. It tests separately that
  the circle commutes with the Hamiltonian only at s = j and that the quarter turn and
  the invariant projection of the circle commute with it at every coupling.
- A test of characters shows invariant pure states whose vectors the average annihilates.

The abstract is not changed on this branch. Review of `b08bad6` supplied two sentences
for it, on the response result of Section 9 and on the geometric phase of Section 10,
to be added when the section is integrated and the manuscript PDF is rebuilt.
