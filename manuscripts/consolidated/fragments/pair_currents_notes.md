# `pair_currents.tex`: the one-body currents of the pair model and of the trap

Source: LG lines 7229 to 7830, the parts of "A finite pair-transfer realization" and "An
orbital source preserving the full ground" that the manuscript did not yet state. Base:
`067ce96`. Planned home: outline IV, as a new Section 13 after "Four maps with one
algebra". It could also be divided between Sections 1, 8 and 11, where the pair model and
the trap already appear. I wrote it as one section because the three currents organize
all of it.

    python manuscripts/consolidated/fragments/build_fragment.py pair_currents
    python checks/run_all.py -k pair_currents
    python consolidation/inventory/proposed/build_pair_currents_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger.

## Source results a reader can now find

| LG | Fragment |
|---|---|
| `eq:pair-parent`: the coupling in Singh's convention | sentence after `eq:pair-currents-hamiltonian` |
| `eq:pair-mixer-selection`: currents, `T² = 2P_c*P_d`, `Xg = m`, `Yg = 0`, `Zg = e` | `eq:pair-currents`, `eq:pair-currents-square`, `eq:pair-currents-action` |
| the same passage: the rotation as `exp(-2iθY)`, the drives that keep the ground | `eq:pair-currents-rotation`, `eq:pair-currents-drive` |
| "Which operations the pair component supplies": the parity restriction | second paragraph of `subsec:pair-currents-drives` |
| `eq:pair-channel`: restriction and compression of the input generator | `eq:pair-currents-compression` |
| `eq:pair-readout-action`, `eq:pair-defect-measurement` | `eq:pair-currents-readout`, `eq:pair-currents-zero` and the text after it |
| `eq:pair-general-mixer`, `eq:pair-mixer-ground-certificate` | `eq:pair-currents-mixer` to `eq:pair-currents-quaternion` |
| `prop:pair-invariant-mixer` and the paragraph after it | `prop:pair-currents-invariant-mixer` and the paragraph after it |
| `prop:pair-ground-defect`: cyclicity, and the paragraph after its proof | first paragraph of `subsec:pair-currents-algebras` |
| `eq:pair-two-frequency-response`, first line of `eq:pair-completed-response` | `eq:pair-currents-correlation` and the text after it |
| `eq:orbital-parent`, `eq:orbital-bias-parameters`: gauge, minimal coupling, spin | `eq:pair-currents-minimal-coupling` and the text after it |
| `eq:orbital-derived-current` | `eq:pair-currents-angular-momentum`, `eq:pair-currents-shell` |
| `eq:orbital-quadrupoles` | `eq:pair-currents-quadrupoles`, `eq:pair-currents-shell` |
| `eq:orbital-low-carrier`: the third basis vector | sentence before `eq:pair-currents-orbital-preparation` |
| `eq:orbital-typed-preparation` | `eq:pair-currents-orbital-preparation`, `eq:pair-currents-orbital-compression` |
| `eq:orbital-filtered-readout`: generated algebras | paragraph before `eq:pair-currents-orbital-correlation` |
| `eq:orbital-low-response`: ground correlation | `eq:pair-currents-orbital-correlation` |
| `eq:orbital-unfiltered-response`: departure probability | `eq:pair-currents-unfiltered` |

## Change in the clause count

Before writing rows I searched the whole ledger for clauses that name the mixers, the
currents, the quadrupoles, the gauge or the spin term. All of them lie in LG lines 7200
to 7830, which hold 57 clauses.

`consolidation/inventory/proposed/pair_currents_rows.csv` has 28 rows, all proposed as
explicit: 11 that were owed and 17 that were recoverable. With them every clause of those
two LG sections is explicit. If the rows are accepted as they stand, the ledger goes from
496 rows (296 explicit, 95 recoverable, 105 owed) to 496 (324, 78, 94), and LG from 153,
78, 74 to 181, 61, 63. No row is added.

## What the proofs need

From the consolidated manuscript: the pair model of Section 1 (`g`, `e`, `m`, `W`, `H_s`,
`A_θ`, `b`, `d`, `eq:pair-invariance-amplitude`, `eq:pair-ambient-zero`); the rotation
`κ(φ)` of Section 11; the trap of Section 8 (`eq:compensated-orbital-parent`, the shell,
the exact branches, `Z_0`, `e_Q`, `C_0`); strong compatibility
(`prop:ground-polar-compatibility`); R2.

Used without proof: the canonical anticommutation relations; that a `2×2` unitary of
determinant one has the form `eq:pair-currents-quaternion`; that a one-body unitary acts
on the coefficient matrix of a two-particle state by a congruence; the Gaussian integrals
behind the three quadrupole matrix elements.

The convention for the coupling is taken from Singh's Eq. (3.11). I checked the wording
against a copy of the paper outside the repository: the constant is `J = 4κ` for the
scalar trace and matrix mode specified there, and `4κ/N` for the other mode named there.

## Mathematics added to the source

Offered for review as additions.

1. The three vectors `g`, `e`, `m` carry the spin-one representation of the currents, and
   `W` is three invariant vectors. LG gives the action on `g` only.
2. The pair Hamiltonian in the currents, `H_s = -sZ² - j(X² - Y²)` on the two-particle
   space, and `[Y, H_s] = i(j - s)(XZ + ZX)`. LG gives the commutator applied to `e`.
3. The imbalance as a rotated current, `Z_θ = u_θ Z + v_θ X`, from which LG's three
   actions follow.
4. The explicit image `eq:pair-currents-mixer-image` behind LG's norm formula.
5. Part (c) of the proposition as a statement with its proof. LG states the determinant
   condition in the paragraph after its proposition.
6. The origin of the quaternion form: the one-body unitaries that keep the ground ray
   are the phases times the unitary symplectic group of the form `ε ⊕ ε`, with a
   ten-dimensional Lie algebra. This is new.
7. The first-order position of the lowest triplet component with a spin moment. LG says
   only that a spin Zeeman term must be included in the test.
8. The coupling out of the prepared plane in the trap, `-2iħΩ_B m_λ`, and the contrast
   between the two models: the pair input generator is a restriction that differs from
   the compression, and the trap input Hamiltonian is a compression that is not a
   restriction.

## Sentences that interpret

These are the ones to read as proposals.

- "The answer separates what the pair Hamiltonian supplies from what the apparatus has
  to add." and "Singh's canonical component supplies `H_0` and its ground. It does not
  supply the drive." LG says that Singh's example leaves those interventions to a further
  physical model.
- "A charged mediator or reference system could give an invariant joint interaction. It
  would then belong to the preparation and has to be kept in the model of the
  measurement." This follows LG.
- "Cyclicity does not select the grounds. The Hamiltonians do." This follows LG.
- "*The Luminous Ground* takes this compression as the input Hamiltonian. That is a
  choice in addition to the prepared plane and the output dynamics." The consolidated
  manuscript already names LG in this way in Section 8.

## Not carried

LG's last paragraph of the pair section says that adding `W` to the output while keeping
only the current algebra would destroy cyclicity of `g`. It has no ledger row. The
fragment states that the currents annihilate `W`, from which this follows, and does not
draw the conclusion.

## Names

`T`, `X`, `Y`, `Z` for the currents and `T_K`, `Y_K` for a general mixer, as in LG.
`Π_c` for the parity of `N_c`. `R_c`, `R_d` for the two representations, with group
elements `h`. `Ω` and `ϵ` for the antisymmetric forms and `𝖱` for a one-particle unitary;
I avoided `w`, which is an entry of `K`. `m_λ` for the third exact vector. `p_⊥(t)` for
the weight outside the prepared plane, as in Section 2. `g_s` and `M_S` for the spin
factor and projection.

## Left for integration

- `checks/fragments/test_pair_currents.py` has nine tests. The statements of LG in this
  block were already checked in `checks/pair/` and `checks/orbital/`; the new file checks
  the additions and the forms in which the text states the results. On this branch the
  suite gives 143 passed and 2 skipped, and 145 with `--full`.
- The check cites existing `UM` labels only. The fragment's labels can be added to its
  `LABELS`, and to those of the existing pair and orbital checks, once they are in the
  manuscript.
- The section goes before the line that opens the small type for the bibliography.
