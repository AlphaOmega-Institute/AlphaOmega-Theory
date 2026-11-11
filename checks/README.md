# Verification checks

Reproducible numerical checks of the pair, orbital and response calculations in the two
working manuscripts in this repository:

- **LG** — `manuscripts/lg/LG_Master_Unified.tex`
- **TAUM** — `manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex`

Each check rebuilds the objects a passage defines as explicit matrices and asserts the
identities that passage displays. A check either passes or fails; nothing is fitted.

## Running

```bash
pip install -r checks/requirements.txt      # NumPy and SciPy
python checks/run_all.py                    # about 45 seconds
python checks/run_all.py --full             # adds the slow convergence check, about 70 seconds
python checks/run_all.py -k orbital         # only checks whose name contains "orbital"
```

The standard run should end with

```
36 passed, 0 failed, 1 skipped in ...s
```

and the full run with `37 passed, 0 failed, 0 skipped`. The exit status is nonzero if any
check fails. The same files also run under `pytest checks` (set `CHECKS_FULL=1` for the
slow check).

Tested with Python 3.11, NumPy 2.4 and SciPy 1.17.
Timings vary with hardware and BLAS settings.

## What the checks establish

Labels are the `\label{...}` names in the manuscript sources. `test_labels.py` fails if a
check cites a label that no longer exists. It detects removed or renamed targets;
changes to the mathematics under an unchanged label still require review.

`test_linalg.py` checks that nearby resolved eigenvalues belong to disjoint spectral
blocks, so their port weights sum to the total spectral mass without double-counting.
It also checks that numerically coincident eigenvalues are merged at the declared
rounding precision.

### Pair model — LG, `sec:finite-pair-realization`

Four fermion modes with the pair-transfer Hamiltonian `H = -j (P_c^* P_d + P_d^* P_c)`,
built on the 16-dimensional Fock space (`pair/model.py`).

| Check (`pair/test_pair_realization.py`) | Establishes | Labels |
| --- | --- | --- |
| `test_spectrum_and_unique_ground` | Two-particle spectrum `{-j, 0, 0, 0, 0, +j}`; `g` is the unique ground | `eq:pair-parent`, `eq:pair-energy-vectors` |
| `test_preparation_channel` | `A g_Ψ = g`, `A e_Ψ = b`, `A^*A = I`, `Ker A^* = C d` | `eq:pair-preparation-map`, `eq:pair-channel`, `prop:pair-ground-defect` |
| `test_readout_on_carrier`, `test_readout_on_full_two_particle_space` | `Z_θ^2 = AA^*` on the carrier; on the full sector `1 - Z_θ^2 = |d><d| + P_W` | `eq:pair-readout-action`, `eq:pair-defect-measurement`, `eq:pair-ambient-zero` |
| `test_zero_outcome_probability_for_ambient_states` | `Pr(Z_θ = 0)` = line weight + weight in `W`, for arbitrary states | `eq:pair-ambient-zero` |
| `test_port_weights_recover_the_carrier` | The three spectral weights of the preparation port; rank count gives a one-dimensional defect | `eq:pair-port-weights` |
| `test_defect_weight_and_two_line_correlation` | Defect weight `sin^2(4θ) sin^2(jt/2)`; two-frequency ground correlation | `eq:pair-defect-evolution`, `eq:pair-two-frequency-response` |
| `test_generated_algebras_are_full_matrix_algebras` | Both generated algebras are full matrix algebras | `prop:pair-ground-defect` |
| `test_mixer_identities_and_rotation` | `T^2 = 2P_c^*P_d`, `Yg = 0`; `exp(-2iθY)` implements the orbital rotation | `eq:pair-mixer-selection`, `eq:pair-orbital-rotation` |
| `test_endpoint_angles_have_no_defect`, `test_doubled_preparation_channel` | No defect at `θ = 0, π/4`; the doubled channel is onto with a one-dimensional kernel | `sec:finite-pair-realization` |

| Check (`pair/test_pair_controls.py`) | Establishes | Labels |
| --- | --- | --- |
| `test_exchange_completion_spectrum_and_response` | Spectrum, correlation and defect weight with symmetric self-channels `s` | `eq:pair-completed-hamiltonian`, `eq:pair-completed-spectrum`, `eq:pair-completed-response` |
| `test_minimal_completion_removes_the_accessible_defect` | At `s = j` the prepared range is invariant and the response defect vanishes | `eq:pair-completed-spectrum` |
| `test_asymmetric_completion_moves_the_ground` | Unequal self-channels move the ground off `g` | `eq:pair-completed-hamiltonian` |
| `test_pair_controls_preserve_orbital_parity` | Pair-only controls cannot reach `m` from the paired plane | `eq:pair-general-mixer` |
| `test_mixer_ground_certificate`, `test_invariant_mixer_needs_equivalent_representations`, `test_unitary_intertwiner_gives_a_ground_preserving_mixer` | The mixer certificate and the representation criterion | `eq:pair-mixer-ground-certificate`, `prop:pair-invariant-mixer` |
| `test_unused_space_is_not_protected_by_number_conservation` | A number-conserving operator couples `W` to `m` | `eq:pair-ambient-zero` |

### Orbital model — LG, `sec:orbital-ground-source`

Four spin-1/2 fermions in a two-dimensional harmonic trap with Gaussian attraction and an
orbital field term (`orbital/model.py`), in units `ħ = M = ω = 1`.

| Check | Establishes | Labels |
| --- | --- | --- |
| `test_orbital_first_order.py::test_projected_interaction_matrix` | The first-order matrix on `s²p²`, from exact Gaussian integrals; gaps `λ d_a` and `2λ d_a` | `eq:orbital-projected-energies` |
| `test_orbital_first_order.py::test_angular_momentum_on_valence_shell`, `::test_quadrupole_compressions` | `L_tot = 2Y` on the valence shell; `Q_x`, `Q_y` compress to `Z`, `X` | `eq:orbital-derived-current`, `eq:orbital-quadrupoles` |
| `test_orbital_parent.py::test_ground_and_excitation_structure` | Beyond first order: unique `L = 0` singlet ground; gaps approach `λ d_a`, `λ d_a`, `2λ d_a`; next level one trap quantum up | `prop:orbital-full-ground` |
| `test_orbital_parent.py::test_bias_shifts_and_window` | Field shifts are exact; ground unchanged for `|χ| < 1`, replaced at `χ = 1.2` | `prop:orbital-full-ground`, `eq:orbital-minimal-coupling` |
| `test_orbital_parent.py::test_filtered_quadrupole_readout` | On the selected carrier the filtered quadrupole is `|G><e| + |e><G|`, so `Z^2 = AA^*` | `eq:orbital-low-carrier`, `eq:orbital-typed-preparation`, `eq:orbital-filtered-readout` |
| `test_orbital_parent.py::test_low_carrier_dynamics` | Defect weight `sin^2(2Ω_B t)` and the low-carrier correlation | `eq:orbital-low-response` |
| `test_orbital_parent.py::test_unfiltered_quadrupole_response` | `C_B(t) = C_0(t) cos(2Ω_B t)`; only part of the unfiltered state lies in the low carrier | `eq:orbital-unfiltered-response` |
| `test_orbital_parent.py::test_basis_convergence_full` (slow) | The spectral checks repeated with one more oscillator shell | `prop:orbital-full-ground` |

The parent checks diagonalize the Hamiltonian in a truncated oscillator basis (shells
`n_x + n_y ≤ 3`, dimension 2025; `≤ 4`, dimension 11025 with `--full`). The truncation
keeps rotation, spin and time-reversal symmetry exactly. These runs are numerical evidence
at specific parameters (`ξ = 1`, `λ ≤ 0.2`); the proof of the proposition is the argument
in the manuscript, not this computation.

### Port response — TAUM

| Check (`response/test_response.py`) | Establishes | Labels |
| --- | --- | --- |
| `test_finite_response_dimension_count`, `test_relations_among_ranges_without_pairwise_overlap` | `dim H_min = Σ rank R_j`; the defect count; the ambient kernel decomposition | `cor:finite-response-defect`, `eq:finite-response-defect-count`, `eq:ambient-response-defect` |
| `test_permanent_invisibility_is_the_complement_of_the_response_space` | Vectors invisible to the port at all times are exactly `H_min^⊥` | `prop:permanent-port-invisibility` |
| `test_symmetry_resolved_response` | The factorized response, departure and fixed-line weights; overlapping branches | `prop:opposite-sector-response` and its equations |
| `test_three_state_example` | Weights, minimal space, defect line and transfer formula of the example | `ex:split-excitation-response` and its equations |
| `test_intertwining_criterion_and_detuning_cap` | Prepared states never leave the range exactly when `(1 - AA^*) H A = 0`; detuning caps the reachable defect weight | not yet tied to a label |

## Scope

These are mathematical realizations. The checks show that the displayed identities hold
for the stated finite or truncated models, on positive Hilbert spaces. They do not

- derive those models from microscopic generalized trace dynamics;
- determine what fixes the preparation interface or which operations are available;
- depend on, or support, any interpretation of the output carrier.

The pair Hamiltonian is the canonical component reported in T. P. Singh, *Interacting
Generalised Trace Dynamics* (2026), §3.4, Eq. (3.11),
[doi:10.5281/zenodo.23853851](https://doi.org/10.5281/zenodo.23853851). The preparation,
readout, self-channel analysis and the orbital model are constructions of the manuscripts
in this repository.

## Layout

```
checks/
  run_all.py            runner (no dependencies beyond NumPy and SciPy)
  test_labels.py        cited labels must exist in the manuscripts
  test_linalg.py        nearby eigenvalues must not double-count spectral mass
  lib/                  Fock-space operators and small linear-algebra helpers
  pair/                 pair model and its checks
  orbital/              orbital model and its checks
  response/             finite-dimensional response checks
```
