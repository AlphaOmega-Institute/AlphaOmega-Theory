# Verification checks

Reproducible numerical checks of the mathematical results and worked examples in the
working manuscripts and the consolidated draft:

- **LG** — `manuscripts/lg/LG_Master_Unified.tex`
- **TAUM** — `manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex`
- **UM** — `manuscripts/consolidated/Universal_Model.tex`

Each check rebuilds the objects a passage defines as explicit matrices and asserts the
identities that passage displays. A check either passes or fails; nothing is fitted.

## Running

```bash
pip install -r checks/requirements.txt      # NumPy and SciPy
python checks/run_all.py                    # standard run
python checks/run_all.py --full             # includes the orbital convergence and higher-cutoff checks
python checks/run_all.py -k orbital         # only checks whose name contains "orbital"
```

The standard run should end with

```
194 passed, 0 failed, 2 skipped in ...s
```

and the full run with `196 passed, 0 failed, 0 skipped`. The exit status is nonzero if any
check fails. The same files also run under `pytest checks` (set `CHECKS_FULL=1` for the
slow checks).

Tested with Python 3.12 using NumPy 2.3.5 and SciPy 1.17.0, and using NumPy 2.0.2 and SciPy 1.14.1.
Timings vary with hardware and BLAS settings.

## What the checks establish

Labels are the `\label{...}` names in the manuscript sources. `test_labels.py` fails if a
check cites a label that no longer exists. It detects removed or renamed targets;
changes to the mathematics under an unchanged label still require review.

`test_linalg.py` checks that nearby resolved eigenvalues belong to disjoint spectral
blocks, so their port weights sum to the total spectral mass without double-counting.
It also checks that numerically coincident eigenvalues are merged at the declared
absolute spectral tolerance. A regression covers square roots and inverse square roots
of positive Hermitian matrices with real and complex floating-point inputs.

### Signed orientations and operator domains

Six checks recover the channel of an odd Krein-self-adjoint matrix, verify both ordered
return signs, and factor the two branches for a complex channel with a kernel and a
defect. Two diagonal examples have exact cutoff sums. For one, the first step has
finite norm while the norms of its truncated returns diverge. For the other, the generalized
inverse is unbounded while the polar phase extends to the whole active support.
A proper isometric embedding checks the difference between norm preservation on
the input sector and nilpotency on the direct sum.
[Signed branch checks](test_odd_branches.py)

### Signed Schur reduction

Seven checks cover the factorization, its four inverse blocks, and the imaginary
and below-threshold signs of the return. A diagonal unbounded example uses
exact cutoff sums to test which triangular factor preserves the operator domain.
In the balanced pair, reversing the lower coupling moves the real poles off
the real axis. A three-level example checks how a scalar port can cancel a
negative spectral weight while the full retained resolvent still detects it.
The spectral-limit check compares the resolvent near an eigenvalue with its
growth for a continuous uniform spectral measure.
[Signed Schur checks](test_signed_reduction.py)

### Ground correspondence

| Check (`test_ground_bridge.py`) | Establishes | Labels |
| --- | --- | --- |
| `test_raw_ground_bridge_and_polar_ground_measurement` | Two choices of ground for the same weighted channel satisfy the forward bridge. The output-ground measurement under polar transport gives probability 9/10 or 1. The five strong conditions fail or hold together in those cases. | UM `prop:ground-polar-compatibility`, LG `thm:ground-channel-polar`, `thm:strong-ground-equivalences` |
| `test_normalizing_polar_image_can_hide_a_kernel_component` | For a projection, a ground with a kernel component maps onto the designated output ray with squared norm 1/2. Normalizing the image discards that loss. | UM `prop:ground-polar-compatibility` |
| `test_path_ground_measurement_depends_on_cell_sizes` | Path Hamiltonians select constant grounds. Pullback across cells of sizes 1 and 4 gives ground-test probability 9/10 and excess energy 1/10. Equal cell sizes give probability 1 and zero excess energy. | UM `prop:ground-polar-compatibility`, `eq:ground-return-energy` |

Seven checks for §§2–4 cover spectral measurements on the active and full output spaces,
the graph projections and density operators, and the compression of the full polar
operators to the ground plane. They compare the two graph weights with the ground
fidelity, including when the input has a component in the kernel.
[Polar graph checks](test_polar_graphs.py)

Seven checks for §4.5 cover the energy–fidelity factorization, the gap bound and the
state left by a Lüders ground test. They include the number-operator example of a
finite-energy state outside the operator domain and the trap compression that preserves
mean energy while changing its distribution.
[Fidelity and energy checks](test_fidelity_energy.py)

Seven checks across §§4–7 distinguish operators that preserve the ground line from
unitary groups that return to it at isolated times. They compare the shifted input
return on its active space and on the whole input, identify the output return's
zero eigenspace, and verify that channels with the same preparation can move or
preserve the ground although the differential squares to zero.
[Ground preservation and return checks](test_ground_preservation.py)

### Deconditioning across two moments

`test_deconditioning.py` checks the relaxation law on the two coarse sites by comparing its action on all four matrix units with deconditioning on the five fine sites. The target passes the coarse ground test with probability 9/10. At half relaxation an orthogonal pure input becomes maximally mixed. The tests also verify the exponential decay of excess weight for a paired target and the maximum leakage of 1/6 for a target with one third of its weight outside the paired space (UM `prop:two-moment-deconditioning`).

Eight checks for the general relaxation law compare its generator with the channel at finite times,
verify the sharp trace-norm bound on pure inputs, and exercise ground-carrying channels and
conditional outcomes. They also test a one-step unitary dilation, successive steps with varying
ground spaces, and return instruments on the unequal-cell path.
[Relaxation checks](test_relaxation.py)

### Index and the common defect

`test_index.py` checks the two input columns and three output coordinates that contain the entire change to LG's unilateral shift. Rotations of one microradian and π/2 each preserve the ground and the isometry, change the map by rank one, and leave zero common defect with the original shift. A second check extends a map with a one-dimensional kernel and a two-dimensional defect. It verifies how many kernel and defect dimensions remain when the new input columns cover one or both defect directions (UM `eq:index-shift-family`, `eq:index-extension-counts`).

### Determination on the five-site path

`test_determination.py` checks the centered cell and site orbits (dimensions 1 and 4), their covariance Gram matrices, and the cell orbit's equality to the preparation range. Two unit states have identical site probabilities and respective preparation-range probabilities 1 and 0 (UM `prop:cyclic-determination-content`, `eq:path-readout-states`; LG `thm:concrete-GNS-appearance`, `prop:covariance-content-inner-product`).

On the same path, a determination at the second site takes the cyclic ground outside the space generated by the cell preparations (UM `prop:cyclic-common-excess`).

A free wait under the path Hamiltonian changes the site-3 readout. The second derivatives of the two probabilities at zero are +1/8 and −1/8 in units of τ = t/ℏ. Direct evolution at a short time also checks the resulting τ²/8 coefficient in their difference.

Six checks for §5 reconstruct the cyclic representation, compare the generated and
dark subspaces, and test the independence of uniqueness, cyclicity, and a gap.
Purity is tested on the generated space, including the scalar algebra acting on
two or three dimensions.
[Cyclicity and ground-state checks](test_cyclic_ground_data.py)

Four checks for §5.4 compare the algebra and ground orbit generated by the path
Hamiltonian and each single-site projection. They test word order, including
vanishing outcomes and different vectors on the same ray, agreement on a
chosen readout algebra, and the ground character of a family that preserves
the ground line under both its generators and their adjoints.
[Generator family checks](test_generator_families.py)

Five checks for §5.5 verify the grading, the weight outside the sharp frame, and
the spectral construction. They compare ground and differentiated states at zero
frame weight, test the support bound for unsharp effects, and compute the
conditional states of two operations with the same outcome effect.
[Duality frame checks](test_duality_frame.py)

### Pair model — LG, `sec:finite-pair-realization`

Four fermion modes with the pair-transfer Hamiltonian `H = -j (P_c^* P_d + P_d^* P_c)`,
built on the 16-dimensional Fock space (`pair/model.py`).

| Check (`pair/test_pair_realization.py`) | Establishes | Labels |
| --- | --- | --- |
| `test_spectrum_and_unique_ground` | Two-particle spectrum `{-j, 0, 0, 0, 0, +j}`; `g` is the unique ground | `eq:pair-parent`, `eq:pair-energy-vectors` |
| `test_preparation_channel` | `A g_Ψ = g`, `A e_Ψ = b`, `A^*A = I`, `Ker A^* = C d` | `eq:pair-preparation-map`, `eq:pair-channel`, `prop:pair-ground-defect` |
| `test_readout_on_carrier`, `test_readout_on_full_two_particle_space` | `Z_θ^2 = AA^*` on the carrier; on the full sector `1 - Z_θ^2 = |d><d| + P_W` | `eq:pair-readout-action`, `eq:pair-defect-measurement`, `eq:pair-ambient-zero` |
| `test_zero_outcome_probability_for_ambient_states` | `Pr(Z_θ = 0)` = line weight + weight in `W`, for arbitrary states | `eq:pair-ambient-zero` |
| `test_matched_occupation_probabilities_before_and_after_waiting` | The three imbalance projections give the two transported occupations and the defect outcome. At the transfer time the two signed outcomes each pull back to half the ground projection and the zero outcome to the excitation projection. The 36:64 mixture has probabilities (0.25, 0.25, 0.50). A coherent input with the same initial occupations gives (0.49, 0.49, 0.02). | UM `eq:pair-matched-probabilities`, `eq:pair-transfer-povm` |
| `test_port_weights_recover_the_carrier` | The three spectral weights of the preparation port; rank count gives a one-dimensional defect | `eq:pair-port-weights` |
| `test_defect_weight_and_two_line_correlation` | Two-frequency ground correlation for bare transfer. Across symmetric self-channel strengths, including the minimal completion, the input effect of a zero count is the excitation projection times `sin²(4θ) sin²((j−s)t/2)`. | `eq:pair-defect-evolution`, `eq:pair-two-frequency-response`, UM `eq:pair-detector-pullback` |
| `test_generated_algebras_are_full_matrix_algebras` | Both generated algebras are full matrix algebras | `prop:pair-ground-defect` |
| `test_mixer_identities_and_rotation` | `T^2 = 2P_c^*P_d`, `Yg = 0`; `exp(-2iθY)` implements the orbital rotation | `eq:pair-mixer-selection`, `eq:pair-orbital-rotation` |
| `test_endpoint_angles_have_no_defect`, `test_doubled_preparation_channel` | Each endpoint fills its two-dimensional minimal space. The π/8 defect vector projects with weight 1/2 onto each endpoint plane. Together the planes span the three-dimensional carrier. The coherent map has a one-dimensional kernel. | LG `sec:finite-pair-realization`, UM `eq:pair-two-setting-ranges` |

| Check (`pair/test_pair_controls.py`) | Establishes | Labels |
| --- | --- | --- |
| `test_exchange_completion_spectrum_and_response` | Spectrum, correlation and defect weight with symmetric self-channels `s` | `eq:pair-completed-hamiltonian`, `eq:pair-completed-spectrum`, `eq:pair-completed-response` |
| `test_minimal_completion_removes_the_accessible_defect` | At `s = j` the prepared range is invariant and the response defect vanishes | `eq:pair-completed-spectrum` |
| `test_projected_generator_on_ground_and_excitation` | The projected generator vanishes on the ground and gives the signed vector `−(j−s) cos(2θ) sin(2θ) d` on the excitation across five angles and four self-couplings | `eq:finite-preparation-invariance`, `eq:pair-completed-spectrum`, UM `eq:pair-invariance-amplitude` |
| `test_asymmetric_completion_moves_the_ground` | With unequal self-channels the prepared input ground has a nonzero `m` component and fidelity below one with the output ground | `eq:pair-completed-hamiltonian`, UM `subsec:pair-grounds-excitations` |
| `test_pair_controls_preserve_orbital_parity` | Pair-only controls cannot reach `m` from the paired plane | `eq:pair-general-mixer` |
| `test_mixer_ground_certificate`, `test_invariant_mixer_needs_equivalent_representations`, `test_unitary_intertwiner_gives_a_ground_preserving_mixer` | The mixer certificate and the representation criterion | `eq:pair-mixer-ground-certificate`, `prop:pair-invariant-mixer` |
| `test_unused_space_is_not_protected_by_number_conservation` | A number-conserving operator couples `W` to `m` | `eq:pair-ambient-zero` |
| `test_leakage_during_wait_changes_zero_readout` | Starting with zero weight in `W`, evolves the perturbed Hamiltonian and measures the imbalance. At the ideal transfer time its zero-count excess equals the acquired `W` weight. The line, leakage and detector probabilities agree with the exact solution at perturbations of 0, 0.03 and 0.10 times `j`. | `eq:pair-ambient-zero`, `eq:pair-defect-evolution`, UM `eq:pair-leakage-readout` |

The eleven §13 checks verify the orbital currents, the ground-preserving mixers and
the imbalance readout on the three-state carrier and on the full two-particle space.
They also check the signed triplet gap and compare the two-dimensional cyclic space
of the trap excitation with the three-dimensional space reached by the preparation.
[Pair-current checks](pair/test_pair_currents.py)

| Check (`pair/test_phase_graphs.py`) | Establishes | Consolidated labels |
| --- | --- | --- |
| `test_phase_graph_projectors_and_sector_population` | Graph projections, Hilbert normalization, Krein neutrality, the sine and cosine readouts, and paired effects across sector populations for fixed content | `eq:polar-graph-projector`, `prop:phase-readout` |
| `test_phase_readout_depends_on_active_cross_block_and_input` | An ambient cross block can vanish on the active support. A nonzero active cross block can detect phase for one input and give a constant expectation for another. | `prop:phase-readout` |
| `test_pi_over_four_pulse_and_pi_over_two_swap` | The π/4 pulse gives a signed mean of cos φ and an input-sector probability of (1 + cos φ)/2. The π/2 swap reverses the grading and preserves the block-diagonal readouts. | `eq:rotated-sector-imbalance`, `eq:operational-phase-blindness` |
| `test_pair_phase_graph_evolution_and_departure` | The minimal completion preserves the graphs; an input energy shift rotates their phase. At bare complete transfer the evolved excitation has squared projection 1/4 on every original phase graph. | `eq:graph-evolution-phase`, `eq:pair-graph-departure` |

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

[The passive-response checks](response/test_passive_response.py) exercise the Herglotz and Schur kernels, the large-imaginary-energy bound, and the two-level signed witness. They compare spectral weights with time evolutions, recover the unitary between two minimal realizations, and test spectral grouping with nearby levels and with a genuine degeneracy. The symmetry and heat-smoothing checks work on finite models.

### Moving subspaces and curvature

The nine checks for §10 compute projector derivatives, curvature, changes of frame and
parallel transport on finite-dimensional families. A rank-two loop compares the
determinant of its transport with the integrated trace of its curvature. Surface and
four-sphere examples test the characteristic forms. In the pair model the checks recover
the phase from the measured line weight, integrate the flux over the sphere of excited
lines, and verify its cancellation over a full period of the apparatus setting.
[Moving-subspace checks](geometry/test_moving_subspace.py)

### Compact symmetry averaging

The seven checks for §11 compute invariant projections for permutations, circles and
spin representations of SU(2). They compare the cube and full rotation groups on
quartics, and test symmetry-compatible maps and compressed port responses. In the
fermion pair, quarter turns give the full circle's fixed-vector projection at every
coupling. The character-weighted average recovers the three-dimensional space whose
rays are invariant and whose vectors change sign under orbital exchange.
[Symmetry-averaging checks](geometry/test_symmetry_average.py)

### Four maps with one algebra

Eight checks exercise the scalar, translation, oscillator and twisted-interval maps of §12.
They include the three-state transfer and its finite-window bound, the multiplicity in the
oscillator heat trace and Mellin transform, and the response and smoothed mass for opposite twists.
[Four-map checks](response/test_four_maps.py)

### Numerical tolerance for spectral grouping

The response checks group eigenvalues into levels with an absolute tolerance of `1e-9` (`SPECTRAL_TOLERANCE` in `lib/linalg.py`). Groups are disjoint and exhaustive, and no group spans more than the tolerance. Each group starts at its lowest member and contains only energies within that tolerance; nearby neighbours may fall in different groups when joining them would exceed the group's total-span limit. The group's mean represents its energy. Merging genuinely distinct energies approximates the exact response, which needs a separate direction for each distinct energy the port couples to. The reported response-generated space then belongs to this grouped Hamiltonian and need not reduce the original one. The `tol` argument of `minimal_space` controls the separate SVD rank cutoff. `test_linalg.py` covers the boundary cases.

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
  test_ground_bridge.py raw and polar transport of ground states
  test_ground_preservation.py ground line, unitary returns and the two return operators
  test_polar_graphs.py  polar transport, graph states and the ground plane
  test_fidelity_energy.py excitation energy, ground fidelity and the spectral gap
  test_deconditioning.py paired relaxation and leakage on the five-site path
  test_relaxation.py    relaxation semigroup, successive moments and return instruments
  test_odd_branches.py  signed extensions, natural domains and polar factors
  test_signed_reduction.py  Schur reduction, effective signs and port visibility
  test_index.py         rotated shift and dimensions after input enlargement
  test_determination.py cyclic generation and readout on a finite path
  test_cyclic_ground_data.py pointed states, the cyclic representation and ground conditions
  test_generator_families.py generator families, word order and content equivalence
  test_duality_frame.py sharp and unsharp frame effects, grading and conditional states
  lib/                  Fock-space operators and small linear-algebra helpers
  pair/                 pair model, controls, §13 currents and readout
  orbital/              orbital model and its checks
  response/             port responses and worked examples in §12
  geometry/             moving subspaces, characteristic forms and compact symmetry averaging
```
