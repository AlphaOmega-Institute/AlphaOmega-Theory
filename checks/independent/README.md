# Independent checks

These files recompute results already checked elsewhere in `checks/`, by a second route
that shares no model code with the first. Agreement between the two routes is the point.
They were written during review of the consolidated manuscript and then cleaned up so
that they run with the rest of the suite.

```bash
python checks/run_all.py -k independent/         # about two seconds
python checks/run_all.py --full -k independent/  # adds the ten-quanta orbital cutoff
```

| File | Second route | What it adds to the committed checks |
| --- | --- | --- |
| `test_pair_shell_jordan_wigner.py` | Four fermion modes built by Jordan–Wigner in this file. No import from `pair/model.py` or `orbital/model.py`. | The reduction onto the prepared plane for several angles and self-couplings; peak leakage `sin²(4θ)` across the symmetric family; the six symmetry labels of the shell; the detector reading at a general detector angle. |
| `test_shell_bounds_random.py` | The hypotheses of UM `prop:shell-accuracy` sampled directly: a six-fold lowest level, a gap, and an interaction of given spectral width placed anywhere. | The three estimates on 600 random models with widths up to 1.5 gaps (the committed checks use two-level systems and one truncated trap); the energy estimate for widths up to three gaps; the gap of an isolated branch. |
| `test_orbital_four_fermion.py` | Slater determinants of Cartesian oscillator orbitals with the interaction applied to the six shell determinants; Gauss–Hermite quadrature in centre and relative coordinates. | Cross-implementation agreement with `orbital/test_orbital_shell_reduction.py` on the parity rule, the table of second-order coefficients and the level strengths; the absolute first-order ground energy at three values of ξ. |
| `test_covariant_channels_random.py` | Random Kraus operators within each energy-shift class, several per class. | Block-diagonality across the missing line for channels beyond the single maximal-support channel of the committed check; the two-factor coherence bound on 1,500 random channels with no covariance imposed. |

Random numbers use fixed seeds. Each test that samples also asserts that the sample comes
reasonably close to the bound it tests, so a vacuous sample fails.
