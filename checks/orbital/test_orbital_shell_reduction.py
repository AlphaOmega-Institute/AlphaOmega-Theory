"""Four-fermion shell reduction in oscillator units (hbar = M = omega = 1)."""
from collections import defaultdict
from itertools import combinations
from math import pi

import numpy as np
from scipy.sparse.linalg import expm_multiply

from orbital.model import build, gaussian_integrals_1d


LABELS = {
    "UM": ["subsec:orbital-apparatus", "eq:orbital-shell-accuracy", "eq:orbital-second-order"],
    "LG": ["eq:orbital-parent", "eq:orbital-minimal-coupling"],
}


def _shell_columns(cutoff, xi):
    """Apply the attraction to six shell determinants without forming the full Hamiltonian.

    A spin orbital is indexed by its spatial index plus spin times the count of orbitals.
    A row is a sorted tuple of four such indices. Each pair interaction replaces two
    entries, retaining their spins and the two spectator orbitals.
    """
    orbitals = [(x, n - x) for n in range(cutoff + 1) for x in range(n + 1)]
    n_orb = len(orbitals)
    x, y = np.array(orbitals).T
    valence = (1, 2, n_orb + 1, n_orb + 2)
    shell = [tuple(sorted((0, n_orb, *pair))) for pair in combinations(valence, 2)]
    integrals = gaussian_integrals_1d(cutoff, xi, screen_parity=False)
    rows = defaultdict(lambda: np.zeros(len(shell)))

    for column, ket in enumerate(shell):
        for i, j in combinations(range(4), 2):
            c, d = ket[i], ket[j]
            cx, cy = orbitals[c % n_orb]
            dx, dy = orbitals[d % n_orb]
            gaussian = integrals[x[:, None], x[None, :], cx, dx] * integrals[y[:, None], y[None, :], cy, dy]
            matrix = -gaussian / (2 * pi * xi)
            spectators = [ket[k] for k in range(4) if k != i and k != j]
            for a, b in zip(*np.nonzero(matrix)):
                a_spin = int(a) + (c // n_orb) * n_orb
                b_spin = int(b) + (d // n_orb) * n_orb
                if a_spin == b_spin or a_spin in spectators or b_spin in spectators:
                    continue
                # Remove c then d, then insert b and a with fermionic signs.
                swaps = i + j - 1 + sum(s < b_spin for s in spectators)
                swaps += sum(s < a_spin for s in spectators) + int(b_spin < a_spin)
                bra = tuple(sorted((*spectators, a_spin, b_spin)))
                rows[bra][column] += (-1)**swaps * matrix[a, b]

    orbital_quanta = x + y
    quanta = orbital_quanta[np.array(list(rows)) % n_orb]
    excitations = quanta.sum(axis=1) - 2
    _, shell_basis = np.linalg.eigh(np.array([rows[row] for row in shell]))
    return excitations, quanta.max(axis=1), np.array(list(rows.values())), shell_basis


def test_parity_and_second_order_orbital_energies():
    excitations, highest, coupling, basis = _shell_columns(10, xi=2.0)
    assert np.max(np.abs(coupling[excitations % 2 == 1])) < 1e-14
    # Project the same quadrature onto the single-particle cutoffs used in the table.
    for cutoff, (ground, chiral, triplet) in [
        (6, (-0.00162813, -0.00170517, -0.00177006)),
        (10, (-0.00162846, -0.00170561, -0.00177062)),
    ]:
        included = (highest <= cutoff) & (excitations > 0)
        amplitudes = coupling[included] @ basis
        correction = -(amplitudes.T / excitations[included]) @ amplitudes
        expected = np.diag([ground, chiral, chiral, triplet, triplet, triplet])
        np.testing.assert_allclose(correction, expected, atol=6e-9, rtol=0)

    # These strengths sum squared amplitudes over all six shell vectors.
    # Pauli exclusion leaves at least one oscillator quantum in the three spectators
    # of a highest occupied orbital: a row k quanta above the shell needs cutoff ≤ k+1.
    for step, expected in [(2, 0.01972668), (4, 0.00179189), (6, 0.00021752)]:
        level = excitations == step
        assert np.max(highest[level]) == step + 1
        assert np.isclose(np.linalg.norm(coupling[level])**2, expected, atol=5e-9, rtol=0)
    six_up_to_six = (excitations == 6) & (highest <= 6)
    assert np.isclose(np.linalg.norm(coupling[six_up_to_six])**2, 0.00020630, atol=5e-9, rtol=0)


def test_compensated_field_preserves_bare_shell_leakage():
    ops = build(2, xi=2.0)
    shell = np.flatnonzero(np.isclose(ops["Hho"].diagonal(), 6.0))
    initial = np.eye(ops["dim"])[:, shell]
    outside = np.ones(ops["dim"], dtype=bool)
    outside[shell] = False
    zero_field = ops["Hho"] - 0.4 * ops["W"]
    biased = zero_field - 0.17 * ops["L"]
    for time in (0.3, 1.1):
        evolved = expm_multiply(-1j * time * zero_field, initial)
        rotated = expm_multiply(-1j * time * biased, initial)
        assert np.linalg.norm(evolved - rotated) > 0.01
        # The Gram matrices give leakage for every initial shell superposition.
        gram = evolved[outside].conj().T @ evolved[outside]
        assert np.linalg.norm(gram) > 1e-6
        np.testing.assert_allclose(rotated[outside].conj().T @ rotated[outside], gram, atol=1e-13, rtol=0)


def test_three_shell_error_bounds_in_truncated_trap():
    """Compare the shell predictions with the full cutoff-two trap in S_z=0.

    Every spin multiplet has a representative here. Energies are evaluated at
    zero field, and the probabilities and amplitude errors use all shell inputs.
    """
    xi, lam = 2.0, 0.1034
    width = 3 * lam / (pi * xi)
    ops = build(2, xi)
    shell = np.flatnonzero(np.isclose(ops["Hho"].diagonal(), 6))
    outside = np.ones(ops["dim"], dtype=bool)
    outside[shell] = False
    h = (ops["Hho"] - lam * ops["W"]).toarray()
    ritz, ritz_vectors = np.linalg.eigh(h[np.ix_(shell, shell)])
    energies, vectors = np.linalg.eigh(h)
    errors = ritz - energies[:len(shell)]
    assert np.all(errors > 1e-8)
    assert np.all(errors <= width**2 / 8)

    # Compare the two propagators on the entire bare shell, including a long
    # wait, by evolving each shell basis vector with the full spectral resolution.
    for t in (0.4, pi / 2, 2.0, 2 * pi, 52 * pi):
        evolved = vectors @ (np.exp(-1j * energies * t)[:, None] * vectors[shell].conj().T)
        shell_evolution = (ritz_vectors * np.exp(-1j * ritz * t)) @ ritz_vectors.conj().T
        leakage = np.linalg.norm(evolved[outside], ord=2)**2
        drift = np.linalg.norm(evolved[shell] - shell_evolution, ord=2)
        assert leakage <= (width / 2)**2
        assert drift <= width**2 * t / 4
        if t == pi / 2:
            assert leakage > 1e-8 and drift > 1e-8


def test_sector_isolation_and_seventh_level_in_truncated_trap():
    """Resolve angular momentum and spin and count levels with their spin multiplicities."""
    xi, lam = 2.0, 0.1034
    width = 3 * lam / (pi * xi)
    omega_b = 0.04
    ops = build(2, xi)
    h = (ops["Hho"] - lam * ops["W"]).toarray()
    angular, angular_vectors = np.linalg.eigh(ops["L"].toarray())
    angular = np.rint(angular).astype(int)
    shell_sectors = {(0, 0), (0, 1), (-2, 0), (2, 0)}
    verified = set()
    biased_spectrum = []
    for m in np.unique(angular):
        space = angular_vectors[:, angular == m]
        spin_squared, spin_vectors = np.linalg.eigh(space.conj().T @ (ops["S2"] @ space))
        spin_squared = np.rint(spin_squared).astype(int)
        for s in (0, 1, 2):
            columns = spin_squared == s * (s + 1)
            if not columns.any():
                continue
            basis = space @ spin_vectors[:, columns]
            levels = np.linalg.eigvalsh(basis.conj().T @ h @ basis)
            if (m, s) in shell_sectors:
                assert levels[0] <= 6
                assert levels[1] >= 8 - width
                assert levels[1] - levels[0] >= 2 - width
                verified.add((m, s))
            biased_spectrum.extend(np.repeat(levels - omega_b * m, 2 * s + 1))
    assert verified == shell_sectors
    sorted_levels = np.sort(biased_spectrum)
    assert len(sorted_levels) == 495  # four fermions in twelve spin orbitals
    assert sorted_levels[5] < 7 - 3 * abs(omega_b) - width <= sorted_levels[6]
