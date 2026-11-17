"""Four-fermion shell reduction in oscillator units (hbar = M = omega = 1)."""
from collections import defaultdict
from itertools import combinations
from math import pi

import numpy as np
from scipy.sparse.linalg import expm_multiply

from orbital.model import build, gaussian_integrals_1d


LABELS = {
    "UM": ["eq:orbital-shell-accuracy", "eq:orbital-second-order"],
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
