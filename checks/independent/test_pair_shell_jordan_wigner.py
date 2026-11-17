"""Pair model and orbital shell rebuilt from a Jordan-Wigner representation (hbar = j = 1).

This file shares no code with pair/model.py or orbital/model.py.  It constructs four
fermion modes directly, restricts to the two-particle space, and recomputes the
statements of Section 8 of the consolidated manuscript that concern the six-dimensional
space: the reduction onto the prepared plane, the two spectral lines and their weights,
the three silent directions, the angular momentum of the shell, and the two detector
readings after a rotation of the preparation.
"""
import numpy as np
from scipy.linalg import expm, null_space

LABELS = {
    "UM": ["eq:prepared-schur", "eq:pair-invariance-amplitude", "subsec:orbital-apparatus",
           "eq:rotated-orbital-readout"],
    "LG": ["eq:pair-parent", "eq:pair-energy-vectors", "eq:pair-preparation-map",
           "eq:orbital-derived-current"],
}

_I2, _SZ = np.eye(2), np.diag([1.0, -1.0])
_LOWER = np.array([[0.0, 1.0], [0.0, 0.0]])


def _annihilator(k, n=4):
    out = np.array([[1.0]])
    for i in range(n):
        out = np.kron(out, _SZ if i < k else _LOWER if i == k else _I2)
    return out


def _dag(x):
    return x.conj().T


# Modes: c = first orbital (p_x in the trap), d = second orbital (p_y); index 1, 2 = the two internal (spin) labels.
C1, C2, D1, D2 = (_annihilator(k) for k in range(4))
_NUMBER = sum(_dag(a) @ a for a in (C1, C2, D1, D2))
P6 = np.eye(16)[:, np.isclose(np.diag(_NUMBER), 2)]          # isometry onto the two-particle space
VAC = np.zeros(16)
VAC[0] = 1.0
PC, PD = C2 @ C1, D2 @ D1


def _on6(operator):
    return P6.T @ operator @ P6


def _to6(vector):
    return P6.T @ vector


def hamiltonian(s, j=1.0):
    return _on6(-s * (_dag(PC) @ PC + _dag(PD) @ PD) - j * (_dag(PC) @ PD + _dag(PD) @ PC))


_CC, _DD = _dag(PC) @ VAC, _dag(PD) @ VAC
G = _to6((_CC + _DD) / np.sqrt(2))
E = _to6((_CC - _DD) / np.sqrt(2))
M = _to6((_dag(C1) @ _dag(D2) + _dag(D1) @ _dag(C2)) @ VAC / np.sqrt(2))
W = null_space(np.column_stack([G, E, M]).T)
# Orbital angular momentum, total spin and its z component on the shell (hbar = 1).
L = _on6(1j * sum(_dag(dk) @ ck - _dag(ck) @ dk for ck, dk in ((C1, D1), (C2, D2))))
_SPLUS = _on6(_dag(C1) @ C2 + _dag(D1) @ D2)
SZ = _on6(0.5 * (_dag(C1) @ C1 - _dag(C2) @ C2 + _dag(D1) @ D1 - _dag(D2) @ D2))
S2 = _dag(_SPLUS) @ _SPLUS + SZ @ SZ + SZ


def preparation(theta):
    """Columns A g_Psi and A e_Psi for the rotated pair preparation."""
    f = [np.cos(theta) * _dag(a) + np.sin(theta) * _dag(b) for a, b in ((C1, D1), (C2, D2))]
    h = [-np.sin(theta) * _dag(a) + np.cos(theta) * _dag(b) for a, b in ((C1, D1), (C2, D2))]
    first, second = _to6(f[0] @ f[1] @ VAC), _to6(h[0] @ h[1] @ VAC)
    return np.column_stack([(first + second) / np.sqrt(2), (first - second) / np.sqrt(2)])


def _line_weights(h, vector):
    values, vectors = np.linalg.eigh(h)
    weights = {}
    for value, column in zip(values, vectors.T):
        weight = abs(np.vdot(column, vector)) ** 2
        if weight > 1e-12:
            key = round(float(value), 9)
            weights[key] = weights.get(key, 0.0) + weight
    return weights


def test_spectrum_and_rank_one_minimal_completion():
    for s, expected in ((0.0, (-1.0, 1.0, 0.0)), (0.4, (-1.4, 0.6, 0.0)), (1.0, (-2.0, 0.0, 0.0))):
        h = hamiltonian(s)
        for vector, energy in zip((G, E, M), expected):
            np.testing.assert_allclose(h @ vector, energy * vector, atol=1e-13)
        np.testing.assert_allclose(W.T @ h @ W, 0, atol=1e-13)
    np.testing.assert_allclose(hamiltonian(1.0), -2 * np.outer(G, G), atol=1e-13)


def test_reduction_onto_the_prepared_plane():
    for theta, s in ((np.pi / 8, 0.0), (np.pi / 18, 0.0), (np.pi / 8, 0.4), (0.7, 0.25), (np.pi / 8, 1.0)):
        a = preparation(theta)
        h = hamiltonian(s)
        b = a[:, 1]
        u, v = np.cos(2 * theta), np.sin(2 * theta)
        d = -v * E + u * M
        np.testing.assert_allclose(b, u * E + v * M, atol=1e-13)
        delta = 1.0 - s
        coupling = (np.eye(6) - a @ a.T) @ h @ a @ a.T
        np.testing.assert_allclose(coupling, -delta * u * v * np.outer(d, b), atol=1e-13)
        complement = null_space(a.T)
        h_q = complement.T @ h @ complement
        for z in (0.37 + 0.2j, -1.4 + 0.05j, 2.2 - 0.3j):
            sigma = a.T @ h @ complement @ np.linalg.solve(z * np.eye(4) - h_q, complement.T @ h @ a)
            measured = a.T @ np.linalg.solve(z * np.eye(6) - h, a)
            np.testing.assert_allclose(measured, np.linalg.inv(z * np.eye(2) - a.T @ h @ a - sigma), atol=1e-11)
            # (H - z)^{-1} convention of the manuscript: weights u^2 at delta and v^2 at zero
            np.testing.assert_allclose(-measured[1, 1], u * u / (delta - z) + v * v / (-z), atol=1e-11)
        weights = _line_weights(h, b)
        if abs(delta) > 1e-12:
            expected = {round(delta, 9): u * u, 0.0: v * v}
            expected = {k: w for k, w in expected.items() if w > 1e-12}
        else:
            expected = {0.0: 1.0}
        assert set(weights) == set(expected)
        for key in expected:
            assert np.isclose(weights[key], expected[key], atol=1e-12)


def test_silent_directions_and_peak_leakage():
    rng = np.random.default_rng(20261117)
    a = preparation(np.pi / 8)
    h = hamiltonian(0.0)
    block = rng.normal(size=(3, 3))
    moved = h + W @ (block + block.T) @ W.T
    assert not np.allclose(np.linalg.eigvalsh(h), np.linalg.eigvalsh(moved))
    for t in np.linspace(0.3, 30, 25):
        np.testing.assert_allclose(a.T @ expm(-1j * t * moved) @ a, a.T @ expm(-1j * t * h) @ a, atol=1e-12)
    for theta in np.deg2rad([5.0, 10.0, 22.5, 40.0]):
        a = preparation(theta)
        b = a[:, 1]
        outside = np.eye(6) - a @ a.T
        for s in (0.0, 0.4):
            h = hamiltonian(s)
            peak_time = np.pi / (1.0 - s)
            leak = np.linalg.norm(outside @ expm(-1j * peak_time * h) @ b) ** 2
            assert np.isclose(leak, np.sin(4 * theta) ** 2, atol=1e-12)
            for t in np.linspace(0.05, 2 * peak_time, 60):
                value = np.linalg.norm(outside @ expm(-1j * t * h) @ b) ** 2
                assert value <= np.sin(4 * theta) ** 2 + 1e-12
                assert np.isclose(value, np.sin(4 * theta) ** 2 * np.sin((1.0 - s) * t / 2) ** 2, atol=1e-12)


def test_shell_angular_momentum_and_symmetry_labels():
    np.testing.assert_allclose(L @ G, 0, atol=1e-13)
    np.testing.assert_allclose(L @ E, 2j * M, atol=1e-13)
    np.testing.assert_allclose(L @ M, -2j * E, atol=1e-13)
    np.testing.assert_allclose(L @ W, 0, atol=1e-13)
    for theta in (0.2, 0.7, np.pi / 8):
        np.testing.assert_allclose(expm(-1j * theta * L) @ np.column_stack([G, E]), preparation(theta), atol=1e-12)
    for first, second in ((L, S2), (L, SZ), (S2, SZ)):
        np.testing.assert_allclose(first @ second, second @ first, atol=1e-13)
    _, vectors = np.linalg.eigh(L + 10 * S2 + 100 * SZ)
    labels = {tuple(int(round(np.vdot(v, op @ v).real)) for op in (L, S2, SZ)) for v in vectors.T}
    assert labels == {(0, 0, 0), (2, 0, 0), (-2, 0, 0), (0, 2, 1), (0, 2, 0), (0, 2, -1)}


def test_rotation_keeps_orbital_weights_and_changes_pair_weights():
    omega_b = 0.13
    orbital = 1.0 * (np.outer(E, E) + np.outer(M, M)) + 2.0 * W @ W.T - omega_b * L
    pair = hamiltonian(0.0)
    np.testing.assert_allclose(orbital @ L, L @ orbital, atol=1e-13)
    assert np.linalg.norm(pair @ L - L @ pair) > 1.0
    for theta in np.deg2rad([0.0, 10.0, 22.5, 40.0]):
        b = preparation(theta)[:, 1]
        orbital_weights = _line_weights(orbital, b)
        assert set(orbital_weights) == {round(1 - 2 * omega_b, 9), round(1 + 2 * omega_b, 9)}
        np.testing.assert_allclose(list(orbital_weights.values()), [0.5, 0.5], atol=1e-12)
        pair_weights = _line_weights(pair, b)
        assert np.isclose(pair_weights.get(1.0, 0.0), np.cos(2 * theta) ** 2, atol=1e-12)
        assert np.isclose(pair_weights.get(0.0, 0.0), np.sin(2 * theta) ** 2, atol=1e-12)


def test_fixed_and_rotating_detectors_on_the_shell():
    omega_b = 0.13
    orbital = 1.0 * (np.outer(E, E) + np.outer(M, M)) + 2.0 * W @ W.T - omega_b * L
    for theta in np.deg2rad([0.0, 10.0, 22.5, 40.0]):
        prepared = preparation(theta)[:, 1]
        for beta in (theta, 0.0, 0.31):
            plane = preparation(beta)
            for t in np.linspace(0, 30, 61):
                evolved = expm(-1j * t * orbital) @ prepared
                amplitude = np.vdot(plane[:, 1], evolved)
                assert np.isclose(amplitude, np.exp(-1j * t) * np.cos(2 * omega_b * t - 2 * theta + 2 * beta), atol=1e-11)
                zero = 1 - np.linalg.norm(plane.conj().T @ evolved) ** 2
                assert np.isclose(zero, np.sin(2 * omega_b * t - 2 * theta + 2 * beta) ** 2, atol=1e-11)
