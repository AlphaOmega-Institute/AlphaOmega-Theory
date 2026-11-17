"""Signed line weights and the calibration of the two ports (hbar = 1).

Planned for the section on signed response (outline III.3; working title Section 9).
Nothing here is in the consolidated manuscript yet, so the file cites no label.

Statements checked:
  1. A scalar response c(t) = sum_E w_E exp(-i E t) with real weights can be written as
     <v, U(t) u> for an ordinary unitary group exactly when the two port vectors satisfy
     ||u|| ||v|| >= sum_E |w_E|, and the value sum_E |w_E| is attained.
  2. The peak |c(t)| is only a lower bound for that minimum: weights (1, 1, -1) at
     energies (0, 1, 2) have peak sqrt(5) and minimum 3.
  3. The two-level generator with levels 0 and 3 and return coupling of the opposite sign
     has, at coupling sqrt(2), energies 1 and 2 with weights 2 and -1 at the first level;
     c(t) = 2 exp(-i t) - exp(-2 i t), so |c|^2 = 5 - 4 cos t.  At coupling 3/2 the two
     energies meet, and at coupling 2 they are complex.
"""
import numpy as np
from scipy.linalg import expm


def _response(weights, energies, times):
    return (np.asarray(weights)[:, None] * np.exp(-1j * np.outer(energies, times))).sum(axis=0)


def test_minimal_port_norm_product_is_the_sum_of_absolute_weights():
    times = np.linspace(0, 2 * np.pi, 20001)
    for weights, energies, peak in (((2.0, -1.0), (1.0, 2.0), 3.0),
                                    ((1.0, 1.0, -1.0), (0.0, 1.0, 2.0), np.sqrt(5.0)),
                                    ((1.5, -1.0, 0.5), (0.0, 1.0, 2.0), 3.0)):
        weights, energies = np.array(weights), np.array(energies)
        assert np.isclose(np.abs(_response(weights, energies, times)).max(), peak, atol=1e-6)
        u = np.sqrt(np.abs(weights))
        v = np.sign(weights) * u
        assert np.isclose(np.linalg.norm(u) * np.linalg.norm(v), np.abs(weights).sum())
        for t in (0.0, 0.4, 2.2, np.pi):
            assert np.isclose(np.vdot(v, np.exp(-1j * energies * t) * u), _response(weights, energies, [t])[0])
        assert peak <= np.abs(weights).sum() + 1e-12


def test_every_unitary_realization_obeys_the_lower_bound():
    rng = np.random.default_rng(31)
    for _ in range(400):
        n = rng.integers(2, 8)
        raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        values, vectors = np.linalg.eigh(raw + raw.conj().T)
        u = rng.normal(size=n) + 1j * rng.normal(size=n)
        v = rng.normal(size=n) + 1j * rng.normal(size=n)
        weights = (vectors.conj().T @ v).conj() * (vectors.conj().T @ u)      # <v, P_E u> for each simple level
        assert np.abs(weights).sum() <= np.linalg.norm(u) * np.linalg.norm(v) + 1e-12
        for t in (0.3, 1.7):
            unitary = (vectors * np.exp(-1j * values * t)) @ vectors.conj().T
            assert np.isclose(np.vdot(v, unitary @ u), (weights * np.exp(-1j * values * t)).sum())
            assert abs(np.vdot(v, unitary @ u)) <= np.linalg.norm(u) * np.linalg.norm(v) + 1e-12
        # equal ports give nonnegative weights
        assert np.all(((vectors.conj().T @ u).conj() * (vectors.conj().T @ u)).real >= -1e-14)


def test_two_level_generator_with_the_opposite_return_sign():
    port = np.array([1.0, 0.0])
    coupling = np.sqrt(2.0)
    generator = np.array([[0.0, coupling], [-coupling, 3.0]])
    values, vectors = np.linalg.eig(generator)
    order = np.argsort(values.real)
    values, vectors = values[order], vectors[:, order]
    np.testing.assert_allclose(values, [1.0, 2.0], atol=1e-12)
    weights = vectors[0, :] * np.linalg.inv(vectors)[:, 0]
    np.testing.assert_allclose(weights, [2.0, -1.0], atol=1e-12)
    for z in (0.3 + 0.2j, -1.0 + 0.5j, 2.5 - 0.1j):
        measured = port @ np.linalg.solve(generator - z * np.eye(2), port)
        assert np.isclose(measured, (3 - z) / ((1 - z) * (2 - z)))
        assert np.isclose(measured, 2 / (1 - z) - 1 / (2 - z))
        # the same function from the ordinary system diag(1, 2) with two different ports
        prepare, read = np.array([np.sqrt(2.0), 1.0]), np.array([np.sqrt(2.0), -1.0])
        assert np.isclose(measured, sum(read[k] * prepare[k] / (level - z) for k, level in enumerate((1.0, 2.0))))
    for t in np.linspace(0, 4 * np.pi, 81):
        amplitude = (expm(-1j * t * generator))[0, 0]
        assert np.isclose(amplitude, 2 * np.exp(-1j * t) - np.exp(-2j * t), atol=1e-10)
        assert np.isclose(abs(amplitude) ** 2, 5 - 4 * np.cos(t), atol=1e-10)
    assert np.isclose(expm(-1j * np.pi * generator)[0, 0], -3.0, atol=1e-10)
    # exceptional point and the unstable side
    np.testing.assert_allclose(np.linalg.eigvals(np.array([[0.0, 1.5], [-1.5, 3.0]])), [1.5, 1.5], atol=1e-6)
    unstable = np.linalg.eigvals(np.array([[0.0, 2.0], [-2.0, 3.0]]))
    np.testing.assert_allclose(sorted(unstable.imag), [-np.sqrt(7) / 2, np.sqrt(7) / 2], atol=1e-12)
    # with the ordinary sign the same coupling gives positive weights
    values, vectors = np.linalg.eigh(np.array([[0.0, coupling], [coupling, 3.0]]))
    assert np.all(vectors[0, :] ** 2 > 0) and np.isclose((vectors[0, :] ** 2).sum(), 1.0)
