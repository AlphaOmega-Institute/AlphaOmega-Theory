"""Clauses of LG's canonical deconditioning theorem that the consolidated manuscript does
not yet state (ledger rows: strong continuity, the trace-norm estimate with constant 2,
and the generator), together with the three it leaves one step away.

    D_t(rho) = P rho P + r^(1/2) (P rho Q + Q rho P) + r Q rho Q + (1 - r) Tr(Q rho) P,
    r = exp(-gamma t),  P the ground projection,  Q = I - P.
"""
import numpy as np

LABELS = {
    "LG": ["thm:canonical-deconditioning-semigroup", "eq:canonical-deconditioning-channel",
           "eq:deconditioning-exact-contraction", "eq:deconditioning-trace-distance-bound",
           "eq:deconditioning-generator"],
}


def _channel(rho, gamma, t):
    n = rho.shape[0]
    p = np.zeros((n, n))
    p[0, 0] = 1
    q = np.eye(n) - p
    r = np.exp(-gamma * t)
    return p @ rho @ p + np.sqrt(r) * (p @ rho @ q + q @ rho @ p) + r * q @ rho @ q + (1 - r) * np.trace(q @ rho).real * p


def _state(rng, n):
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    rho = raw @ raw.conj().T
    return rho / np.trace(rho).real


def test_semigroup_fixed_ground_and_exact_decay():
    rng = np.random.default_rng(12)
    for _ in range(300):
        n = rng.integers(2, 7)
        gamma, s, t = rng.uniform(0.2, 2), rng.uniform(0, 3), rng.uniform(0, 3)
        rho = _state(rng, n)
        np.testing.assert_allclose(_channel(_channel(rho, gamma, s), gamma, t), _channel(rho, gamma, s + t), atol=1e-12)
        out = _channel(rho, gamma, t)
        assert np.isclose(np.trace(out).real, 1.0) and np.linalg.eigvalsh(out).min() > -1e-12
        ground = np.zeros((n, n))
        ground[0, 0] = 1
        np.testing.assert_allclose(_channel(ground, gamma, t), ground, atol=1e-13)
        deviation = 1 - rho[0, 0].real
        assert np.isclose(1 - out[0, 0].real, np.exp(-gamma * t) * deviation)


def test_trace_norm_estimate_with_constant_two():
    rng = np.random.default_rng(13)
    closest = 0.0
    for _ in range(2000):
        n = rng.integers(2, 7)
        gamma, t = rng.uniform(0.2, 2), rng.uniform(0, 6)
        rho = _state(rng, n)
        ground = np.zeros((n, n))
        ground[0, 0] = 1
        distance = np.abs(np.linalg.eigvalsh(_channel(rho, gamma, t) - ground)).sum()
        bound = 2 * np.exp(-gamma * t / 2) * np.sqrt(1 - rho[0, 0].real)
        assert distance <= bound + 1e-12
        closest = max(closest, distance / bound)
    assert closest > 0.9
    # pure states close to the ground approach the constant
    for angle in (0.3, 0.05, 0.005):
        vector = np.array([np.cos(angle), np.sin(angle)])
        rho = np.outer(vector, vector)
        distance = np.abs(np.linalg.eigvalsh(rho - np.diag([1.0, 0.0]))).sum()
        assert np.isclose(distance, 2 * np.sin(angle))            # equals the bound at t = 0


def test_generator_and_unique_fixed_state():
    rng = np.random.default_rng(14)
    for _ in range(200):
        n = rng.integers(2, 6)
        gamma = rng.uniform(0.2, 2)
        rho = _state(rng, n)
        p = np.zeros((n, n))
        p[0, 0] = 1
        q = np.eye(n) - p
        generator = gamma * (np.trace(q @ rho).real * p - 0.5 * (q @ rho + rho @ q))
        h = 1e-6
        np.testing.assert_allclose((_channel(rho, gamma, h) - _channel(rho, gamma, -h)) / (2 * h), generator, atol=1e-8)
        # continuity at t = 0 in trace norm
        assert np.abs(np.linalg.eigvalsh(_channel(rho, gamma, 1e-9) - rho)).sum() < 1e-7
        # as a linear map on matrices D_t has eigenvalue one exactly once, at any positive time
        t = rng.uniform(0.1, 3)
        basis = np.eye(n * n).reshape(n * n, n, n)
        matrix = np.array([_channel(b.astype(complex), gamma, t).reshape(-1) for b in basis]).T
        assert np.sum(np.isclose(np.linalg.eigvals(matrix), 1.0, atol=1e-9)) == 1
