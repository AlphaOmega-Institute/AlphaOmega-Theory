"""Relations between an earlier and a later cognition (outline, Part V, planned entry
"Grasper-grasped frame and recollection").  Nothing here is in the consolidated manuscript
yet, so the file cites no label.

Statements checked, for an isometry W from the earlier state space into the later one:
  1. Relation between rays: the component of the later unit vector eta in Ran W is nonzero
     and has the ray of W phi.  It does not depend on the phase chosen for phi.
  2. Matched effects: conditional on that component, W E W* has the probability of E.
  3. Mixing equality: if G' W = W G then
         <eta, G' eta> = q <phi, G phi> + (1 - q) r,   q = ||W W* eta||^2,   0 <= r <= 1.
     An effect that differs from W G W* only outside Ran W still satisfies it.  Without the
     intertwining relation the later attribute can vanish while the presented one is certain.
  4. An aspect of a paired state: for a unit vector phi in D and a projection P with
     w = ||P phi||^2, the part of P phi in D has norm at least w, the part outside has
     squared norm at most w (1 - w) with equality when D is the line through phi, and the
     outside parts of P phi and (I - P) phi cancel.
"""
import numpy as np


def _vector(rng, n):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


def _unit(v):
    return v / np.linalg.norm(v)


def _isometry(rng, m, n):
    q, _ = np.linalg.qr(rng.normal(size=(m, n)) + 1j * rng.normal(size=(m, n)))
    return q


def _effect(rng, n):
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    positive = raw @ raw.conj().T
    return positive / (1.05 * np.linalg.eigvalsh(positive).max())


def test_ray_relation_and_matched_effects():
    rng = np.random.default_rng(99)
    for _ in range(500):
        n = rng.integers(2, 6)
        m = n + rng.integers(1, 6)
        w = _isometry(rng, m, n)
        phi = _unit(_vector(rng, n))
        rest = (np.eye(m) - w @ w.conj().T) @ _vector(rng, m)
        eta = _unit(_vector(rng, 1)[0] * w @ phi + rest)
        branch = w @ w.conj().T @ eta
        q = np.vdot(branch, branch).real
        assert q > 0
        for phase in (1.0, np.exp(1.1j), -1.0):
            assert np.isclose(abs(np.vdot(_unit(branch), w @ (phase * phi))), 1.0)
        effect = _effect(rng, n)
        assert np.isclose(np.vdot(eta, w @ effect @ w.conj().T @ eta).real / q, np.vdot(phi, effect @ phi).real)


def test_mixing_equality_under_intertwining():
    rng = np.random.default_rng(314)
    for _ in range(500):
        n = rng.integers(2, 5)
        extra = rng.integers(2, 5)
        m = n + extra
        w = _isometry(rng, m, n)
        g = _effect(rng, n)
        full, _ = np.linalg.qr(np.hstack([w, _vector(rng, m * extra).reshape(m, extra)]))
        outside = full[:, n:]
        g_later = w @ g @ w.conj().T + outside @ _effect(rng, extra) @ outside.conj().T
        np.testing.assert_allclose(g_later @ w, w @ g, atol=1e-12)
        assert np.linalg.norm(g_later - w @ g @ w.conj().T, 2) > 1e-3
        phi = _unit(_vector(rng, n))
        eta = _unit(_vector(rng, 1)[0] * w @ phi + outside @ _vector(rng, extra))
        q = np.linalg.norm(w.conj().T @ eta) ** 2
        remainder = _unit((np.eye(m) - w @ w.conj().T) @ eta)
        r = np.vdot(remainder, g_later @ remainder).real
        assert -1e-12 <= r <= 1 + 1e-12
        assert np.isclose(np.vdot(eta, g_later @ eta).real, q * np.vdot(phi, g @ phi).real + (1 - q) * r)


def test_without_intertwining_the_two_attributes_separate():
    g = np.diag([1.0, 0, 0])
    phi = np.array([1.0, 0, 0])
    w = np.zeros((7, 3))
    w[3:6, :] = np.eye(3)
    g_later = np.diag([1.0, 0, 0, 0, 0, 0, 0])
    eta = np.sqrt(0.9) * w @ phi + np.sqrt(0.1) * np.eye(7)[:, 6]
    q = np.linalg.norm(w.T @ eta) ** 2
    assert np.isclose(q, 0.9)
    assert np.isclose(eta @ g_later @ eta, 0.0)
    assert np.isclose(eta @ w @ g @ w.T @ eta / q, 1.0)
    assert np.linalg.norm(g_later @ w - w @ g) > 0.5
    assert eta @ g_later @ eta < q * (phi @ g @ phi)            # the mixing equality would require at least 0.9


def test_aspect_of_a_paired_state():
    rng = np.random.default_rng(42)
    for _ in range(1500):
        n = 7
        k_d, k_p = rng.integers(1, 6), rng.integers(1, 6)
        basis = _isometry(rng, n, k_d)
        projection_d = basis @ basis.conj().T
        frame = _isometry(rng, n, k_p)
        p = frame @ frame.conj().T
        phi = basis @ _unit(_vector(rng, k_d))
        weight = np.vdot(p @ phi, p @ phi).real
        assert np.isclose(np.vdot(phi, p @ phi).real, weight)
        assert np.linalg.norm(projection_d @ p @ phi) >= weight - 1e-12
        outside = (np.eye(n) - projection_d) @ p @ phi
        assert np.vdot(outside, outside).real <= weight * (1 - weight) + 1e-12
        if k_d == 1:
            assert np.isclose(np.vdot(outside, outside).real, weight * (1 - weight))
        np.testing.assert_allclose(outside + (np.eye(n) - projection_d) @ (np.eye(n) - p) @ phi, 0, atol=1e-12)
        # the weight is reproduced by an effect on the paired side
        coordinates = basis.conj().T @ phi
        assert np.isclose(np.vdot(coordinates, basis.conj().T @ p @ basis @ coordinates).real, weight)
