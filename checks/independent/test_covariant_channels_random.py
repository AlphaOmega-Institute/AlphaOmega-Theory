"""Random time-covariant channels in the shift model of UM Section 7 (hbar = 1).

Input levels: x = f0 at energy 0 and the remaining inputs at 1.  Output levels: the
missing line z = g0 at epsilon, g1 at 0 and the remaining outputs at 1.  A covariant
channel is one whose Kraus operators each have a single energy shift.  The committed
check uses one channel with a positive amplitude on every allowed transition.  Here the
Kraus operators are random within each shift class, with several operators per class.
"""
import numpy as np
from scipy.linalg import expm

LABELS = {
    "UM": ["eq:index-channel-coherence"],
}

N = 5                                                          # input dimension; the output has N + 1


def _levels(epsilon):
    return np.r_[0.0, np.ones(N - 1)], np.r_[epsilon, 0.0, np.ones(N - 1)]


def _covariant_kraus(rng, epsilon, preserve_ground, per_class=3):
    e_in, e_out = _levels(epsilon)
    shift = e_out[:, None] - e_in[None, :]
    operators = []
    for value in np.unique(np.round(shift, 9)):
        mask = np.isclose(shift, value)
        if preserve_ground:
            mask[:, 0] &= np.arange(N + 1) == 1                # x may only go to g1
        if mask.any():
            for _ in range(per_class):
                operators.append((rng.normal(size=mask.shape) + 1j * rng.normal(size=mask.shape)) * mask)
    total = sum(k.conj().T @ k for k in operators)
    values, vectors = np.linalg.eigh(total)
    fix = vectors @ np.diag(values**-0.5) @ vectors.conj().T    # commutes with the input Hamiltonian
    return [k @ fix for k in operators]


def _apply(kraus, operator):
    return sum(k @ operator @ k.conj().T for k in kraus)


def _state(rng):
    raw = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    rho = raw @ raw.conj().T
    return rho / np.trace(rho).real


def test_random_covariant_channels_and_the_missing_line():
    rng = np.random.default_rng(11)
    for epsilon, preserve, expect_coherence in ((0.37, False, False), (0.37, True, False), (2.0, False, True),
                                                (2.0, True, False), (1.0, False, True), (1.0, True, True)):
        e_in, e_out = _levels(epsilon)
        largest = 0.0
        for _ in range(25):
            kraus = _covariant_kraus(rng, epsilon, preserve)
            np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(N), atol=1e-12)
            rho = _state(rng)
            out = _apply(kraus, rho)
            t = rng.uniform(0, 7)
            u_in, u_out = expm(-1j * t * np.diag(e_in)), expm(-1j * t * np.diag(e_out))
            np.testing.assert_allclose(u_out @ out @ u_out.conj().T, _apply(kraus, u_in @ rho @ u_in.conj().T), atol=1e-12)
            largest = max(largest, np.abs(out[0, 1:]).max())
            # the line weight depends only on the two diagonal energy blocks of the input
            blocks = rho.copy()
            blocks[0, 1:] = 0
            blocks[1:, 0] = 0
            assert np.isclose(out[0, 0], _apply(kraus, blocks)[0, 0], atol=1e-13)
            if preserve:
                ground = np.zeros((N, N))
                ground[0, 0] = 1
                target = np.zeros((N + 1, N + 1))
                target[1, 1] = 1
                np.testing.assert_allclose(_apply(kraus, ground), target, atol=1e-12)
        assert (largest > 1e-3) == expect_coherence


def test_two_factor_bound_without_covariance():
    rng = np.random.default_rng(42)
    environment = 3
    q = np.eye(N + 1)
    q[0, 0] = 0
    x = np.eye(N)[:, 0].astype(complex)
    closest = 0.0
    for _ in range(1500):
        raw = rng.normal(size=((N + 1) * environment, N)) + 1j * rng.normal(size=((N + 1) * environment, N))
        isometry, _ = np.linalg.qr(raw)
        kraus = [isometry.reshape(N + 1, environment, N)[:, j, :] for j in range(environment)]
        f = np.zeros(N, complex)
        f[1:] = rng.normal(size=N - 1) + 1j * rng.normal(size=N - 1)
        f /= np.linalg.norm(f)
        coherence = np.linalg.svd(q @ _apply(kraus, np.outer(f, x.conj())) @ q, compute_uv=False).sum()
        p_f = _apply(kraus, np.outer(f, f.conj()))[0, 0].real
        p_x = _apply(kraus, np.outer(x, x.conj()))[0, 0].real
        bound = np.sqrt((1 - p_f) * (1 - p_x))
        assert coherence <= bound + 1e-12
        closest = max(closest, coherence / bound)
    assert closest > 0.3
