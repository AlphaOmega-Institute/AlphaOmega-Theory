"""Excitation energy and distance from a ground ray, for every state of finite energy.

    E_exc(psi) = || (H - E_0)^(1/2) psi ||^2 = (1 - F_x(psi)) Z(psi),
    F_x(psi) = |<x, psi>|^2,   Z(psi) = E_exc(eta_psi),   eta_psi = Q_x psi / sqrt(1 - F_x(psi)).

With a gap Delta above the ground on the complement of x:  E_exc >= Delta sin^2 d_FS([psi], [x]).
"""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import polar
from scipy.special import zeta

LABELS = {
    "TAUM": ["thm:reference-fidelity", "def:conditional-excitation-scale", "cor:fidelity-gap",
             "eq:Eexc-form-domain", "eq:FS-gap-bound", "eq:cyclic-reference"],
    "LG": ["eq:fidelity-state-decomposition", "prop:ground-fidelity-factorization", "cor:ground-gap-fidelity",
           "cor:BRLC-FS-isometry", "thm:projective-ground-recovery", "prop:no-state-below-ground"],
    "UM": ["eq:ground-return-energy", "eq:actual-ground-partner", "prop:ground-polar-compatibility",
           "prop:matched-born-probabilities", "prop:cyclic-determination-content", "sec:neutral-graphs-phase"],
}


def _unit(rng, n):
    v = rng.normal(size=n) + 1j * rng.normal(size=n)
    return v / np.linalg.norm(v)


def _hamiltonian(rng, n, degeneracy=1):
    """A Hermitian matrix with prescribed spectrum; returns it with its eigenvalues and eigenvectors."""
    levels = np.sort(rng.uniform(-2, 3, size=n - degeneracy + 1))
    levels = np.concatenate([np.full(degeneracy, levels[0]), levels[1:] + 0.05])
    basis, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    return (basis * levels) @ basis.conj().T, levels, basis


def _fidelity(x, psi):
    return abs(np.vdot(x, psi))**2


def _distance(x, psi):
    return np.arccos(min(1.0, abs(np.vdot(x, psi))))


def _excitation(H, E0, psi):
    return np.vdot(psi, (H - E0 * np.eye(len(psi))) @ psi).real


def test_factorization_gap_bound_and_decomposition():
    rng = np.random.default_rng(51)
    for _ in range(300):
        n = int(rng.integers(2, 8))
        H, levels, basis = _hamiltonian(rng, n)
        E0, x = levels[0], basis[:, 0]
        gap, width = levels[1] - E0, levels[-1] - E0
        psi = _unit(rng, n)
        F = _fidelity(x, psi)
        Q = np.eye(n) - np.outer(x, x.conj())
        energy = _excitation(H, E0, psi)
        assert energy >= -1e-12                                              # nothing lies below the ground
        root = (basis * np.sqrt(levels - E0)) @ basis.conj().T               # (H - E_0)^(1/2), zero on the ground
        assert_allclose(root @ x, 0, atol=1e-12)
        assert_allclose(np.linalg.norm(root @ psi)**2, energy, atol=1e-12)   # the form-domain expression
        assert_allclose(_excitation(H, E0, Q @ psi), energy, atol=1e-12)     # the ground component costs nothing
        eta = Q @ psi / np.sqrt(1 - F)
        assert_allclose(np.linalg.norm(eta), 1, atol=1e-12)
        assert_allclose(np.vdot(x, eta), 0, atol=1e-12)
        scale = _excitation(H, E0, eta)
        assert_allclose(energy, (1 - F) * scale, atol=1e-12)
        assert gap - 1e-12 <= scale <= width + 1e-12
        assert energy >= gap * (1 - F) - 1e-12
        assert_allclose(1 - F, np.sin(_distance(x, psi))**2, atol=1e-12)
        assert energy <= width * (1 - F) + 1e-12
        # reconstruction from the fidelity, one phase and the conditional vector
        phase = np.vdot(x, psi) / abs(np.vdot(x, psi))
        assert_allclose(phase * np.sqrt(F) * x + np.sqrt(1 - F) * eta, psi, atol=1e-12)
        # a state orthogonal to the ground needs no phase, and a state on the ground ray has no energy
        orthogonal = Q @ _unit(rng, n)
        orthogonal = orthogonal / np.linalg.norm(orthogonal)
        assert_allclose(_fidelity(x, orthogonal), 0, atol=1e-12)
        assert_allclose(_excitation(H, E0, orthogonal), _excitation(H, E0, Q @ orthogonal), atol=1e-12)
        assert_allclose(_excitation(H, E0, np.exp(0.7j) * x), 0, atol=1e-12)
    # with a degenerate ground level the identity holds for either ground vector and no gap bound does
    H, levels, basis = _hamiltonian(rng, 5, degeneracy=2)
    x, other = basis[:, 0], basis[:, 1]
    psi = _unit(rng, 5)
    F = _fidelity(x, psi)
    eta = (psi - x * np.vdot(x, psi)) / np.sqrt(1 - F)
    assert_allclose(_excitation(H, levels[0], psi), (1 - F) * _excitation(H, levels[0], eta), atol=1e-12)
    assert_allclose(_excitation(H, levels[0], other), 0, atol=1e-12)         # zero energy off the ray of x
    assert _fidelity(x, other) < 1e-12


def test_a_state_of_finite_energy_outside_the_operator_domain():
    # number operator, probabilities p_n = (n + 1)^(-5/2) / zeta(5/2)
    norm = zeta(2.5)
    fidelity = 1 / norm
    mean = zeta(1.5) / norm - 1
    assert_allclose(fidelity, 0.74544, atol=1e-5)
    assert_allclose(mean, 0.94737, atol=1e-5)
    n = np.arange(0, 2_000_000)
    p = (n + 1.0)**-2.5 / norm
    # the tail of the mean beyond M is about 2 / (zeta(5/2) sqrt M)
    assert abs(np.sum(n * p) - mean) < 2.5 / (norm * np.sqrt(len(n)))
    # the second moment has no limit: its partial sums grow like the square root of the cutoff
    second = np.cumsum(n**2 * p)
    assert 1.8 < second[1_999_999] / second[499_999] < 2.2
    assert second[1_999_999] > 1000
    # gap one: the bound, and the conditional scale
    assert mean >= 1 - fidelity
    assert_allclose(mean / (1 - fidelity), 3.7216, atol=1e-4)
    # the same identity on every truncation
    for cutoff in (10, 1000):
        amplitudes = np.sqrt(p[:cutoff] / p[:cutoff].sum())
        F = amplitudes[0]**2
        eta = amplitudes.copy()
        eta[0] = 0
        eta = eta / np.linalg.norm(eta)
        assert_allclose(np.sum(n[:cutoff] * amplitudes**2), (1 - F) * np.sum(n[:cutoff] * eta**2), atol=1e-12)


def test_fidelity_near_one_does_not_bound_an_unbounded_energy():
    # levels 0 (ground) and k^2 for k = 1, 2, ...; the state sqrt(1 - 1/k) x + phi_k / sqrt k
    for k in (2, 10, 100, 1000):
        fidelity = 1 - 1 / k
        energy = (1 / k) * k**2
        assert_allclose(energy, k)
        assert_allclose(energy, (1 - fidelity) * k**2)                      # the scale Z is k^2
        assert energy >= 1 * (1 - fidelity)                                 # the gap is one
    # with a gap, small energy forces a small distance from the ground ray
    rng = np.random.default_rng(52)
    H, levels, basis = _hamiltonian(rng, 6)
    E0, x = levels[0], basis[:, 0]
    gap = levels[1] - E0
    for size in (0.3, 0.03, 0.003):
        psi = x + size * _unit(rng, 6)
        psi = psi / np.linalg.norm(psi)
        energy = _excitation(H, E0, psi)
        assert np.sin(_distance(x, psi))**2 <= energy / gap + 1e-12


def test_path_examples():
    def path(n):
        H = np.zeros((n, n))
        for k in range(n - 1):
            H[k, k] += 0.5
            H[k + 1, k + 1] += 0.5
            H[k, k + 1] -= 0.5
            H[k + 1, k] -= 0.5
        return H

    H2, H5 = path(2), path(5)
    assert_allclose(np.linalg.eigvalsh(H2), [0, 1], atol=1e-12)
    x = np.ones(2) / np.sqrt(2)
    r = np.array([1, 2]) / np.sqrt(5)
    assert_allclose(_fidelity(x, r), 9 / 10)
    assert_allclose(_excitation(H2, 0, r), 1 / 10)
    rng = np.random.default_rng(53)
    for _ in range(20):                                   # one excited direction: the scale is the gap for every state
        psi = _unit(rng, 2)
        assert_allclose(_excitation(H2, 0, psi), 1 - _fidelity(x, psi), atol=1e-12)
    spectrum = np.linalg.eigvalsh(H5)
    assert_allclose(spectrum, 1 - np.cos(np.arange(5) * np.pi / 5), atol=1e-12)
    assert_allclose(spectrum[1], 0.19098, atol=1e-5)
    assert_allclose(spectrum[-1], 1.80902, atol=1e-5)
    y = np.ones(5) / np.sqrt(5)
    site = np.eye(5)[1]                                   # the second site
    F = _fidelity(y, site)
    energy = _excitation(H5, 0, site)
    assert_allclose(F, 1 / 5)
    assert_allclose(energy, 1)
    assert_allclose(energy / (1 - F), 5 / 4)
    assert spectrum[1] <= energy / (1 - F) <= spectrum[-1]


def test_fidelity_is_transported_by_the_polar_map_and_energy_is_not():
    rng = np.random.default_rng(54)
    for _ in range(100):
        m, n = int(rng.integers(2, 5)), int(rng.integers(5, 8))
        A = rng.normal(size=(n, m)) + 1j * rng.normal(size=(n, m))
        V, modulus = polar(A)                                               # A = V |A|, V an isometry
        x = _unit(rng, m)
        y = A @ x / np.linalg.norm(A @ x)                                   # weak bridge
        r = V.conj().T @ y                                                  # the partner of y
        assert_allclose(r, modulus @ x / np.linalg.norm(A @ x), atol=1e-10)
        psi, chi = _unit(rng, m), _unit(rng, m)
        assert_allclose(_fidelity(y, V @ psi), _fidelity(r, psi), atol=1e-12)
        # Fubini-Study distances between the rays, their images and their graph vectors
        distance = _distance(psi, chi)
        assert_allclose(_distance(V @ psi, V @ chi), distance, atol=1e-10)
        for phase in (0.0, 1.3, np.pi):
            graph = lambda v: np.concatenate([v, np.exp(1j * phase) * V @ v]) / np.sqrt(2)
            assert_allclose(np.linalg.norm(graph(psi)), 1, atol=1e-12)
            assert_allclose(_distance(graph(psi), graph(chi)), distance, atol=1e-10)
    # energies are transported only by an intertwining map: the unequal-cell path
    copying = np.repeat(np.eye(2), [1, 4], axis=0) / 2
    V, _ = polar(copying)
    H2 = np.array([[0.5, -0.5], [-0.5, 0.5]])
    H5 = np.zeros((5, 5))
    for k in range(4):
        H5[k:k + 2, k:k + 2] += np.array([[0.5, -0.5], [-0.5, 0.5]])
    x = np.ones(2) / np.sqrt(2)
    assert not np.allclose(H5 @ V, V @ H2)
    assert_allclose(x @ H2 @ x, 0, atol=1e-12)
    assert_allclose((V @ x) @ H5 @ (V @ x), 1 / 16, atol=1e-12)             # the image of the coarse ground has energy


def test_recovery_of_the_ground_ray_along_a_sequence():
    rng = np.random.default_rng(55)
    n = 5
    y = _unit(rng, n)
    B = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    previous = np.inf
    for size in (1.0, 0.1, 0.01, 0.001):
        X = np.eye(n) + size * B                                            # determinations approaching the identity
        mean = np.vdot(y, X @ y)                                            # omega(X)
        centered = X @ y - mean * y                                         # X° y, orthogonal to y
        assert_allclose(np.vdot(y, centered), 0, atol=1e-12)
        ratio = np.linalg.norm(centered) / abs(mean)
        state = X @ y / np.linalg.norm(X @ y)
        F = _fidelity(y, state)
        assert_allclose(F, 1 / (1 + ratio**2), atol=1e-12)
        assert_allclose(F, abs(mean)**2 / np.vdot(X @ y, X @ y).real, atol=1e-12)
        assert_allclose(np.cos(_distance(y, state)), np.sqrt(F), atol=1e-12)
        assert _distance(y, state) < previous
        previous = _distance(y, state)
    assert previous < 5e-3
    # a sequence whose ratio does not tend to zero stays away from the ground ray
    for scale in (1.0, 1e-3, 1e-6):
        X = scale * (np.eye(n) + B)
        mean = np.vdot(y, X @ y)
        ratio = np.linalg.norm(X @ y - mean * y) / abs(mean)
        assert_allclose(_fidelity(y, X @ y / np.linalg.norm(X @ y)), 1 / (1 + ratio**2), atol=1e-12)
        assert ratio > 0.1


def test_the_return_operator_as_the_hamiltonian():
    rng = np.random.default_rng(56)
    for _ in range(100):
        m = int(rng.integers(2, 6))
        n = m + int(rng.integers(0, 3))                                     # at least m - 1 rows, so the kernel is the line of x
        x = _unit(rng, m)
        A = rng.normal(size=(n, m)) + 1j * rng.normal(size=(n, m))
        A = A @ (np.eye(m) - np.outer(x, x.conj()))                         # kernel spanned by x
        assert np.linalg.matrix_rank(A) == m - 1
        c = np.linalg.svd(A, compute_uv=False)[m - 2]                       # lower bound on the complement of the kernel
        psi = _unit(rng, m)
        energy = np.linalg.norm(A @ psi)**2                                 # the form of A*A, with E_0 = 0
        assert_allclose(energy, _excitation(A.conj().T @ A, 0, psi), atol=1e-10)
        assert energy >= c**2 * np.sin(_distance(x, psi))**2 - 1e-10
    # under strong compatibility x is an eigenvector of A*A, not necessarily its lowest
    A = np.diag([2.0, 1.0])
    x = np.array([1.0, 0.0])
    assert_allclose(A.T @ A @ x, 4 * x)
    assert np.linalg.eigvalsh(A.T @ A).min() < 4
