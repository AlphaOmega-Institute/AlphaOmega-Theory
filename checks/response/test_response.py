"""Checks for the port-response statements of TAUM: the finite response dimension
count, permanent invisibility to a fixed port, the symmetry-resolved response, and the
three-state example.  Finite-dimensional instances only; hbar = 1."""
import numpy as np
from scipy.linalg import expm, null_space, orth

from lib.linalg import rank, response_weights, minimal_space, spectral_blocks

LABELS = {
    "TAUM": [
        "def:response-generated-subspace", "eq:minimal-subspace", "thm:canonical-minimal-reduction",
        "cor:finite-response-defect", "eq:finite-response-weights", "eq:indirect-response-space",
        "eq:finite-response-defect-count", "eq:ambient-response-defect",
        "prop:permanent-port-invisibility", "sec:symmetry-response", "prop:opposite-sector-response",
        "eq:opposite-sector-measure", "eq:opposite-sector-port", "eq:opposite-sector-factorization",
        "eq:opposite-sector-departure", "eq:opposite-sector-fixed-line",
        "sec:three-state-response", "ex:split-excitation-response", "eq:three-state-data",
        "eq:three-state-port", "eq:three-state-response", "eq:three-state-weights",
        "eq:three-state-minimal-defect", "eq:three-state-transfer",
    ],
}


def _random_realization(rng, dim, multiplicities, ports, port_rank):
    """Hermitian H with prescribed eigenvalue multiplicities and a port of given rank."""
    assert sum(multiplicities) == dim
    eigenvalues = np.repeat(np.arange(len(multiplicities), dtype=float) * 0.7 - 1.0, multiplicities)
    Q, _ = np.linalg.qr(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    H = Q @ np.diag(eigenvalues) @ Q.conj().T
    left = rng.normal(size=(dim, port_rank)) + 1j * rng.normal(size=(dim, port_rank))
    right = rng.normal(size=(port_rank, ports)) + 1j * rng.normal(size=(port_rank, ports))
    return (H + H.conj().T) / 2, left @ right


def test_finite_response_dimension_count():
    """dim H_min = sum_j rank R_j, dim(H_min minus Ran Gamma) = sum_j rank R_j - rank(sum_j R_j),
    Gamma^* Gamma = sum_j R_j, and Ker Gamma^* = D_Gamma + (X minus H_min)."""
    rng = np.random.default_rng(7)
    cases = [(6, (1, 2, 3), 2, 2), (7, (2, 2, 3), 3, 2), (8, (1, 1, 3, 3), 2, 1), (9, (3, 3, 3), 4, 3)]
    for dim, multiplicities, ports, port_rank in cases:
        H, Gamma = _random_realization(rng, dim, multiplicities, ports, port_rank)
        weights = [R for _, R in response_weights(H, Gamma)]
        Hmin = minimal_space(H, Gamma)
        assert Hmin.shape[1] == sum(rank(R) for R in weights)
        assert np.allclose(sum(weights), Gamma.conj().T @ Gamma)
        range_dim = rank(Gamma)
        assert range_dim == rank(sum(weights))
        defect = Hmin.shape[1] - range_dim
        assert defect == sum(rank(R) for R in weights) - rank(sum(weights))
        # H_min reduces H and contains the range of Gamma
        P = Hmin @ Hmin.conj().T
        assert np.allclose(P @ H, H @ P) and np.allclose(P @ Gamma, Gamma)
        # ambient kernel: dim Ker Gamma^* = defect + (dim X - dim H_min)
        assert null_space(Gamma.conj().T).shape[1] == defect + (dim - Hmin.shape[1])


def test_relations_among_ranges_without_pairwise_overlap():
    """Three weights whose ranges are three distinct lines in a plane: every pairwise
    intersection is zero, yet the count is 3 - 2 = 1."""
    lines = [np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([1.0, 1.0]) / np.sqrt(2)]
    weights = [np.outer(x, x) for x in lines]
    assert sum(rank(R) for R in weights) - rank(sum(weights)) == 1
    # a realization with these weights: H = diag(0, 1, 2), Gamma rows = the three lines
    H, Gamma = np.diag([0.0, 1.0, 2.0]), np.array(lines)
    for (lam, R), expected in zip(response_weights(H, Gamma), weights):
        assert np.allclose(R, expected)
    assert minimal_space(H, Gamma).shape[1] == 3 and rank(Gamma) == 2


def test_permanent_invisibility_is_the_complement_of_the_response_space():
    """{x : Gamma^* exp(-itH) x = 0 for all t} = H_min^perp.  In finite dimensions this is
    the null space of the stacked matrices Gamma^* H^k."""
    rng = np.random.default_rng(11)
    for dim, multiplicities, ports, port_rank in [(7, (2, 2, 3), 2, 2), (8, (1, 3, 4), 3, 2)]:
        H, Gamma = _random_realization(rng, dim, multiplicities, ports, port_rank)
        Hmin = minimal_space(H, Gamma)
        stacked = np.vstack([Gamma.conj().T @ np.linalg.matrix_power(H, k) for k in range(dim)])
        invisible = null_space(stacked, rcond=1e-9)
        assert invisible.shape[1] == dim - Hmin.shape[1]
        assert np.allclose(Hmin.conj().T @ invisible, 0, atol=1e-8)
        for t in (0.3, 2.0, 17.0):
            assert np.allclose(Gamma.conj().T @ expm(-1j * t * H) @ invisible, 0, atol=1e-8)
        # a vector with a component in H_min is seen at some time
        x = Hmin[:, 0]
        seen = max(np.linalg.norm(Gamma.conj().T @ expm(-1j * t * H) @ x) for t in np.linspace(0, 10, 101))
        assert seen > 1e-3


def _opposite_sectors(branch_energies, nu, ell=1.0, E0=0.0):
    """g with L = 0, and q_+ / q_- in the L = +-ell sectors with the same H_0 spectral measure."""
    k = len(branch_energies)
    H0 = np.diag([E0] + list(branch_energies) * 2)
    L = np.diag([0.0] + [ell] * k + [-ell] * k)
    weights = np.random.default_rng(1).uniform(0.5, 1.5, k)
    weights /= np.linalg.norm(weights)
    qp, qm = np.zeros(1 + 2 * k), np.zeros(1 + 2 * k)
    qp[1:1 + k], qm[1 + k:] = weights, weights
    g = np.zeros(1 + 2 * k)
    g[0] = 1
    e = (qp + qm) / np.sqrt(2)
    d = (qp - qm) / (1j * np.sqrt(2))
    return H0 - nu * L, g, e, d, weights


def test_symmetry_resolved_response():
    """C_nu(t) = C_0(t) cos(nu ell t); weight on the fixed line d is |C_0|^2 sin^2;
    weight outside the port range within H_min is 1 - |C_0|^2 cos^2.  When the shifted
    branches overlap, d need not lie in H_min."""
    nu = 0.2
    for energies, d_reachable in (([1.0, 3.0], True), ([1.0, 1.4], False)):
        H, g, e, d, weights = _opposite_sectors(energies, nu)
        Gamma = np.column_stack([g, e])
        Hmin = minimal_space(H, Gamma)
        outside = Hmin @ Hmin.conj().T - Gamma @ Gamma.conj().T
        for t in (0.3, 1.1, 2.9, np.pi / (2 * nu)):
            psi = expm(-1j * t * H) @ e
            C0 = np.sum(weights**2 * np.exp(-1j * t * np.array(energies)))
            assert np.isclose(np.vdot(e, psi), C0 * np.cos(nu * t))
            assert np.isclose(abs(np.vdot(d, psi)) ** 2, abs(C0) ** 2 * np.sin(nu * t) ** 2)
            assert np.isclose(np.vdot(psi, outside @ psi).real, 1 - abs(C0) ** 2 * np.cos(nu * t) ** 2)
        assert np.isclose(np.linalg.norm(Hmin.conj().T @ d), 1) == d_reachable


def test_three_state_example():
    """Three atoms with weights diag(1,0), diag(0,1/2), diag(0,1/2); H_min = X; the
    response defect is the line through d; exp(-it(H-E0)) e = exp(-it(E2-E0)) (cos e - sin d);
    the compression Gamma^* H Gamma is diag(E0, E2); at zero bias the defect vanishes."""
    E0, E2, ell, nu = 0.0, 1.0, 1.0, 0.3
    g, fp, fm = np.eye(3)
    H = np.diag([E0, E2 - nu * ell, E2 + nu * ell])
    e = (fp + fm) / np.sqrt(2)
    d = (fp - fm) / (1j * np.sqrt(2))
    Gamma = np.column_stack([g, e])
    atoms = response_weights(H, Gamma)
    assert [round(lam, 9) for lam, _ in atoms] == [E0, E2 - nu * ell, E2 + nu * ell]
    assert np.allclose(atoms[0][1], np.diag([1.0, 0.0]))
    assert np.allclose(atoms[1][1], np.diag([0.0, 0.5])) and np.allclose(atoms[2][1], np.diag([0.0, 0.5]))
    assert np.allclose(sum(R for _, R in atoms), np.eye(2))
    assert minimal_space(H, Gamma).shape[1] == 3
    kernel = null_space(Gamma.conj().T)
    assert kernel.shape[1] == 1 and np.isclose(abs(np.vdot(kernel[:, 0], d)), 1)
    for t in (0.8, 4.0):
        psi = expm(-1j * t * (H - E0 * np.eye(3))) @ e
        assert np.allclose(psi, np.exp(-1j * t * (E2 - E0)) * (np.cos(nu * ell * t) * e - np.sin(nu * ell * t) * d))
    assert np.allclose(Gamma.conj().T @ H @ Gamma, np.diag([E0, E2]))
    assert not np.allclose(H @ Gamma, Gamma @ (Gamma.conj().T @ H @ Gamma))     # port not invariant
    H_zero_bias = np.diag([E0, E2, E2])
    assert minimal_space(H_zero_bias, Gamma).shape[1] == 2


def test_intertwining_criterion_and_detuning_cap():
    """For an isometric port A with P = AA^*: prepared states never leave the range exactly
    when (1 - P) H A = 0, equivalently H A = A (A^* H A).  With coupling kappa between a
    prepared line and a defect line detuned by Delta, the maximal defect weight is
    4 kappa^2 / (4 kappa^2 + Delta^2).  (Supports the correspondence-versus-intertwining
    remark; that remark is not yet in the manuscripts.)"""
    rng = np.random.default_rng(0)
    A = np.array([[1, 0], [0, 1], [0, 0]], dtype=complex)
    Q = np.eye(3) - A @ A.conj().T
    for trial in range(200):
        X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        H = (X + X.conj().T) / 2
        if trial % 2:
            H[2, :2] = 0
            H[:2, 2] = 0
        decoupled = np.linalg.norm(Q @ H @ A) < 1e-12
        assert decoupled == np.allclose(H @ A, A @ (A.conj().T @ H @ A))
        if decoupled:
            for t in (0.7, 3.0):
                assert np.allclose(Q @ expm(-1j * t * H) @ A, 0, atol=1e-10)
    kappa = 0.14
    times = np.linspace(0, 400, 40001)
    for Delta in (0.0, 0.1, 0.28, 1.0):
        H = np.array([[0, 0, 0], [0, 1.0, 1j * kappa], [0, -1j * kappa, 1.0 + Delta]])
        w, V = np.linalg.eigh(H)
        start = V.conj().T @ np.array([0, 1, 0], dtype=complex)
        defect_amplitude = (V[2, :] * np.exp(-1j * np.outer(times, w))) @ start
        assert np.isclose(np.max(np.abs(defect_amplitude) ** 2),
                          4 * kappa**2 / (4 * kappa**2 + Delta**2), atol=1e-4)
