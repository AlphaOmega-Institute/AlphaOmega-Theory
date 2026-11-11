"""Checks for the exchange self-channels and the control/mixer analysis of the
pair model (LG, the two unnumbered subsections following the pair realization)."""
import numpy as np
from scipy.linalg import expm, null_space

from lib.fock import dag, number
from lib.linalg import rank, response_weights, in_span
from pair.model import PairModel

LABELS = {
    "LG": [
        "eq:pair-completed-hamiltonian", "eq:pair-completed-spectrum",
        "eq:pair-completed-response", "eq:pair-general-mixer",
        "eq:pair-mixer-ground-certificate", "prop:pair-invariant-mixer",
    ],
}

M = PairModel(theta=0.3, j=1.0)


def _mixer(K):
    TK = sum(K[i, k] * (dag(M.c[i]) @ M.d[k]) for i in range(2) for k in range(2))
    return (TK - dag(TK)) / 2j


def test_exchange_completion_spectrum_and_response():
    """H_s g = -(s+j) g, H_s e = (j-s) e, H_s m = 0 = H_s|_W; correlation
    u^2 exp(-2ijt) + v^2 exp(-i(s+j)t); defect weight sin^2(4 theta) sin^2((j-s)t/2);
    [Y, H_s] e = i (j-s) m."""
    j = M.j
    for s in (0.0, 0.4, 1.0, 1.7):
        Hs = M.exchange_hamiltonian(s)
        assert np.allclose(Hs @ M.g, -(s + j) * M.g)
        assert np.allclose(Hs @ M.e, (j - s) * M.e)
        assert np.allclose(Hs @ M.m, 0) and np.allclose(Hs @ M.W, 0)
        spectrum = np.linalg.eigvalsh(M.F2.conj().T @ Hs @ M.F2)
        assert np.sum(np.isclose(spectrum, spectrum[0])) == 1          # unique ground
        for t in (0.5, 1.9):
            U = expm(-1j * Hs * t)
            corr = M.g @ dag(U) @ M.Ztheta @ U @ M.Ztheta @ M.g
            assert np.isclose(corr, M.u**2 * np.exp(-2j * j * t) + M.v**2 * np.exp(-1j * (s + j) * t))
            assert np.isclose(abs(M.dvec @ (U @ M.b)) ** 2,
                              np.sin(4 * M.theta) ** 2 * np.sin((j - s) * t / 2) ** 2)
        assert np.allclose((M.Y @ Hs - Hs @ M.Y) @ M.e, 1j * (j - s) * M.m)


def test_minimal_completion_removes_the_accessible_defect():
    """At s = j: H = -j (P_c + P_d)^*(P_c + P_d), spectrum {-2j, 0 x 5}; the prepared
    range is invariant, the response-generated defect vanishes, and H commutes with Y.
    For s != j the defect count is one."""
    j = M.j
    Hmin = M.exchange_hamiltonian(j)
    S = M.Pc + M.Pd
    assert np.allclose(Hmin, -j * dag(S) @ S)
    assert np.allclose(np.linalg.eigvalsh(M.F2.conj().T @ Hmin @ M.F2), [-2 * j, 0, 0, 0, 0, 0], atol=1e-12)
    assert np.allclose(M.Y @ Hmin, Hmin @ M.Y)
    range_basis = np.column_stack([M.g, M.b])
    for vec in (M.g, M.b):
        assert in_span(Hmin @ vec, range_basis)
    for s, expected in ((0.0, 1), (j, 0), (1.5 * j, 1)):
        H2 = M.F2.conj().T @ M.exchange_hamiltonian(s) @ M.F2
        weights = [R for _, R in response_weights(H2, M.F2.conj().T @ M.A)]
        assert sum(rank(R) for R in weights) - rank(sum(weights)) == expected
    # on the fixed three-dimensional carrier the ground is no longer cyclic at s = j
    Hc, Zc = M.on_carrier(Hmin), M.on_carrier(M.Ztheta)
    g = M.Phi.T @ M.g
    orbit = np.column_stack([g, Zc @ g, Hc @ Zc @ g, Zc @ Hc @ Zc @ g, Hc @ Hc @ Zc @ g])
    assert rank(orbit) == 2


def test_asymmetric_completion_moves_the_ground():
    """For self-channels s_c, s_d the paired ground is alpha C + beta D with
    Y (alpha C + beta D) = i (alpha - beta) m / sqrt 2; alpha = beta exactly when s_c = s_d."""
    for s_c, s_d in ((1.5, 1.5), (1.5, 1.2), (2.0, 1.1)):
        _, V = np.linalg.eigh(-np.array([[s_c, M.j], [M.j, s_d]]))
        alpha, beta = V[:, 0] * np.sign(V[0, 0])
        ground = alpha * M.C + beta * M.D
        assert np.allclose(M.exchange_hamiltonian(s_c, s_d) @ ground,
                           (ground @ M.exchange_hamiltonian(s_c, s_d) @ ground) * ground)
        assert np.allclose(M.Y @ ground, 1j * (alpha - beta) * M.m / np.sqrt(2))
        assert np.isclose(alpha, beta) == (s_c == s_d)


def test_pair_controls_preserve_orbital_parity():
    """The pair Hamiltonian, both quadratures of P_c^*P_d, the imbalance Z and the
    self-channels commute with (-1)^{N_c} and preserve span{C, D}; m has the opposite parity."""
    Nc = number(M.c[0]) + number(M.c[1])
    parity = np.diag(np.where(np.round(np.diag(Nc)) % 2 == 0, 1.0, -1.0))
    hop = dag(M.Pc) @ M.Pd
    controls = [M.H, hop + dag(hop), 1j * (hop - dag(hop)), M.Z, dag(M.Pc) @ M.Pc, dag(M.Pd) @ M.Pd]
    paired = np.column_stack([M.C, M.D])
    for op in controls:
        assert np.allclose(op @ parity, parity @ op)
        assert in_span(op @ M.C, paired) and in_span(op @ M.D, paired)
    assert np.allclose(parity @ M.C, M.C) and np.allclose(parity @ M.D, M.D)
    assert np.allclose(parity @ M.m, -M.m)
    assert not np.allclose(M.Y @ parity, parity @ M.Y)       # the one-body mixer is an extra operation


def test_mixer_ground_certificate():
    """||Y_K g||^2 = (|K11 - conj K22|^2 + |K12 + conj K21|^2)/4 and Y_K g is orthogonal
    to g; it vanishes exactly for K = [[z, w], [-conj w, conj z]]."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        K = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        YK = _mixer(K)
        expected = (abs(K[0, 0] - np.conj(K[1, 1])) ** 2 + abs(K[0, 1] + np.conj(K[1, 0])) ** 2) / 4
        assert np.isclose(np.linalg.norm(YK @ M.g) ** 2, expected)
        assert np.isclose(np.vdot(M.g, YK @ M.g), 0)
    z, w = 0.3 + 0.4j, -0.2 + 0.7j
    assert np.allclose(_mixer(np.array([[z, w], [-np.conj(w), np.conj(z)]])) @ M.g, 0)
    assert np.allclose(_mixer(np.eye(2)), M.Y)


def test_invariant_mixer_needs_equivalent_representations():
    """U(1) charges (+1,-1) on c and (+2,-2) on d: the pair term is neutral but no
    nonzero one-body mixer between the blocks is invariant.  With equal charges the
    aligned mixer K = I is invariant."""
    def charge(qc, qd):
        return sum(q * number(a) for q, a in zip(qc + qd, M.c + M.d))

    hop = dag(M.Pc) @ M.Pd
    units = [np.eye(2)[:, [i]] @ np.eye(2)[[k], :] for i in range(2) for k in range(2)]

    def invariant_mixer_dimension(Q):
        cols = []
        for E in units:
            TK = sum(E[i, k] * (dag(M.c[i]) @ M.d[k]) for i in range(2) for k in range(2))
            cols.append((Q @ TK - TK @ Q).ravel())
        return 4 - rank(np.array(cols).T)

    Q = charge([1, -1], [2, -2])
    assert np.allclose(Q @ hop, hop @ Q)                 # determinant characters agree
    assert invariant_mixer_dimension(Q) == 0
    Q_equal = charge([1, -1], [1, -1])
    assert invariant_mixer_dimension(Q_equal) == 2       # diagonal K, including K = I
    assert np.allclose(Q_equal @ M.T, M.T @ Q_equal)


def test_unitary_intertwiner_gives_a_ground_preserving_mixer():
    """If R_d = V^* R_c V for a unitary V, then K = V / sqrt(det V) intertwines and Y_K g = 0."""
    rng = np.random.default_rng(1)
    for _ in range(50):
        X = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        V, _ = np.linalg.qr(X)
        Y_, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        Rc = Y_
        Rd = V.conj().T @ Rc @ V
        K = V / np.sqrt(np.linalg.det(V))
        assert np.allclose(Rc @ K, K @ Rd)
        assert np.allclose(_mixer(K) @ M.g, 0)


def test_unused_space_is_not_protected_by_number_conservation():
    """X, Y, Z annihilate W, but the number-conserving operator
    n_c1 - n_c2 - n_d1 + n_d2 maps w_0 to 2m (and m to 2 w_0)."""
    for op in (M.X, M.Y, M.Z):
        assert np.allclose(op @ M.W, 0)
    O = number(M.c[0]) - number(M.c[1]) - number(M.d[0]) + number(M.d[1])
    assert np.allclose(O @ M.w0, 2 * M.m) and np.allclose(O @ M.m, 2 * M.w0)
    assert np.allclose(O @ M.N, M.N @ O)
