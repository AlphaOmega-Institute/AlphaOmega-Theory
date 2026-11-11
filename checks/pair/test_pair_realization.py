"""Checks for the finite pair-transfer realization (LG, sec:finite-pair-realization).

Each check rebuilds the stated objects as explicit matrices on the four-mode
Fock space and asserts the identities the manuscript displays.
"""
import numpy as np
from scipy.linalg import expm, null_space

from lib.fock import dag, satisfies_car
from lib.linalg import rank, response_weights, minimal_space, generated_algebra_dim
from pair.model import PairModel

LABELS = {
    "LG": [
        "sec:finite-pair-realization", "eq:pair-parent", "eq:pair-energy-vectors",
        "eq:pair-orbital-rotation", "eq:pair-preparation-map", "eq:pair-imbalance",
        "eq:pair-bright-defect", "eq:pair-channel", "eq:pair-port-weights",
        "prop:pair-ground-defect", "eq:pair-readout-action", "eq:pair-defect-evolution",
        "eq:pair-defect-measurement", "eq:pair-ambient-zero",
        "eq:pair-two-frequency-response", "eq:pair-mixer-selection",
    ],
}

M = PairModel(theta=0.3, j=1.0)


def test_canonical_anticommutation_relations():
    assert satisfies_car(M.c + M.d)


def test_spectrum_and_unique_ground():
    """H has spectrum {-j, 0, 0, 0, 0, +j} on the two-particle sector; g, e, m are
    orthonormal with Hg = -jg, He = je, Hm = 0; the ground is unique there."""
    H2 = M.F2.conj().T @ M.H @ M.F2
    assert np.allclose(np.linalg.eigvalsh(H2), [-M.j, 0, 0, 0, 0, M.j], atol=1e-12)
    for vec, energy in ((M.g, -M.j), (M.e, M.j), (M.m, 0.0)):
        assert np.isclose(np.linalg.norm(vec), 1)
        assert np.allclose(M.H @ vec, energy * vec)
    assert np.allclose([M.g @ M.e, M.g @ M.m, M.e @ M.m], 0)
    assert np.allclose(M.H @ M.W, 0)


def test_preparation_channel():
    """A g_Psi = g, A e_Psi = b, A^*A = I, and Ker A^* inside span{g,e,m} is the line C d."""
    assert np.allclose(M.A @ M.gPsi, M.g)
    assert np.allclose(M.A @ M.ePsi, M.b)
    assert np.allclose(M.A.T @ M.A, np.eye(2))
    Ac = M.Phi.T @ M.A                      # channel with codomain span{g, e, m}
    assert np.allclose(M.Phi @ Ac, M.A)     # the range lies in the carrier
    kernel = null_space(Ac.T)
    assert kernel.shape[1] == 1
    assert np.isclose(abs(kernel[:, 0] @ (M.Phi.T @ M.dvec)), 1)


def test_readout_on_carrier():
    """Z_theta g = b, Z_theta b = g, Z_theta d = 0, hence Z_theta^2 = AA^* on the carrier."""
    assert np.allclose(M.Ztheta @ M.g, M.b)
    assert np.allclose(M.Ztheta @ M.b, M.g)
    assert np.allclose(M.Ztheta @ M.dvec, 0)
    Zc, Ac = M.on_carrier(M.Ztheta), M.Phi.T @ M.A
    assert np.allclose(Zc @ Zc, Ac @ Ac.T)


def test_readout_on_full_two_particle_space():
    """On F_2: Z_theta^2 is the projection onto the prepared range, and
    1 - Z_theta^2 = |d><d| + P_W; the zero eigenspace is four-dimensional."""
    Z2 = M.P2 @ M.Ztheta @ M.Ztheta @ M.P2
    assert np.allclose(Z2, M.A @ M.A.T)
    PW = M.P2 - M.Phi @ M.Phi.T
    assert np.allclose(PW, M.W @ M.W.T)
    assert np.allclose(M.P2 - Z2, np.outer(M.dvec, M.dvec) + PW)
    spectrum = np.linalg.eigvalsh(M.F2.T @ M.Ztheta @ M.F2)
    assert np.allclose(spectrum, [-1, 0, 0, 0, 0, 1], atol=1e-12)
    assert np.allclose(M.Ztheta @ M.W, 0) and np.allclose(M.H @ M.W, 0)


def test_zero_outcome_probability_for_ambient_states():
    """Pr(Z=0) = weight on the line C d + weight in W, for arbitrary states on F_2.
    Relative to the state conditioned on the carrier, the excess is eps (1 - w) <= eps."""
    rng = np.random.default_rng(5)
    P0 = M.P2 - M.P2 @ M.Ztheta @ M.Ztheta @ M.P2
    PW = M.W @ M.W.T
    Pd = np.outer(M.dvec, M.dvec)
    for _ in range(500):
        X = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
        s6 = X @ X.conj().T
        sigma = M.F2 @ (s6 / np.trace(s6)) @ M.F2.conj().T
        pr0, eps, line = (np.trace(sigma @ P).real for P in (P0, PW, Pd))
        assert np.isclose(pr0, line + eps)
        w = line / (1 - eps)
        assert np.isclose(pr0 - w, eps * (1 - w)) and -1e-12 <= pr0 - w <= eps + 1e-12


def test_port_weights_recover_the_carrier():
    """R_{-j} = P_g, R_0 = v^2 P_e, R_{+j} = u^2 P_e; the three ranges span {g, e, m},
    and the rank count gives a one-dimensional response defect."""
    H2 = M.F2.conj().T @ M.H @ M.F2
    A2 = M.F2.conj().T @ M.A
    Pg, Pe = np.outer(M.gPsi, M.gPsi), np.outer(M.ePsi, M.ePsi)
    weights = dict((round(lam, 9), R) for lam, R in response_weights(H2, A2))
    assert np.allclose(weights[-M.j], Pg)
    assert np.allclose(weights[0.0], M.v**2 * Pe)
    assert np.allclose(weights[M.j], M.u**2 * Pe)
    Hmin = M.F2 @ minimal_space(H2, A2)
    assert Hmin.shape[1] == 3
    assert np.allclose(Hmin @ Hmin.conj().T, M.Phi @ M.Phi.T)
    count = sum(rank(R) for R in weights.values()) - rank(sum(weights.values()))
    assert count == 1


def test_defect_weight_and_two_line_correlation():
    """Defect weight sin^2(4 theta) sin^2(jt/2); ground correlation
    u^2 exp(-2ijt) + v^2 exp(-ijt)."""
    for t in (0.4, 1.7, np.pi):
        U = expm(-1j * M.H * t)
        assert np.isclose(abs(M.dvec @ (U @ M.b)) ** 2,
                          np.sin(4 * M.theta) ** 2 * np.sin(M.j * t / 2) ** 2)
        Zt = dag(U) @ M.Ztheta @ U
        assert np.isclose(M.g @ Zt @ M.Ztheta @ M.g,
                          M.u**2 * np.exp(-2j * M.j * t) + M.v**2 * np.exp(-1j * M.j * t))
    full = PairModel(theta=np.pi / 8)
    assert np.isclose(abs(full.dvec @ (expm(-1j * full.H * np.pi) @ full.b)) ** 2, 1.0)


def test_generated_algebras_are_full_matrix_algebras():
    """C*(I, H_Psi, Z_0) = M_2 and C*(I, H_Phi, Z_theta) = M_3, so every nonzero
    vector is cyclic; the grounds are selected by the Hamiltonians."""
    HPsi = np.array([[0.0, -M.j], [-M.j, 0.0]])       # restriction of H to span{C, D}
    Z0 = np.diag([1.0, -1.0])
    assert generated_algebra_dim([HPsi, Z0]) == 4
    assert generated_algebra_dim([M.on_carrier(M.H), M.on_carrier(M.Ztheta)]) == 9


def test_mixer_identities_and_rotation():
    """T^2 = 2 P_c^* P_d, Xg = m, Yg = 0, Zg = e; U = exp(-2i theta Y) conjugates
    c_i^*, d_i^* to f_i^*, h_i^*, fixes the vacuum and g, and restricts to the preparation map."""
    assert np.allclose(M.T @ M.T, 2 * dag(M.Pc) @ M.Pd)
    assert np.allclose(M.X @ M.g, M.m)
    assert np.allclose(M.Y @ M.g, 0)
    assert np.allclose(M.Z @ M.g, M.e)
    U = expm(-2j * M.theta * M.Y)
    for a, f in zip(M.c, M.fdag):
        assert np.allclose(U @ dag(a) @ dag(U), f)
    for a, h in zip(M.d, M.hdag):
        assert np.allclose(U @ dag(a) @ dag(U), h)
    assert np.allclose(U @ M.vac, M.vac) and np.allclose(U @ M.g, M.g)
    assert np.allclose(U @ M.C, M.A[:, 0]) and np.allclose(U @ M.D, M.A[:, 1])


def test_endpoint_angles_have_no_defect():
    """At theta = 0 or pi/4 one excitation line vanishes: the minimal output is
    two-dimensional and the channel is onto it."""
    for theta, survivor in ((0.0, "e"), (np.pi / 4, "m")):
        E = PairModel(theta=theta)
        H2, A2 = E.F2.conj().T @ E.H @ E.F2, E.F2.conj().T @ E.A
        Hmin = minimal_space(H2, A2)
        assert Hmin.shape[1] == 2 == rank(A2)
        assert np.isclose(np.linalg.norm(E.b - getattr(E, survivor)), 0, atol=1e-12)


def test_doubled_preparation_channel():
    """L(x, y) = A_0 x + A_{pi/4} y has range span{g, e, m} (no final defect) and
    kernel C (g_Psi, -g_Psi)."""
    L = M.Phi.T @ np.hstack([M.preparation(0.0), M.preparation(np.pi / 4)])
    assert rank(L) == 3
    kernel = null_space(L)
    assert kernel.shape[1] == 1
    expected = np.concatenate([M.gPsi, -M.gPsi]) / np.sqrt(2)
    assert np.isclose(abs(kernel[:, 0] @ expected), 1)
