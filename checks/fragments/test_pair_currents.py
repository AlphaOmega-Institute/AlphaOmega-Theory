"""The one-body currents of the pair model: the carrier as their spin-one representation,
the drives and mixers that keep the ground, the readout as a rotated current, and the same
currents in the orbital trap (angular momentum and the two quadrupoles).

The statements of LG that this block carries are checked in pair/test_pair_controls.py,
pair/test_pair_realization.py and the orbital checks.  This file checks what the text adds
to them and the forms in which it states them.  hbar = 1.
"""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm, null_space

from lib.fock import dag, number
from lib.linalg import generated_algebra_dim, minimal_space, rank
from orbital.model import build, first_order, lowest
from pair.model import PairModel

LABELS = {
    "LG": ["eq:pair-parent", "eq:pair-channel", "prop:pair-ground-defect", "eq:pair-readout-action",
           "eq:pair-defect-measurement", "eq:pair-two-frequency-response", "eq:pair-mixer-selection",
           "eq:pair-completed-response", "eq:pair-general-mixer", "eq:pair-mixer-ground-certificate",
           "prop:pair-invariant-mixer", "eq:orbital-parent", "eq:orbital-bias-parameters",
           "eq:orbital-minimal-coupling", "eq:orbital-derived-current", "eq:orbital-quadrupoles",
           "eq:orbital-low-carrier", "eq:orbital-typed-preparation", "eq:orbital-filtered-readout",
           "eq:orbital-low-response", "eq:orbital-unfiltered-response"],
    "UM": ["subsec:pair-grounds-excitations", "eq:pair-ambient-zero", "eq:pair-invariance-amplitude",
           "eq:pair-detector-pullback", "subsec:symmetry-pair", "eq:symmetry-pair-rotation",
           "subsec:orbital-apparatus", "eq:compensated-orbital-parent", "eq:rotated-orbital-readout"],
}

M = PairModel(theta=0.3, j=1.0)
MODES = M.c + M.d                                   # c1, c2, d1, d2


def _commutator(a, b):
    return a @ b - b @ a


def _projection(v):
    return np.outer(v, v.conj())


def _one_body(h):
    """Second quantization of a 4 x 4 matrix in the mode order c1, c2, d1, d2."""
    return sum(h[a, b] * dag(MODES[a]) @ MODES[b] for a in range(4) for b in range(4))


def _mixer(K):
    TK = sum(K[i, k] * dag(M.c[i]) @ M.d[k] for i in range(2) for k in range(2))
    return (TK - dag(TK)) / 2j


def test_currents_act_as_spin_one_on_the_carrier():
    X, Y, Z, T = M.X, M.Y, M.Z, M.T
    assert_allclose(_commutator(X, Y), 1j * Z, atol=1e-12)
    assert_allclose(_commutator(Y, Z), 1j * X, atol=1e-12)
    assert_allclose(_commutator(Z, X), 1j * Y, atol=1e-12)
    assert_allclose(T @ T, 2 * dag(M.Pc) @ M.Pd, atol=1e-12)           # on the whole Fock space
    g, e, m = M.g, M.e, M.m
    table = {"X": (X, m, 0 * g, g), "Y": (Y, 0 * g, 1j * m, -1j * e), "Z": (Z, e, g, 0 * g)}
    for op, on_g, on_e, on_m in table.values():
        assert_allclose(op @ g, on_g, atol=1e-12)
        assert_allclose(op @ e, on_e, atol=1e-12)
        assert_allclose(op @ m, on_m, atol=1e-12)
        assert_allclose(op @ M.W, 0, atol=1e-12)                       # three singlets
    casimir = X @ X + Y @ Y + Z @ Z
    assert_allclose(casimir @ M.Phi, 2 * M.Phi, atol=1e-12)            # spin one: 1 (1 + 1)
    assert_allclose(casimir @ M.W, 0, atol=1e-12)
    # the pair transfer and the symmetric family written in the currents
    assert_allclose(dag(M.Pc) @ M.Pd + dag(M.Pd) @ M.Pc, X @ X - Y @ Y, atol=1e-12)       # on the Fock space
    for s in (0.0, 0.4, 1.0, 1.7):
        H = M.exchange_hamiltonian(s)
        assert_allclose(H @ M.P2, (-s * Z @ Z - (X @ X - Y @ Y)) @ M.P2, atol=1e-12)
        assert_allclose(_commutator(Y, H) @ M.P2, 1j * (1 - s) * (X @ Z + Z @ X) @ M.P2, atol=1e-12)
    assert_allclose((X @ Z + Z @ X) @ e, m, atol=1e-12)
    assert_allclose(_commutator(Y, dag(M.c[0])), 0.5j * dag(M.d[0]), atol=1e-12)
    assert_allclose(_commutator(Y, dag(M.d[0])), -0.5j * dag(M.c[0]), atol=1e-12)
    assert_allclose(Y @ M.vac, 0, atol=1e-12)
    parity = np.diag(np.where(np.round(np.diag(number(M.c[0]) + number(M.c[1]))) % 2 == 0, 1.0, -1.0))
    for current in (X, Y):
        assert_allclose(current @ parity, -parity @ current, atol=1e-12)
    # each of g, e, m is the zero vector of one current, and Z has C, D, m as eigenvectors
    assert_allclose(Z @ M.C, M.C, atol=1e-12)
    assert_allclose(Z @ M.D, -M.D, atol=1e-12)


def test_rotation_generator_and_the_readout_as_a_rotated_current():
    for phi in (0.3, 1.1, -0.7):
        kappa = expm(-2j * phi * M.Y)
        f, h = M.rotated_creators(phi)
        for c, fdag in zip(M.c, f):
            assert_allclose(kappa @ dag(c) @ dag(kappa), fdag, atol=1e-12)
        for d, hdag in zip(M.d, h):
            assert_allclose(kappa @ dag(d) @ dag(kappa), hdag, atol=1e-12)
        assert_allclose(kappa @ M.vac, M.vac, atol=1e-12)
        assert_allclose(kappa @ M.g, M.g, atol=1e-12)
        assert_allclose(kappa @ M.preparation(0.0), M.preparation(phi), atol=1e-12)
        u, v = np.cos(2 * phi), np.sin(2 * phi)
        assert_allclose(kappa @ M.e, u * M.e + v * M.m, atol=1e-12)
        assert_allclose(kappa @ M.Z @ dag(kappa), u * M.Z + v * M.X, atol=1e-12)
    # the measured imbalance at the model's angle
    assert_allclose(M.Ztheta, M.u * M.Z + M.v * M.X, atol=1e-12)
    assert_allclose(M.Ztheta @ M.g, M.b, atol=1e-12)
    assert_allclose(M.Ztheta @ M.b, M.g, atol=1e-12)
    assert_allclose(M.Ztheta @ M.dvec, 0, atol=1e-12)
    square = M.Ztheta @ M.Ztheta
    assert_allclose(M.on_carrier(square), np.eye(3) - M.on_carrier(_projection(M.dvec)), atol=1e-12)
    assert_allclose(M.P2 - square @ M.P2, _projection(M.dvec) + M.W @ dag(M.W), atol=1e-12)
    assert_allclose(square @ M.P2, M.A @ dag(M.A), atol=1e-12)


def test_common_label_drives_that_keep_the_ground():
    # Hermitian one-body operators acting identically on both labels: span of N, X, Y, Z
    samples = []
    rng = np.random.default_rng(41)
    for _ in range(6):
        a, b = rng.normal(size=2)
        t = rng.normal() + 1j * rng.normal()
        block = np.array([[a, t], [np.conj(t), b]])
        drive = _one_body(np.kron(block, np.eye(2)))                   # the same 2 x 2 block for both labels
        assert_allclose(drive, (a + b) / 2 * M.N + (a - b) * M.Z + 2 * t.real * M.X - 2 * t.imag * M.Y, atol=1e-12)
        samples.append(drive.ravel())
    assert rank(np.array(samples).T) == 4
    # the ground stays on its ray for every pulse length exactly when the X and Z parts vanish
    for _ in range(40):
        alpha, beta, gamma = rng.normal(size=3)
        if rng.random() < 0.5:
            alpha = gamma = 0.0
        drive = alpha * M.X + beta * M.Y + gamma * M.Z
        leaves = any(abs(np.vdot(M.g, expm(-1j * t * drive) @ M.g)) < 1 - 1e-9 for t in (0.3, 1.0, 2.2))
        assert leaves == (alpha != 0 or gamma != 0)
        assert_allclose(np.linalg.norm(drive @ M.g), np.hypot(alpha, gamma), atol=1e-12)
    # a pulse with a time-dependent Y drive only
    state = M.g.astype(complex)
    for amplitude in rng.normal(size=20):
        state = expm(-0.1j * amplitude * M.Y) @ state
    assert_allclose(state, M.g, atol=1e-12)


def test_one_body_unitaries_that_keep_the_ground_ray():
    # a real basis of the Hermitian 4 x 4 matrices
    hermitian = []
    for a in range(4):
        for b in range(a, 4):
            E = np.zeros((4, 4), dtype=complex)
            if a == b:
                E[a, a] = 1
                hermitian.append(E)
            else:
                E[a, b] = E[b, a] = 1
                hermitian.append(E)
                F = np.zeros((4, 4), dtype=complex)
                F[a, b], F[b, a] = 1j, -1j
                hermitian.append(F)
    assert len(hermitian) == 16
    columns = []
    for h in hermitian:
        image = _one_body(h) @ M.g
        image = image - M.g * np.vdot(M.g, image)
        columns.append(np.concatenate([image.real, image.imag]))
    conditions = np.array(columns).T
    assert 16 - rank(conditions) == 11
    # they are the identity and the Hermitian matrices of the compact symplectic algebra of eps (+) eps
    eps = np.array([[0, 1], [-1, 0]])
    form = np.kron(np.eye(2), eps)
    keepers = null_space(conditions)
    for column in keepers.T:
        h = sum(x * E for x, E in zip(column, hermitian))
        h = h - np.trace(h) / 4 * np.eye(4)
        assert_allclose(h @ form + form @ h.T, 0, atol=1e-9)
    # a finite unitary of that group keeps the ray, and its off-diagonal generators are the quaternionic mixers
    z, w = 0.3 + 0.4j, -0.2 + 0.7j
    K = np.array([[z, w], [-np.conj(w), np.conj(z)]])
    assert_allclose(_mixer(K) @ M.g, 0, atol=1e-12)
    assert_allclose(expm(-0.8j * _mixer(K)) @ M.g, M.g, atol=1e-12)
    internal = np.array([[0.2, 0.5 - 0.1j], [0.5 + 0.1j, -0.2]])      # traceless: an internal rotation
    for block in (np.kron(np.diag([1, 0]), internal), np.kron(np.diag([0, 1]), internal)):
        assert_allclose(_one_body(block) @ M.g, 0, atol=1e-12)
    # counting: three internal rotations on each orbital, four mixers, and the particle number
    assert 3 + 3 + 4 + 1 == 11


def test_a_mixer_can_keep_the_ground_and_leave_the_detector_space():
    K = np.diag([1j, -1j])                                        # of the form that annihilates the ground
    YK = _mixer(K)
    opposite = sum(sign * (dag(c) @ d + dag(d) @ c) / 2 for sign, c, d in zip((1, -1), M.c, M.d))
    assert_allclose(YK, opposite, atol=1e-12)                    # the current X with opposite signs on the labels
    assert_allclose(YK @ M.g, 0, atol=1e-12)
    assert_allclose(YK @ M.e, -M.w0, atol=1e-12)
    assert_allclose(YK @ M.w0, -M.e, atol=1e-12)
    for t in (0.3, 1.0, np.pi / 2):
        pulse = expm(-1j * t * YK)
        assert_allclose(pulse @ M.g, M.g, atol=1e-12)
        assert_allclose(pulse @ M.e, np.cos(t) * M.e + 1j * np.sin(t) * M.w0, atol=1e-12)
    for theta in (0.2, 0.3, np.pi / 8):
        model = PairModel(theta=theta, j=1.0)
        out = expm(-0.5j * np.pi * YK) @ model.e                 # entirely in W
        assert_allclose(np.linalg.norm(dag(model.W) @ out), 1.0, atol=1e-12)
        assert_allclose(np.vdot(model.dvec, out), 0, atol=1e-12)
        assert_allclose(1 - np.vdot(out, model.Ztheta @ model.Ztheta @ out).real, 1.0, atol=1e-12)   # zero with certainty
    # among the mixers that annihilate the ground, only real multiples of the aligned one keep span{g, e, m}
    family = [np.eye(2, dtype=complex), 1j * np.array([[0, 1], [1, 0]]),
              np.array([[0, 1], [-1, 0]], dtype=complex), 1j * np.diag([1.0, -1.0])]
    outside = M.P2 - M.Phi @ dag(M.Phi)
    columns = []
    for matrix in family:
        assert_allclose(_mixer(matrix) @ M.g, 0, atol=1e-12)
        leak = outside @ _mixer(matrix) @ M.Phi
        columns.append(np.concatenate([leak.real.ravel(), leak.imag.ravel()]))
    assert 4 - rank(np.array(columns).T) == 1
    assert_allclose(outside @ _mixer(family[0]) @ M.Phi, 0, atol=1e-12)
    # from the ground, the common-label currents and H_s generate exactly span{g, e, m}
    for s in (0.0, 0.4):
        generators = [dag(M.F2) @ op @ M.F2 for op in (M.X, M.Y, M.Z, M.exchange_hamiltonian(s))]
        orbit = [dag(M.F2) @ M.g]
        for _ in range(3):
            orbit = orbit + [op @ vector for op in generators for vector in orbit]
        assert rank(np.array(orbit).T) == 3                       # not cyclic on the six-dimensional space


def test_mixer_expansion_in_the_four_cross_states():
    rng = np.random.default_rng(42)
    cross = {(i, k): dag(M.c[i]) @ dag(M.d[k]) @ M.vac for i in range(2) for k in range(2)}
    for _ in range(50):
        K = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        TK = sum(K[i, k] * dag(M.c[i]) @ M.d[k] for i in range(2) for k in range(2))
        assert_allclose(TK @ M.C, 0, atol=1e-12)
        assert_allclose(dag(TK) @ M.D, 0, atol=1e-12)
        expected = ((K[0, 0] - np.conj(K[1, 1])) * cross[(0, 1)] - (K[1, 1] - np.conj(K[0, 0])) * cross[(1, 0)]
                    - (K[0, 1] + np.conj(K[1, 0])) * cross[(0, 0)] + (K[1, 0] + np.conj(K[0, 1])) * cross[(1, 1)])
        assert_allclose(_mixer(K) @ M.g, expected / (2j * np.sqrt(2)), atol=1e-12)
        leakage = (abs(K[0, 0] - np.conj(K[1, 1]))**2 + abs(K[0, 1] + np.conj(K[1, 0]))**2) / 4
        assert_allclose(np.linalg.norm(_mixer(K) @ M.g)**2, leakage, atol=1e-12)
    # invariance of the mixer under c -> R_c c, d -> R_d d is R_c K = K R_d
    for _ in range(20):
        Rc, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        Rd, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        K = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        transformed = Rc.conj().T @ K @ Rd
        assert np.allclose(transformed, K) == np.allclose(Rc @ K, K @ Rd)
        # the pair term picks up conj(det R_c) det R_d
        new_c = [sum(Rc[i, k] * M.c[k] for k in range(2)) for i in range(2)]
        new_d = [sum(Rd[i, k] * M.d[k] for k in range(2)) for i in range(2)]
        assert_allclose(dag(new_c[1] @ new_c[0]) @ (new_d[1] @ new_d[0]),
                        np.conj(np.linalg.det(Rc)) * np.linalg.det(Rd) * dag(M.Pc) @ M.Pd, atol=1e-12)


def test_input_generator_is_a_restriction_and_not_the_compression():
    for theta in (0.2, 0.5):
        model = PairModel(theta=theta, j=1.0)
        frame = np.column_stack([model.gPsi, model.ePsi])
        for s in (0.0, 0.4, 1.0, 1.6):
            H = model.exchange_hamiltonian(s)
            pair_plane = model.preparation(0.0) @ frame
            prepared = model.A @ frame
            assert_allclose(dag(pair_plane) @ H @ pair_plane, np.diag([-(s + 1), 1 - s]), atol=1e-12)
            assert_allclose(dag(prepared) @ H @ prepared, np.diag([-(s + 1), (1 - s) * model.u**2]), atol=1e-12)


def test_readout_moments_correlation_and_generated_algebras():
    # the mean of the imbalance does not determine the zero outcome
    plus = (M.g + M.b) / np.sqrt(2)                   # Z_theta = +1
    minus = (M.g - M.b) / np.sqrt(2)                  # Z_theta = -1
    assert_allclose(M.Ztheta @ plus, plus, atol=1e-12)
    assert_allclose(M.Ztheta @ minus, -minus, atol=1e-12)
    mixture = (_projection(plus) + _projection(minus)) / 2
    on_line = _projection(M.dvec)
    square = M.Ztheta @ M.Ztheta
    for state, zero_probability in ((mixture, 0.0), (on_line, 1.0)):
        assert_allclose(np.trace(state @ M.Ztheta), 0, atol=1e-12)
        assert_allclose(1 - np.trace(state @ square).real, zero_probability, atol=1e-12)
    # ground correlation of the readout for the symmetric family
    for s in (0.0, 0.4, 1.0):
        H = M.exchange_hamiltonian(s)
        for t in (0.6, 2.3):
            U = expm(-1j * t * H)
            correlation = np.vdot(M.g, dag(U) @ M.Ztheta @ U @ M.Ztheta @ M.g)
            assert_allclose(correlation, M.u**2 * np.exp(-2j * t) + M.v**2 * np.exp(-1j * (s + 1) * t), atol=1e-12)
    # the excitation alone reaches two dimensions; the preparation, which also prepares the ground, three
    H2 = dag(M.F2) @ M.exchange_hamiltonian(0.0) @ M.F2
    assert minimal_space(H2, (dag(M.F2) @ M.b)[:, None]).shape[1] == 2
    assert minimal_space(H2, dag(M.F2) @ M.A).shape[1] == 3
    # the algebras generated by the Hamiltonian and the imbalance on the two carriers
    frame = np.column_stack([M.gPsi, M.ePsi])
    H_in = np.diag([-1.0, 1.0])
    Z_in = np.array([[0.0, 1.0], [1.0, 0.0]])        # the input imbalance exchanges ground and excitation
    assert_allclose(dag(M.A @ frame) @ M.Ztheta @ (M.A @ frame), Z_in, atol=1e-12)
    assert generated_algebra_dim([H_in, Z_in]) == 4
    assert generated_algebra_dim([M.on_carrier(M.H), M.on_carrier(M.Ztheta)]) == 9
    # at an endpoint angle one line is lost and the generated algebra on the carrier is smaller
    endpoint = PairModel(theta=0.0, j=1.0)
    assert generated_algebra_dim([endpoint.on_carrier(endpoint.H), endpoint.on_carrier(endpoint.Ztheta)]) < 9


class _Gaussian:
    """A polynomial in x, y times exp(-(x^2 + y^2) / 2); coefficients[a, b] multiplies x^a y^b."""

    def __init__(self, coefficients):
        self.c = np.array(coefficients, dtype=complex)

    @staticmethod
    def monomial(a, b, size=10):
        c = np.zeros((size, size), dtype=complex)
        c[a, b] = 1
        return _Gaussian(c)

    def __add__(self, other):
        return _Gaussian(self.c + other.c)

    def __sub__(self, other):
        return _Gaussian(self.c - other.c)

    def scale(self, factor):
        return _Gaussian(factor * self.c)

    def times_x(self):
        out = np.zeros_like(self.c)
        out[1:, :] = self.c[:-1, :]
        return _Gaussian(out)

    def times_y(self):
        out = np.zeros_like(self.c)
        out[:, 1:] = self.c[:, :-1]
        return _Gaussian(out)

    def dx(self):
        out = np.zeros_like(self.c)
        out[:-1, :] = self.c[1:, :] * np.arange(1, self.c.shape[0])[:, None]
        return _Gaussian(out) - self.times_x()

    def dy(self):
        out = np.zeros_like(self.c)
        out[:, :-1] = self.c[:, 1:] * np.arange(1, self.c.shape[1])[None, :]
        return _Gaussian(out) - self.times_y()

    def momentum(self, axis):
        return (self.dx() if axis == "x" else self.dy()).scale(-1j)

    def angular_momentum(self):
        return self.momentum("y").times_x() - self.momentum("x").times_y()


def test_orbital_currents_and_the_minimal_coupling_square():
    s = _Gaussian.monomial(0, 0)
    px = _Gaussian.monomial(1, 0)
    py = _Gaussian.monomial(0, 1)
    # one-particle angular momentum on the three orbitals (oscillator length one)
    assert_allclose(s.angular_momentum().c, 0, atol=1e-14)
    assert_allclose(px.angular_momentum().c, 1j * py.c, atol=1e-14)
    assert_allclose(py.angular_momentum().c, -1j * px.c, atol=1e-14)
    # in second quantization on the valence modes this is 2 Y
    lifted = sum(1j * dag(d) @ c - 1j * dag(c) @ d for c, d in zip(M.c, M.d))
    assert_allclose(lifted, 2 * M.Y, atol=1e-12)
    for phi in (0.4, 1.3):
        assert_allclose(expm(-1j * phi * lifted), expm(-2j * phi * M.Y), atol=1e-12)
    # the spin component S_z annihilates the three singlets and w_0, and is +-1 on the other two triplet states
    spin = (number(M.c[0]) + number(M.d[0]) - number(M.c[1]) - number(M.d[1])) / 2
    for vector in (M.g, M.e, M.m, M.w0):
        assert_allclose(spin @ vector, 0, atol=1e-12)
    assert_allclose(spin @ M.W[:, 0], M.W[:, 0], atol=1e-12)
    assert_allclose(spin @ M.W[:, 1], -M.W[:, 1], atol=1e-12)
    # minimal coupling in the symmetric gauge, with M = 1 and charge times field = 2 Omega
    omega0, Omega = 0.8, 0.35
    for a, b in ((0, 0), (1, 0), (0, 1), (2, 1), (1, 3)):
        f = _Gaussian.monomial(a, b)

        def kinetic_x(h):                            # p_x - q A_x = p_x + Omega y
            return h.momentum("x") + h.times_y().scale(Omega)

        def kinetic_y(h):                            # p_y - q A_y = p_y - Omega x
            return h.momentum("y") - h.times_x().scale(Omega)

        r2 = f.times_x().times_x() + f.times_y().times_y()
        coupled = (kinetic_x(kinetic_x(f)) + kinetic_y(kinetic_y(f))).scale(0.5) + r2.scale(omega0**2 / 2)
        free = (f.momentum("x").momentum("x") + f.momentum("y").momentum("y")).scale(0.5)
        expanded = free + r2.scale((omega0**2 + Omega**2) / 2) - f.angular_momentum().scale(Omega)
        assert_allclose(coupled.c, expanded.c, atol=1e-12)


def test_triplet_gap_with_a_spin_moment_has_either_sign():
    # first-order shell: singlet excitations at lam d_a, triplet at 2 lam d_a, field -Omega (L + g_s S_z), L = 2 Y
    spin = (number(M.c[0]) + number(M.d[0]) - number(M.c[1]) - number(M.d[1])) / 2
    excited = _projection(M.e) + _projection(M.m)
    triplet = M.W @ dag(M.W)
    scale = 0.01                                                   # lam d_a
    for chi, g_s in ((0.5, 2.0), (0.9, 2.0), (0.75, 8.0), (0.5, -10.0), (0.8, 5.0)):
        Omega = chi * scale / 2
        shell = scale * excited + 2 * scale * triplet - Omega * (2 * M.Y + g_s * spin)
        energies, vectors = np.linalg.eigh(dag(M.F2) @ shell @ M.F2)
        states = M.F2 @ vectors
        spins = np.array([np.vdot(v, spin @ v).real for v in states.T])
        ground_energy = np.vdot(M.g, shell @ M.g).real
        assert_allclose(shell @ M.g, ground_energy * M.g, atol=1e-14)             # g stays an eigenvector
        lowest_triplet = energies[np.abs(spins) > 0.5].min()
        gap = lowest_triplet - ground_energy
        assert_allclose(gap, scale * (2 - abs(g_s * chi) / 2), atol=1e-14)
        if np.isclose(abs(g_s * chi), 4):
            assert_allclose(gap, 0, atol=1e-14)                    # the first order does not decide
            continue
        assert (gap > 0) == (abs(g_s * chi) < 4)
        below = energies.min() < ground_energy - 1e-12             # something lies below g
        assert below == (abs(g_s * chi) > 4)
    # in the truncated trap the triplet gap tends to 2 lam d_a, so the shifted gap tends to the formula
    xi, d_a = 1.0, first_order(1.0, 1.0)["d_a"]
    ops = build(3, xi)
    errors = []
    for lam in (0.05, 0.2):
        energies, vectors = lowest((ops["Hho"] - lam * ops["W"]).astype(complex), 4)
        spin_squared = [np.vdot(vectors[:, k], ops["S2"] @ vectors[:, k]).real for k in range(4)]
        assert_allclose(spin_squared, [0, 0, 0, 2], atol=1e-6)                    # ground, two chiral singlets, triplet
        triplet_gap = energies[3] - energies[0]
        for chi, g_s in ((0.5, 2.0), (0.75, 8.0)):
            # the component with spin projection +-1 moves by |g_s| Omega, with 2 Omega = chi lam d_a
            shifted = triplet_gap - abs(g_s * chi) * lam * d_a / 2
            predicted = lam * d_a * (2 - abs(g_s * chi) / 2)
            assert abs(shifted - predicted) < 0.2 * lam * lam * d_a               # the remainder is second order
            assert (shifted > 0) == (abs(g_s * chi) < 4)
        errors.append(abs(triplet_gap / (lam * d_a) - 2))
    assert errors[0] < errors[1]


def test_typed_orbital_preparation_on_three_exact_levels():
    E0, E2, Omega = -1.3, 0.4, 0.17
    # basis G, F_plus, F_minus with energies E0, E2 - 2 Omega, E2 + 2 Omega
    H = np.diag([E0, E2 - 2 * Omega, E2 + 2 * Omega]).astype(complex)
    G = np.array([1, 0, 0], dtype=complex)
    e = np.array([0, 1, 1]) / np.sqrt(2)
    m = np.array([0, 1, -1]) / (1j * np.sqrt(2))
    A = np.column_stack([G, e])
    assert_allclose(dag(A) @ A, np.eye(2), atol=1e-14)
    assert_allclose(dag(A) @ H @ A, np.diag([E0, E2]), atol=1e-14)        # the compression
    assert_allclose((np.eye(3) - A @ dag(A)) @ H @ A[:, 1], -2j * Omega * m, atol=1e-14)
    assert_allclose(null_space(dag(A))[:, 0] * np.vdot(null_space(dag(A))[:, 0], m), m, atol=1e-12)
    Z = np.outer(G, e.conj()) + np.outer(e, G.conj())
    assert_allclose(Z @ Z, A @ dag(A), atol=1e-14)
    for t in (0.4, 3.1):
        U = expm(-1j * t * H)
        assert_allclose(abs(np.vdot(m, U @ e))**2, np.sin(2 * Omega * t)**2, atol=1e-12)
        assert_allclose(np.vdot(G, dag(U) @ Z @ U @ Z @ G), np.exp(-1j * (E2 - E0) * t) * np.cos(2 * Omega * t), atol=1e-12)
    assert generated_algebra_dim([dag(A) @ H @ A, dag(A) @ Z @ A]) == 4
    assert generated_algebra_dim([H, Z]) == 9
    assert generated_algebra_dim([np.diag([E0, E2, E2]).astype(complex), Z]) < 9       # without the field
    assert minimal_space(H, e[:, None]).shape[1] == 2            # the excitation alone: the two chiral states
    assert minimal_space(H, A).shape[1] == 3                     # the preparation also supplies the ground
    # an unfiltered port with several radial levels, paired by time reversal at zero field
    rng = np.random.default_rng(43)
    levels = rng.uniform(0.5, 3.0, size=4)
    weights = rng.dirichlet(np.ones(4))
    energies = np.concatenate([levels - 2 * Omega, levels + 2 * Omega])
    start = np.concatenate([np.sqrt(weights / 2), np.sqrt(weights / 2)])
    for t in (0.7, 2.9):
        zero_field = np.sum(weights * np.exp(-1j * levels * t))
        with_field = np.vdot(start, np.exp(-1j * energies * t) * start)
        assert_allclose(with_field, zero_field * np.cos(2 * Omega * t), atol=1e-12)
        # the evolving state stays orthogonal to the ground, so this is the weight outside the prepared plane
        assert_allclose(1 - abs(with_field)**2, 1 - abs(zero_field)**2 * np.cos(2 * Omega * t)**2, atol=1e-12)
