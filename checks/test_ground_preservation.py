"""Checks for the clauses added to Sections 4 to 7 on returns and ground preservation:
the output ground as an eigenvector of AA*, zero modes of the returns, the differential
with no second step, operators and unitary groups that keep the ground line, and the
independence of the vanishing of the differential from any ground-fixing channel."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm, null_space

from lib.fock import dag
from lib.linalg import rank
from pair.model import PairModel

LABELS = {
    "LG": ["prop:strong-polar-consequences", "prop:vacuum-not-return-zero-mode",
           "eq:shifted-physical-return-model", "prop:no-iterated-awareness-tower",
           "id:rangrig-terminal-law", "def:ground-line-preserving-operator",
           "prop:ground-preservation-reduction", "thm:unitary-ground-preservation",
           "prop:terminality-nonalteration-distinct", "def:momentary-phenomenal-system",
           "rem:annihilation-not-destruction"],
    "UM": ["prop:ground-polar-compatibility", "eq:ground-output-return", "eq:index-cohomology",
           "prop:index-harmonic-pairing", "eq:index-no-second-step", "eq:return-kernels",
           "eq:determination-variance", "eq:ground-line-kept", "eq:ground-spectral-projection",
           "subsec:relaxation-unitary", "eq:relaxation-unitary", "eq:relaxation-unitary-group",
           "subsec:relaxation-moments",
           "eq:index-channel-coherence", "eq:graph-evolution-phase", "part:problem-origin"],
}

RNG = np.random.default_rng(6157)


def random_complex(*shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def unit(vector):
    return vector / np.linalg.norm(vector)


def projection(vector):
    return np.outer(vector, vector.conj())


def map_with_kernel(n_in=4, n_out=6, singular=(0.5, 1.3, 2.0)):
    left = np.linalg.qr(random_complex(n_out, n_out))[0]
    right = np.linalg.qr(random_complex(n_in, n_in))[0]
    r = len(singular)
    return left[:, :r] @ np.diag(singular) @ dag(right[:, :r]), left, right


def test_output_ground_is_an_eigenvector_of_the_output_return():
    A, left, right = map_with_kernel()
    # strong compatibility: x is a right singular vector with value a
    x, a = right[:, 1], 1.3
    y = A @ x / a
    assert_allclose(dag(A) @ A @ x, a ** 2 * x, atol=1e-12)
    assert_allclose(dag(A) @ y, a * x, atol=1e-12)
    assert_allclose(A @ dag(A) @ y, a ** 2 * y, atol=1e-12)
    # the unequal-cell path has only the weak bridge, and y is not an eigenvector of AA*
    A = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]]) / 2
    x, y, a = np.ones(2) / np.sqrt(2), np.ones(5) / np.sqrt(5), np.sqrt(5 / 8)
    assert_allclose(A @ x, a * y, atol=1e-12)
    out = A @ A.T @ y
    assert np.linalg.norm(out - np.vdot(y, out) * y) > 0.05


def test_zero_modes_of_the_returns_and_the_returns_as_hamiltonians():
    A, left, right = map_with_kernel()
    R_in, R_out = dag(A) @ A, A @ dag(A)
    for _ in range(5):
        x = random_complex(4)
        assert_allclose(np.vdot(x, R_in @ x).real, np.linalg.norm(A @ x) ** 2, atol=1e-10)
    assert_allclose(null_space(R_in) @ dag(null_space(R_in)), null_space(A) @ dag(null_space(A)), atol=1e-10)
    assert_allclose(null_space(R_out) @ dag(null_space(R_out)), null_space(dag(A)) @ dag(null_space(dag(A))), atol=1e-10)
    # strong compatibility with a = 1.3, which is not the bottom of the returns
    x, a = right[:, 1], 1.3
    y = A @ x / a
    shifted_in = R_in - a ** 2 * np.eye(4)
    shifted_out = R_out - a ** 2 * np.eye(6)
    assert np.linalg.eigvalsh(shifted_in).min() < -1          # the kernel of A lies below
    assert np.linalg.eigvalsh(shifted_out).min() < -1         # the defect lies below
    defect = null_space(dag(A))
    assert defect.shape[1] == 3
    assert_allclose(dag(defect) @ shifted_out @ defect, -a ** 2 * np.eye(3), atol=1e-10)
    assert_allclose(np.vdot(y, R_out @ y).real, a ** 2, atol=1e-12)

    # LG's shifted return acts on the active input S. A = diag(0, 1, 2), x = y = e_1, a = 1:
    # strongly compatible; on the whole input the kernel lies below x, on S x is the ground
    D3 = np.diag([0.0, 1.0, 2.0])
    e1 = np.eye(3)[1]
    assert_allclose(D3 @ e1, e1, atol=1e-12)
    assert_allclose(D3.T @ D3 @ e1, e1, atol=1e-12)
    assert_allclose(np.linalg.eigvalsh(D3.T @ D3 - np.eye(3)), [-1, 0, 3], atol=1e-12)
    S_basis = np.eye(3)[:, 1:]                                 # (ker A)^perp
    values, vectors = np.linalg.eigh(S_basis.T @ (D3.T @ D3 - np.eye(3)) @ S_basis)
    assert_allclose(values, [0, 3], atol=1e-12)
    assert_allclose(abs(vectors[:, 0] @ (S_basis.T @ e1)), 1, atol=1e-12)
    # A weak bridge can give the same value of a² at a vector outside S.
    weak_input = np.array([np.sqrt(3) / 2, 0, 0.5])
    assert_allclose(D3 @ weak_input, [0, 0, 1], atol=1e-12)
    assert_allclose(np.linalg.norm(S_basis.T @ weak_input) ** 2, 0.25, atol=1e-12)
    assert not np.allclose(D3.T @ D3 @ weak_input, weak_input)
    # without strong compatibility x is not an eigenvector of A*A (unequal-cell path)
    Cmap = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]]) / 2
    xc = np.ones(2) / np.sqrt(2)
    image = Cmap.T @ Cmap @ xc
    assert np.linalg.norm(image - np.vdot(xc, image) * xc) > 0.1
    # Output side with R3 and R5 only (weak bridge, defect): y is orthogonal to the whole
    # ground eigenspace ker A* of AA*.
    Bw = random_complex(5, 3)
    xw = unit(random_complex(3))
    amplitude = np.linalg.norm(Bw @ xw)
    yw = Bw @ xw / amplitude
    vals, vecs = np.linalg.eigh(Bw @ dag(Bw))
    ground_space = vecs[:, vals < 1e-10]
    assert ground_space.shape[1] == 2
    assert_allclose(dag(ground_space) @ yw, 0, atol=1e-10)
    shifted_output = Bw @ dag(Bw) - amplitude ** 2 * np.eye(5)
    assert_allclose(shifted_output @ ground_space, -amplitude ** 2 * ground_space, atol=1e-10)
    out = Bw @ dag(Bw) @ yw
    assert np.linalg.norm(out - np.vdot(yw, out) * yw) > 1e-3   # not even an eigenvector here

    # an injective map with a at the bottom of A*A: x is a ground of A*A - a^2
    left2 = np.linalg.qr(random_complex(5, 5))[0]
    right2 = np.linalg.qr(random_complex(3, 3))[0]
    B = left2[:, :3] @ np.diag([0.7, 1.1, 1.9]) @ dag(right2)
    x = right2[:, 0]
    values, vectors = np.linalg.eigh(dag(B) @ B - 0.49 * np.eye(3))
    assert_allclose(values[0], 0, atol=1e-12)
    assert_allclose(abs(np.vdot(vectors[:, 0], x)), 1, atol=1e-12)

    # the unequal-cell path: A*A = diag(1/4, 1) has the first coarse site as its ground,
    # while the two-site path Hamiltonian has the constant ground x
    C = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]]) / 2
    assert_allclose(C.T @ C, np.diag([0.25, 1.0]), atol=1e-12)
    H2 = np.array([[0.5, -0.5], [-0.5, 0.5]])
    assert_allclose(H2 @ (np.ones(2) / np.sqrt(2)), 0, atol=1e-12)

    # the pair preparation is an isometry: A*A = I, and AA* projects onto D with the defect
    # of the six-dimensional two-particle space as its zero eigenspace
    model = PairModel(theta=0.3, j=1.0)
    prep = model.F2.conj().T @ model.A                           # on the two-particle space
    assert_allclose(dag(prep) @ prep, np.eye(2), atol=1e-12)
    out = prep @ dag(prep)
    assert_allclose(np.sort(np.linalg.eigvalsh(out)), [0, 0, 0, 0, 1, 1], atol=1e-12)


def test_the_differential_has_no_second_step():
    A, _, _ = map_with_kernel()
    n_in, n_out = A.shape[1], A.shape[0]
    d = np.block([[np.zeros((n_in, n_in)), np.zeros((n_in, n_out))], [A, np.zeros((n_out, n_out))]])
    assert_allclose(d @ d, 0, atol=1e-12)
    assert_allclose((-d) @ (-d), 0, atol=1e-12)
    for _ in range(5):
        x, y = random_complex(n_in), random_complex(n_out)
        image = d @ np.concatenate((x, y))
        assert_allclose(image[:n_in], 0, atol=1e-12)
        assert_allclose(image[n_in:], A @ x, atol=1e-12)
        assert_allclose(d @ np.concatenate((np.zeros(n_in), y)), 0, atol=1e-12)
    # the Laplacian is diag(A*A, AA*)
    assert_allclose(dag(d) @ d + d @ dag(d), np.block([[dag(A) @ A, np.zeros((n_in, n_out))],
                                                       [np.zeros((n_out, n_in)), A @ dag(A)]]), atol=1e-12)


def shift_channel(N, a, p):
    """Kraus operators of Section 7 on truncated copies: inputs f_0..f_N, outputs g_0..g_{N+1}."""
    shift = np.zeros((N + 2, N + 1))
    for n in range(N + 1):
        shift[n + 1, n] = 1
    x, f = np.eye(N + 1)[0], np.eye(N + 1)[1]
    z = np.eye(N + 2)[0]
    factor = np.eye(N + 1) - a * projection(x) - p * projection(f)
    L0 = shift @ np.diag(np.sqrt(np.diag(factor).real))
    Lx = np.sqrt(a) * np.outer(z, x)
    Lf = np.sqrt(p) * np.outer(z, f)
    return shift, [L0, Lx, Lf], x, z


def test_vanishing_of_the_differential_does_not_constrain_ground_dynamics():
    N = 6
    for a in (0.0, 0.4):
        shift, kraus, x, z = shift_channel(N, a, 0.3)
        assert_allclose(sum(dag(L) @ L for L in kraus), np.eye(N + 1), atol=1e-12)
        out = sum(L @ projection(x) @ dag(L) for L in kraus)
        g1 = np.eye(N + 2)[1]
        assert_allclose(np.vdot(g1, out @ g1).real, 1 - a, atol=1e-12)
        assert_allclose(np.vdot(z, out @ z).real, a, atol=1e-12)
        # the same preparation, and its differential vanishes on the output summand either way
        d = np.block([[np.zeros((N + 1, N + 1)), np.zeros((N + 1, N + 2))],
                      [shift, np.zeros((N + 2, N + 2))]])
        assert_allclose(d[:, N + 1:], 0, atol=1e-12)
    # a = 0 carries the input ground to the output ground g_1
    shift, kraus, x, z = shift_channel(N, 0.0, 0.3)
    out = sum(L @ projection(x) @ dag(L) for L in kraus)
    assert_allclose(out, projection(np.eye(N + 2)[1]), atol=1e-12)


def test_operators_that_keep_the_ground_line():
    n = 5
    y = unit(random_complex(n))
    P, Q = projection(y), np.eye(n) - projection(y)
    omega = lambda X: np.vdot(y, X @ y)

    def keeps(X):
        return np.linalg.norm(Q @ X @ P) < 1e-10

    for c in (1.7 - 0.3j, 0.0):
        X = c * P + random_complex(n, n) @ Q                    # X y = c y
        assert keeps(X)
        assert_allclose(X @ y, c * y, atol=1e-12)
        assert_allclose(omega(dag(X) @ X) - abs(omega(X)) ** 2, 0, atol=1e-10)
        assert_allclose(omega(X), c, atol=1e-12)
    # one-sided: X = |y><w| with w orthogonal to y keeps the line (X y = 0), X* does not
    w = unit(Q @ random_complex(n))
    X = np.outer(y, w.conj())
    assert keeps(X) and not keeps(dag(X))
    assert np.linalg.norm(X @ P - P @ X) > 0.5
    # projective preservation with a one-sided operator: P + |y><w|
    X = P + np.outer(y, w.conj())
    assert keeps(X) and abs(omega(X)) > 0.5 and not keeps(dag(X))
    # two-sided: block diagonal with respect to P, equivalently commuting with P
    X = 0.4j * P + Q @ random_complex(n, n) @ Q
    assert keeps(X) and keeps(dag(X))
    assert_allclose(X @ P, P @ X, atol=1e-12)
    # a unitary that keeps the line: its adjoint does too
    U = expm(-1j * (0.9 * P + Q @ (lambda M: M + dag(M))(random_complex(n, n)) @ Q))
    assert keeps(U) and keeps(dag(U))
    # a self-adjoint operator keeps the line exactly when y is an eigenvector
    G = (lambda M: M + dag(M))(random_complex(n, n))
    assert not keeps(G)


def test_unitary_groups_that_keep_the_ground():
    n = 4
    y = unit(random_complex(n))
    P, Q = projection(y), np.eye(n) - projection(y)
    H = (lambda M: M + dag(M))(random_complex(n, n))
    G = 0.8 * P + Q @ H @ Q                                    # G y = 0.8 y
    for t in (0.3, 1.7, -2.2):
        U = expm(-1j * t * G)
        assert_allclose(U @ y, np.exp(-0.8j * t) * y, atol=1e-12)
        assert_allclose(U @ P @ dag(U), P, atol=1e-12)
    # if y is not an eigenvector, the projection moves at some times
    moved = [np.linalg.norm(expm(-1j * t * H) @ P @ expm(1j * t * H) - P) for t in (0.1, 0.5)]
    assert max(moved) > 1e-3
    # return times are those with |<y, U_t y>| = 1; for G = diag(0, 1), y = (1, 1)/sqrt 2
    # they are the multiples of 2 pi, a closed subgroup that is not all of R
    G2 = np.diag([0.0, 1.0])
    y2 = np.ones(2) / np.sqrt(2)
    P2 = projection(y2)
    returns = lambda G, v, t: abs(np.vdot(v, expm(-1j * t * G) @ v))
    for t in (2 * np.pi, 4 * np.pi, -6 * np.pi):
        assert_allclose(expm(-1j * t * G2) @ P2 @ expm(1j * t * G2), P2, atol=1e-12)
        assert_allclose(returns(G2, y2, t), 1, atol=1e-12)
    for t in (1.0, np.pi, 2 * np.pi * np.sqrt(2)):
        assert returns(G2, y2, t) < 1 - 1e-3
    # with frequencies 1 and sqrt 2 the time 2 pi is no longer a return
    G3 = np.diag([0.0, 1.0, np.sqrt(2)])
    y3 = np.ones(3) / np.sqrt(3)
    assert returns(G3, y3, 2 * np.pi) < 1 - 1e-3
    # the phase is a character: chi(s + t) = chi(s) chi(t), recovered from the derivative
    chi = lambda t: np.vdot(y, expm(-1j * t * G) @ y)
    assert_allclose(chi(0.4 + 1.1), chi(0.4) * chi(1.1), atol=1e-12)
    h = 1e-6
    assert_allclose((chi(h) - chi(-h)) / (2 * h), -0.8j, atol=1e-8)


def test_each_moment_has_a_rank_one_ground_projection():
    # A vector chosen in a degenerate ground space has a rank-one projection, which is not
    # the spectral projection at the bottom; for a unique bottom the two coincide.
    for diagonal, degenerate in (([0.0, 0.0, 1.0], True), ([0.0, 0.5, 1.0], False)):
        H = np.diag(diagonal)
        values, vectors = np.linalg.eigh(H)
        bottom = vectors[:, np.abs(values - values[0]) < 1e-12]
        spectral = bottom @ dag(bottom)
        chosen = unit(bottom @ np.ones(bottom.shape[1]))
        assert_allclose(H @ chosen, values[0] * chosen, atol=1e-12)
        assert rank(projection(chosen)) == 1
        assert rank(spectral) == (2 if degenerate else 1)
        assert np.allclose(projection(chosen), spectral) == (not degenerate)
