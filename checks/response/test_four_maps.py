"""Checks for §12, "Four maps with one algebra", in Universal_Model.tex.

Operators on infinite-dimensional spaces are tested through exact finite pieces (the
oscillator map between the first M and the first M + 1 levels), through their Fourier
multipliers with numerical integration, or through closed forms.
"""
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm, null_space
from scipy.special import erfc, gamma, zeta

from lib.linalg import minimal_space, response_weights

LABELS = {
    "TAUM": ["ex:scalar-channel", "ex:split-excitation-response", "ex:translation-channel", "eq:no-translation-gap",
             "eq:free-response", "thm:oscillator-channel", "cor:oscillator-heat-index", "thm:oscillator-zeta",
             "eq:oscillator-response", "ex:quasiperiodic-channel", "eq:unsmoothed-point-mass", "eq:theta-mass",
             "eq:quasiperiodic-smoothed-response", "eq:normalized-smoothed-point"],
    "UM": ["sec:four-maps", "eq:four-maps-scalar", "eq:four-maps-scalar-response",
           "eq:four-maps-three-state-data", "eq:four-maps-three-state-response", "eq:four-maps-three-state-weights",
           "eq:four-maps-three-state-transfer", "eq:four-maps-translation", "eq:four-maps-translation-returns",
           "eq:four-maps-translation-no-gap", "eq:four-maps-translation-response", "eq:four-maps-ladder",
           "prop:four-maps-oscillator", "eq:four-maps-oscillator-returns", "eq:four-maps-oscillator-index",
           "eq:four-maps-oscillator-heat", "prop:four-maps-mellin", "eq:four-maps-mellin",
           "eq:four-maps-oscillator-port", "eq:four-maps-oscillator-response", "eq:four-maps-interval-operator",
           "eq:four-maps-interval-spectrum", "eq:four-maps-interval-point", "eq:four-maps-interval-port",
           "eq:four-maps-interval-mass", "eq:four-maps-interval-response",
           "eq:index-heat", "prop:index-harmonic-pairing", "cor:passive-sources", "prop:passive-opposite-sectors",
           "prop:passive-two-spans", "prop:passive-uniqueness", "prop:passive-heat-discrete", "cor:finite-response-defect"],
}


def _graded(a):
    """d_A, the symmetry J and the Laplacian for a matrix A from C^p to C^q."""
    q, p = a.shape
    d = np.zeros((p + q, p + q), dtype=complex)
    d[p:, :p] = a
    return d, np.diag(np.r_[np.ones(p), -np.ones(q)]), np.block([[a.conj().T @ a, np.zeros((p, q))], [np.zeros((q, p)), a @ a.conj().T]])


def _common_relations(a):
    """The relations that hold for every map; returns its index."""
    d, symmetry, laplacian = _graded(a)
    hilbert, krein = d + d.conj().T, d.conj().T - d
    np.testing.assert_allclose(d @ d, 0, atol=1e-12)
    np.testing.assert_allclose(hilbert @ hilbert, laplacian, atol=1e-9)
    np.testing.assert_allclose(krein @ krein, -laplacian, atol=1e-9)
    np.testing.assert_allclose(symmetry @ hilbert, -hilbert @ symmetry, atol=1e-12)
    np.testing.assert_allclose(symmetry @ krein.conj().T @ symmetry, krein, atol=1e-12)      # self-adjoint for the form
    first = np.sort(np.linalg.eigvalsh(a.conj().T @ a))[::-1]
    second = np.sort(np.linalg.eigvalsh(a @ a.conj().T))[::-1]
    scale = max(1.0, first[0])
    k = int(np.sum(first > 1e-9 * scale))
    assert k == int(np.sum(second > 1e-9 * scale))
    np.testing.assert_allclose(first[:k], second[:k], atol=1e-8 * scale)                       # paired nonzero spectra
    return (a.shape[1] - k) - (a.shape[0] - k)


def test_scalar_map():
    rng = np.random.default_rng(401)
    for b in list(rng.normal(size=5) + 1j * rng.normal(size=5)) + [0.0]:
        a = np.array([[b]], dtype=complex)
        d, _, laplacian = _graded(a)
        np.testing.assert_allclose(d + d.conj().T, [[0, np.conj(b)], [b, 0]], atol=1e-15)
        np.testing.assert_allclose(d.conj().T - d, [[0, np.conj(b)], [-b, 0]], atol=1e-15)
        np.testing.assert_allclose(laplacian, abs(b) ** 2 * np.eye(2), atol=1e-15)
        assert _common_relations(a) == 0
        if b != 0:
            phase = b / abs(b)
            np.testing.assert_allclose(phase * abs(b), b, atol=1e-15)                        # polar decomposition
            for z in (0.3 + 0.2j, -2.0 - 1.0j):
                np.testing.assert_allclose(1 / (abs(b) ** 2 - z), np.linalg.inv(a.conj().T @ a - z * np.eye(1))[0, 0], atol=1e-12)
            weights = response_weights(a.conj().T @ a, np.eye(1))
            assert len(weights) == 1 and np.isclose(weights[0][0], abs(b) ** 2) and np.isclose(weights[0][1][0, 0], 1.0)
        else:
            assert null_space(a).shape[1] == 1 and null_space(a.conj().T).shape[1] == 1


def test_three_state_port_as_a_preparation():
    e0, e2, ell, hbar = 0.2, 1.5, 0.8, 0.7
    g, f_plus, f_minus = np.eye(3, dtype=complex)
    e, d = (f_plus + f_minus) / np.sqrt(2), (f_plus - f_minus) / (1j * np.sqrt(2))
    port = np.column_stack([g, e])
    charge = np.diag([0.0, ell, -ell])
    for nu in (0.4, -0.9, 0.0):
        h = np.diag([e0, e2, e2]) - nu * charge
        assert abs(nu) * ell < e2 - e0
        assert np.argmin(np.diag(h).real) == 0 and np.sum(np.isclose(np.diag(h), e0)) == 1     # g is the unique ground
        np.testing.assert_allclose(h @ charge, charge @ h, atol=1e-15)
        for z in (0.1 + 0.3j, 2.0 - 0.5j):
            expected = np.diag([1 / (e0 - z), 0.5 / (e2 - nu * ell - z) + 0.5 / (e2 + nu * ell - z)])
            np.testing.assert_allclose(port.conj().T @ np.linalg.inv(h - z * np.eye(3)) @ port, expected, atol=1e-12)
        for t in (0.6, 3.1):
            evolved = expm(-1j * t * (h - e0 * np.eye(3)) / hbar) @ e
            angle = nu * ell * t / hbar
            np.testing.assert_allclose(evolved, np.exp(-1j * t * (e2 - e0) / hbar) * (np.cos(angle) * e - np.sin(angle) * d), atol=1e-12)
            np.testing.assert_allclose(abs(np.vdot(d, evolved)) ** 2, np.sin(angle) ** 2, atol=1e-12)
        weights = [w for _, w in response_weights(h, port) if np.linalg.norm(w) > 1e-12]
        count = sum(np.linalg.matrix_rank(w) for w in weights) - np.linalg.matrix_rank(sum(weights))
        assert count == (1 if nu else 0) and minimal_space(h, port).shape[1] == (3 if nu else 2)
        if nu:
            # complete transfer first at t_* = pi hbar / (2 |nu| ell); within a window the weight is bounded by (|nu| ell T / hbar)^2
            weight = lambda t: abs(np.vdot(d, expm(-1j * t * h / hbar) @ e)) ** 2
            t_star = np.pi * hbar / (2 * abs(nu) * ell)
            np.testing.assert_allclose(weight(t_star), 1.0, atol=1e-12)
            assert max(weight(t) for t in np.linspace(0, 0.98 * t_star, 400)) < 1 - 1e-4
            for window in (0.05, 0.4, 3.0):
                largest = max(weight(t) for t in np.linspace(0, window, 2001))
                rate = abs(nu) * ell * window / hbar
                np.testing.assert_allclose(largest, np.sin(min(rate, np.pi / 2)) ** 2, atol=1e-6)
                assert largest <= rate**2 + 1e-12
        reached = minimal_space(h, port)
        np.testing.assert_allclose(np.linalg.norm(reached.conj().T @ d) ** 2, 1.0 if nu else 0.0, atol=1e-10)
        # the port read as a preparation A = Gamma from C^2 to C^3
        np.testing.assert_allclose(port.conj().T @ port, np.eye(2), atol=1e-15)
        defect = null_space(port.conj().T)
        assert defect.shape[1] == 1 and np.isclose(abs(np.vdot(defect[:, 0], d)), 1.0)
        assert _common_relations(port) == -1
        np.testing.assert_allclose(port.conj().T @ h @ port, np.diag([e0, e2]), atol=1e-15)
        np.testing.assert_allclose(port @ np.array([1.0, 0.0]), g, atol=1e-15)                 # the input ground maps to g
        assert not np.allclose(h, port @ port.conj().T)                                        # H is not A A^*
        leaves = np.linalg.norm((np.eye(3) - port @ port.conj().T) @ h @ port)
        assert (leaves > 1e-3) == (nu != 0)                                                    # the plane is invariant only at zero bias
    # beyond the inequality the ground is no longer g
    h = np.diag([e0, e2, e2]) - 2.0 * charge
    assert np.argmin(np.diag(h).real) != 0


def _packet(sigma, k0, x0):
    """Fourier transform of a unit Gaussian wave packet."""
    return lambda k: (2 * sigma**2 / np.pi) ** 0.25 * np.exp(-sigma**2 * (k - k0) ** 2) * np.exp(-1j * k * x0)


def _complex_quad(function, a, b, **options):
    real = quad(lambda x: function(x).real, a, b, limit=400, **options)[0]
    imag = quad(lambda x: function(x).imag, a, b, limit=400, **options)[0]
    return real + 1j * imag


def test_translation_generator():
    # low-frequency unit vectors: no positive lower bound on the support
    for n in (1, 5, 40, 1000):
        height = np.sqrt(n / 2.0)                                # Fourier transform equal to this constant on [-1/n, 1/n]
        norm = quad(lambda k: height**2, -1 / n, 1 / n)[0]
        derivative = np.sqrt(quad(lambda k: k**2 * height**2, -1 / n, 1 / n)[0])
        np.testing.assert_allclose(norm, 1.0, atol=1e-12)
        np.testing.assert_allclose(derivative, 1 / (np.sqrt(3) * n), rtol=1e-10)
        assert derivative <= 1 / n
    # the response of a centred Gaussian against a closed form
    sigma = 0.8
    g_hat = _packet(sigma, 0.0, 0.0)
    for z in (-0.7, -3.0, 0.5 + 0.4j, 1.2 - 0.3j, -0.2 + 2.0j):
        beta = np.sqrt(-z + 0j)
        beta = beta if beta.real > 0 else -beta
        closed = np.sqrt(2 * sigma**2 / np.pi) * (np.pi / beta) * np.exp(2 * sigma**2 * beta**2) * erfc(beta * np.sqrt(2) * sigma)
        numeric = _complex_quad(lambda k: abs(g_hat(k)) ** 2 / (k**2 - z), -np.inf, np.inf)
        np.testing.assert_allclose(numeric, closed, atol=1e-8)
        assert (numeric.imag > 0) == (np.imag(z) > 0) or np.imag(z) == 0
    # the same response from the measure on the spectrum of the return, for a packet that moves to the right
    g_hat = _packet(0.9, 1.3, 0.4)
    density = lambda lam: (abs(g_hat(np.sqrt(lam))) ** 2 + abs(g_hat(-np.sqrt(lam))) ** 2) / (2 * np.sqrt(lam))
    np.testing.assert_allclose(quad(density, 0, np.inf, limit=400)[0], 1.0, atol=1e-7)
    for z in (-0.6, 1.0 + 0.5j, 2.5 - 1.0j):
        in_k = _complex_quad(lambda k: abs(g_hat(k)) ** 2 / (k**2 - z), -np.inf, np.inf)
        in_lambda = _complex_quad(lambda lam: density(lam) / (lam - z), 0, np.inf)
        np.testing.assert_allclose(in_lambda, in_k, atol=1e-6)
    # the weights of the two directions of motion differ, and the response contains only their sum
    right = quad(lambda k: abs(g_hat(k)) ** 2, 0, np.inf)[0]
    assert right > 0.9 and abs(right - 0.5) > 0.3
    # a unit vector orthogonal to everything the port reaches
    hidden = lambda k: np.sign(k) * np.conj(g_hat(-k))
    np.testing.assert_allclose(quad(lambda k: abs(hidden(k)) ** 2, -np.inf, np.inf)[0], 1.0, atol=1e-8)
    for function in (lambda lam: 1 / (lam - (0.4 + 0.9j)), lambda lam: np.exp(-0.7 * lam), lambda lam: 1.0 * (lam < 2.0),
                     lambda lam: np.exp(-0.3j * lam)):
        overlap = _complex_quad(lambda k: np.conj(function(k**2) * g_hat(k)) * hidden(k), -8.0, 8.0, points=[-np.sqrt(2.0), 0.0, np.sqrt(2.0)])
        assert abs(overlap) < 1e-8
    reached = _complex_quad(lambda k: np.conj(g_hat(k)) * g_hat(k), -np.inf, np.inf)
    np.testing.assert_allclose(reached, 1.0, atol=1e-8)                                        # while g itself is reached


def _creation(levels):
    """a^* from the first `levels` oscillator states to the first `levels` + 1."""
    out = np.zeros((levels + 1, levels))
    for n in range(levels):
        out[n + 1, n] = np.sqrt(n + 1)
    return out


def test_oscillator_map():
    rng = np.random.default_rng(402)
    for levels in (1, 4, 25):
        a_star = _creation(levels)
        np.testing.assert_allclose(a_star.T @ a_star, np.diag(np.arange(1, levels + 1)), atol=1e-12)        # a a^* = N + 1
        np.testing.assert_allclose(a_star @ a_star.T, np.diag(np.arange(0, levels + 1)), atol=1e-12)        # a^* a = N
        assert null_space(a_star).shape[1] == 0
        kernel = null_space(a_star.T)
        assert kernel.shape[1] == 1 and np.isclose(abs(kernel[0, 0]), 1.0)                                  # spanned by e_0
        assert _common_relations(a_star.astype(complex)) == -1
        x = rng.normal(size=levels) + 1j * rng.normal(size=levels)
        np.testing.assert_allclose(np.linalg.norm(a_star @ x) ** 2, np.sum(np.arange(1, levels + 1) * abs(x) ** 2), rtol=1e-12)
        assert np.linalg.norm(a_star @ x) >= np.linalg.norm(x)
        shift = np.zeros((levels + 1, levels))
        shift[np.arange(1, levels + 1), np.arange(levels)] = 1.0
        np.testing.assert_allclose(shift.T @ shift, np.eye(levels), atol=1e-15)                             # an isometry
        np.testing.assert_allclose(shift @ np.diag(np.sqrt(np.arange(1, levels + 1))), a_star, atol=1e-12)  # polar decomposition
    # heat traces: closed forms of the two geometric series, and direct sums
    n = np.arange(0, 6000)
    for t in (0.05, 0.3, 1.0, 4.0):
        np.testing.assert_allclose(1 / np.expm1(t) - 1 / (1 - np.exp(-t)), -1.0, atol=1e-10)
        np.testing.assert_allclose(np.sum(np.exp(-t * (n + 1))) - np.sum(np.exp(-t * n)), -1.0, atol=1e-9)
        np.testing.assert_allclose(np.sum(np.exp(-np.pi * t * (n + 1.0) ** 2)) - np.sum(np.exp(-np.pi * t * n.astype(float) ** 2)), -1.0, atol=1e-9)
    values = rng.normal(size=400) / (1.0 + np.arange(400)) ** 2                # a summable function on the non-negative integers
    np.testing.assert_allclose(np.sum(values[1:]) - np.sum(values), -values[0], atol=1e-12)
    # the isometric shift itself has return 1, whose heat operator has no trace
    assert np.isclose(np.trace(np.exp(-0.5) * np.eye(25)), 25 * np.exp(-0.5))


def _theta_sum(t):
    """sum over n >= 1 of exp(-pi n^2 t), by Jacobi's identity for small t."""
    direct = lambda u: np.sum(np.exp(-np.pi * np.arange(1, 60) ** 2 * u))
    if t >= 1:
        return direct(t)
    return (t**-0.5 * (1 + 2 * direct(1 / t)) - 1) / 2


def _zeta(s, terms=40):
    """Riemann zeta for Re s > 1 by Euler-Maclaurin summation."""
    n = np.arange(1, terms)
    return (np.sum(n ** (-s + 0j)) + terms ** (1 - s) / (s - 1) + terms ** (-s) / 2 + s * terms ** (-s - 1) / 12
            - s * (s + 1) * (s + 2) * terms ** (-s - 3) / 720)


def _mellin(function, s):
    """Integral over t > 0 of t^(s/2 - 1) function(t), with t = v^4 below one to remove the endpoint singularity."""
    lower = _complex_quad(lambda v: 4 * v ** (2 * s - 1) * function(v**4), 0, 1)
    upper = _complex_quad(lambda t: t ** (s / 2 - 1) * function(t), 1, np.inf)
    return lower + upper


def test_mellin_transform_of_the_heat_trace():
    np.testing.assert_allclose(_theta_sum(0.01), np.sum(np.exp(-np.pi * np.arange(1, 4000) ** 2 * 0.01)), rtol=1e-12)
    for s in (2.0, 3.0, 4.5):
        np.testing.assert_allclose(_zeta(s), zeta(s), rtol=1e-11)
    for s in (2.0, 3.0, 4.5, 2.5 + 1.3j, 1.4 - 0.7j):
        completed = np.pi ** (-s / 2) * gamma(s / 2) * _zeta(s)
        # one term: the substitution u = pi n^2 t
        for n in (1, 2, 5):
            single = _mellin(lambda t: np.exp(-np.pi * n**2 * t), s)
            np.testing.assert_allclose(single, np.pi ** (-s / 2) * n ** (-s) * gamma(s / 2), rtol=1e-7)
        # the whole trace, with its small-time part 1/(2 sqrt t) - 1/2 integrated in closed form
        remainder = lambda t: _theta_sum(t) - (t**-0.5 - 1) / 2 if t < 1 else _theta_sum(t)
        total = 1 / (s - 1) - 1 / s + _mellin(remainder, s)
        np.testing.assert_allclose(total, completed, rtol=2e-7)
        if np.isreal(s) and s >= 2:                                 # and by direct integration where the endpoint is mild
            np.testing.assert_allclose(_mellin(_theta_sum, s), completed, rtol=2e-6)
    # two copies of the oscillator have the same set of eigenvalues and twice the trace and twice the transform
    s = 3.0
    doubled = lambda t: 2 * _theta_sum(t)
    np.testing.assert_allclose(_mellin(doubled, s), 2 * np.pi ** (-s / 2) * gamma(s / 2) * _zeta(s), rtol=2e-6)
    # the interchange for complex s rests on the integrals of the absolute values, which depend only on Re s
    for s in (2.5 + 1.3j, 1.4 - 0.7j):
        sigma = s.real
        bound = sum(_mellin(lambda t, n=n: np.exp(-np.pi * n**2 * t), sigma).real for n in range(1, 400))
        tail = np.pi ** (-sigma / 2) * gamma(sigma / 2) * (_zeta(sigma) - np.sum(np.arange(1, 400) ** -sigma)).real
        np.testing.assert_allclose(bound + tail, (np.pi ** (-sigma / 2) * gamma(sigma / 2) * _zeta(sigma)).real, rtol=1e-6)
        assert abs(np.pi ** (-s / 2) * gamma(s / 2) * _zeta(s)) < bound + tail
    # the completed factor has a simple pole at s = 1 with residue one, so it is not entire
    s = 1 + 1e-6
    np.testing.assert_allclose((s - 1) * np.pi ** (-s / 2) * gamma(s / 2) * _zeta(s), 1.0, atol=1e-5)
    # at Re s = 1 the integral diverges at small t, where the trace behaves like 1 / (2 sqrt t)
    np.testing.assert_allclose(_theta_sum(1e-6) * 2 * np.sqrt(1e-6), 1.0, atol=2e-3)


def test_oscillator_response_is_minimal():
    for r in (0.3, 0.6):
        levels = 60
        vector = np.sqrt(1 - r**2) * r ** np.arange(levels)
        tail = r ** (2 * levels)                                    # squared norm of the part beyond the truncation
        np.testing.assert_allclose(np.linalg.norm(vector) ** 2 + tail, 1.0, atol=1e-12)
        h = np.diag(np.arange(1, levels + 1)).astype(complex)
        port = vector.reshape(-1, 1).astype(complex)
        for z in (0.4 + 0.3j, 2.5 - 0.8j, -3.0):
            series = (1 - r**2) * np.sum(r ** (2 * np.arange(4000)) / (np.arange(4000) + 1 - z))
            np.testing.assert_allclose((port.conj().T @ np.linalg.inv(h - z * np.eye(levels)) @ port)[0, 0], series, atol=1e-10)
        small = 12                                                  # every eigenline is reached
        assert minimal_space(h[:small, :small], port[:small]).shape[1] == small
        weights = [w[0, 0].real for _, w in response_weights(h[:small, :small], port[:small])]
        np.testing.assert_allclose(weights, (1 - r**2) * r ** (2 * np.arange(small)), rtol=1e-10)
        assert min(weights) > 0
    # a vector with one vanishing coefficient is not cyclic
    vector = np.ones(6) / np.sqrt(5)
    vector[2] = 0.0
    assert minimal_space(np.diag(np.arange(1.0, 7.0)).astype(complex), vector.reshape(-1, 1).astype(complex)).shape[1] == 5


def _eigenvalue(n, half_length, twist):
    return (2 * np.pi * n + twist) / (2 * half_length)


def test_interval_with_twisted_ends():
    rng = np.random.default_rng(403)
    for half_length, twist in ((1.0, 0.0), (0.7, 1.1), (2.3, 5.0)):
        mode = lambda n: (lambda x: np.exp(1j * _eigenvalue(n, half_length, twist) * (x + half_length)) / np.sqrt(2 * half_length))
        nodes, weights = np.polynomial.legendre.leggauss(80)
        nodes, weights = half_length * nodes, half_length * weights
        for n in (-3, 0, 2):
            np.testing.assert_allclose(mode(n)(half_length), np.exp(1j * twist) * mode(n)(-half_length), atol=1e-12)
            for m in (-3, 0, 2, 5):
                np.testing.assert_allclose(np.sum(weights * np.conj(mode(m)(nodes)) * mode(n)(nodes)), float(m == n), atol=1e-10)
            x0 = rng.uniform(-half_length, half_length)
            np.testing.assert_allclose(abs(mode(n)(x0)) ** 2, 1 / (2 * half_length), atol=1e-12)
        # a Hermitian difference operator with the twisted condition has the sampled eigenfunctions as eigenvectors
        points = 240
        step = 2 * half_length / points
        difference = np.zeros((points, points), dtype=complex)
        for j in range(points):
            difference[j, (j + 1) % points] += -1j / (2 * step) * (np.exp(1j * twist) if j + 1 == points else 1.0)
            difference[j, (j - 1) % points] += 1j / (2 * step) * (np.exp(-1j * twist) if j == 0 else 1.0)
        np.testing.assert_allclose(difference, difference.conj().T, atol=1e-12)
        grid = -half_length + step * np.arange(points)
        for n in (-2, 0, 1, 4):
            lam = _eigenvalue(n, half_length, twist)
            np.testing.assert_allclose(difference @ mode(n)(grid), np.sin(lam * step) / step * mode(n)(grid), atol=1e-9)
            np.testing.assert_allclose(np.sin(lam * step) / step, lam, atol=2e-3 * max(1.0, abs(lam) ** 3))
        # gap of the first-order operator
        values = np.abs(_eigenvalue(np.arange(-50, 51), half_length, twist))
        if twist == 0.0:
            np.testing.assert_allclose(np.sort(values)[:2], [0.0, np.pi / half_length], atol=1e-12)
        else:
            distance = min(twist % (2 * np.pi), 2 * np.pi - twist % (2 * np.pi))
            np.testing.assert_allclose(values.min(), distance / (2 * half_length), atol=1e-12)
            assert values.min() > 0
        # the smoothed point port
        n = np.arange(-4000, 4001)
        lam = _eigenvalue(n, half_length, twist)
        for tau in (0.02, 0.3, 2.0):
            mass = np.sum(np.exp(-tau * lam**2)) / (2 * half_length)
            m = np.arange(-60, 61)
            poisson = np.sum(np.exp(-m**2 * half_length**2 / tau) * np.cos(m * twist)) / np.sqrt(4 * np.pi * tau)
            np.testing.assert_allclose(mass, poisson, rtol=1e-11)
            assert abs(mass * np.sqrt(4 * np.pi * tau) - 1) <= 2 * np.exp(-half_length**2 / tau) / (1 - np.exp(-half_length**2 / tau)) + 1e-12
            for x0 in rng.uniform(-half_length, half_length, size=2):          # the mass does not depend on the point
                coefficients = np.exp(-tau * lam**2 / 2) * np.conj(mode(0)(x0) * np.exp(1j * (lam - lam[4000]) * (x0 + half_length)))
                np.testing.assert_allclose(np.sum(abs(coefficients) ** 2), mass, rtol=1e-12)
            for z in (0.3 + 0.5j, -1.0 - 0.2j):
                response = np.sum(np.exp(-tau * lam**2) / (lam - z)) / (2 * half_length)
                np.testing.assert_allclose(np.sum(abs(coefficients) ** 2 / (lam - z)), response, rtol=1e-11)
                assert (response.imag > 0) == (z.imag > 0)
                assert abs(response) <= mass / abs(z.imag) + 1e-12
                np.testing.assert_allclose(np.sum(abs(coefficients / np.sqrt(mass)) ** 2), 1.0, rtol=1e-12)
            # every line carries the weight exp(-tau lambda^2) / (2 Lambda), which is never zero
            np.testing.assert_allclose(abs(coefficients[3998:4003]) ** 2,
                                       np.exp(-tau * lam[3998:4003] ** 2) / (2 * half_length), rtol=1e-12)
            assert np.all(abs(coefficients[3998:4003]) > 0)
    # opposite twists have the same smoothed mass at every scale; the first-order response tells them apart
    n = np.arange(-4000, 4001)
    for twist in (1.1, 2.5):
        first, second = _eigenvalue(n, 0.7, twist), _eigenvalue(n, 0.7, 2 * np.pi - twist)
        for tau in (0.02, 0.3, 2.0):
            np.testing.assert_allclose(np.sum(np.exp(-tau * first**2)), np.sum(np.exp(-tau * second**2)), rtol=1e-13)
        response = lambda lam, z: np.sum(np.exp(-0.3 * lam**2) / (lam - z)) / 1.4
        z = 0.4 + 0.5j
        assert abs(response(first, z) - response(second, z)) > 0.1
        np.testing.assert_allclose(response(second, z), -response(first, -z), rtol=1e-11)        # the poles are reflected
    for twist, same in ((0.0, True), (np.pi, True), (1.1, False)):
        poles = np.sort(_eigenvalue(np.arange(-20, 21), 0.7, twist))
        reflected = np.sort(-_eigenvalue(np.arange(-21, 20), 0.7, twist) if twist else -poles)
        assert np.allclose(poles, reflected) == same
    # evaluation at a point is bounded on the scale exactly for s > 1/2
    lam = _eigenvalue(np.arange(-10**6, 10**6 + 1), 1.0, 0.4)
    partial = lambda s, cut: np.sum((1 + lam[10**6 - cut:10**6 + cut + 1] ** 2) ** (-s)) / 2.0
    assert partial(0.75, 10**6) - partial(0.75, 10**5) < 2e-3
    np.testing.assert_allclose(partial(0.5, 10**6) - partial(0.5, 10**5), np.log(10.0) / np.pi, rtol=1e-3)
    assert partial(0.0, 10**3) > 1000                                           # the unsmoothed mass grows without bound


def test_common_relations_hold_for_every_map():
    rng = np.random.default_rng(404)
    assert _common_relations(np.array([[0.7 - 0.2j]])) == 0
    assert _common_relations(_creation(9).astype(complex)) == -1
    # the twisted derivative in a band of its eigenbasis, with and without a zero mode
    for twist, zero_modes in ((0.0, 1), (0.9, 0)):
        a = np.diag(_eigenvalue(np.arange(-6, 7), 1.0, twist)).astype(complex)
        assert _common_relations(a) == 0 and null_space(a, rcond=1e-12).shape[1] == zero_modes
    # the generator of translations on a circle of growing length: index zero and a gap that closes
    gaps = []
    for points in (64, 256, 1024):
        step, shift = 0.25, np.roll(np.eye(points), 1, axis=1)
        a = (-1j / (2 * step)) * (shift - shift.T)
        assert _common_relations(a) == 0
        magnitudes = np.sort(np.abs(np.linalg.eigvalsh(a)))
        gaps.append(magnitudes[magnitudes > 1e-9][0])
        np.testing.assert_allclose(gaps[-1], np.sin(2 * np.pi / points) / step, rtol=1e-9)
    assert gaps[0] > gaps[1] > gaps[2] and gaps[2] < 0.03
    # a generic rectangular map: the index is the difference of the dimensions
    a = rng.normal(size=(7, 4)) + 1j * rng.normal(size=(7, 4))
    assert _common_relations(a) == 4 - 7
