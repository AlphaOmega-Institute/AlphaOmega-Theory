"""Finite-dimensional checks of reciprocal port responses in TAUM and the consolidated manuscript.

The models use units with hbar = 1 and fixed seeds for random matrices.
"""
import numpy as np
from scipy.linalg import expm, null_space, orth, sqrtm

from lib.linalg import spectral_blocks
from pair.model import PairModel

LABELS = {
    "UM": ["prop:passive-herglotz", "eq:passive-decay", "prop:passive-kernel", "cor:passive-sources",
           "eq:passive-signed-witness", "prop:passive-cayley", "prop:passive-schur-kernel",
           "prop:passive-spectral-temporal", "eq:passive-finite-lines", "prop:passive-two-spans",
           "prop:passive-minimal-reduction", "prop:passive-uniqueness", "eq:passive-gram",
           "cor:finite-response-defect", "prop:passive-opposite-sectors",
           "prop:passive-heat-distributional", "prop:passive-heat-discrete"],
    "TAUM": ["thm:response-Herglotz", "thm:positive-response-kernel", "cor:canonical-transition-response",
             "cor:Hilbert-block-response", "lem:Cayley-denominator", "thm:Schur-contractivity",
             "thm:positive-Schur-kernel", "def:compressed-spectral-measure", "thm:spectral-temporal-response",
             "def:response-generated-subspace", "thm:canonical-minimal-reduction", "thm:minimal-uniqueness",
             "cor:finite-response-defect", "prop:permanent-port-invisibility", "prop:opposite-sector-response",
             "thm:heat-smoothed-distributional-port", "thm:heat-smoothed-channel", "rem:heat-channel-normalization",
             "rem:positivity-boundary"],
}


def _complex(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _hermitian(rng, n):
    raw = _complex(rng, n, n)
    return (raw + raw.conj().T) / 2


def _resolvent(h, z):
    return np.linalg.inv(h - z * np.eye(len(h)))


def _response(h, port, z):
    return port.conj().T @ _resolvent(h, z) @ port


def _imag(m):
    return (m - m.conj().T) / 2j


def _upper(rng):
    return rng.normal() + 1j * rng.uniform(0.05, 2.0)


def _gram_from_kernel(kernel, points):
    """Block matrix with blocks kernel(z_i, z_j)."""
    return np.block([[kernel(zi, zj) for zj in points] for zi in points])


def test_herglotz_response_and_positive_kernel():
    rng = np.random.default_rng(101)
    for _ in range(40):
        n, p = rng.integers(3, 8), rng.integers(1, 4)
        h, port = _hermitian(rng, n), _complex(rng, n, p)
        z, w = _upper(rng), _upper(rng)
        m = _response(h, port, z)
        np.testing.assert_allclose(_response(h, port, np.conj(z)), m.conj().T, atol=1e-10)
        r = _resolvent(h, z)
        np.testing.assert_allclose(_imag(m), z.imag * port.conj().T @ r @ r.conj().T @ port, atol=1e-10)
        assert np.linalg.eigvalsh(_imag(m)).min() > -1e-10
        kernel = lambda a, b: (_response(h, port, a) - _response(h, port, b).conj().T) / (a - np.conj(b))
        np.testing.assert_allclose(kernel(z, w), port.conj().T @ r @ _resolvent(h, w).conj().T @ port, atol=1e-10)
        points = [_upper(rng) for _ in range(4)]
        gram = _gram_from_kernel(kernel, points)
        vectors = np.hstack([_resolvent(h, zj).conj().T @ port for zj in points])
        np.testing.assert_allclose(gram, vectors.conj().T @ vectors, atol=1e-9)
        assert np.linalg.eigvalsh((gram + gram.conj().T) / 2).min() > -1e-9
        # decay off the real axis, in both half-planes
        for point in (z, np.conj(w), 1j * rng.uniform(0.1, 50.0)):
            assert np.linalg.norm(_response(h, port, point), 2) <= np.linalg.norm(port, 2) ** 2 / abs(point.imag) + 1e-12


def test_positivity_does_not_characterize_bounded_port_responses():
    """M(z) = z I has positive imaginary part and the positive kernel I, yet it grows where a response must decay."""
    rng = np.random.default_rng(110)
    identity = np.eye(2)
    linear = lambda z: z * identity
    kernel = lambda a, b: (linear(a) - linear(b).conj().T) / (a - np.conj(b))
    points = [_upper(rng) for _ in range(5)]
    for a in points:
        assert np.linalg.eigvalsh(_imag(linear(a))).min() > 0
        for b in points:
            np.testing.assert_allclose(kernel(a, b), identity, atol=1e-12)
    assert np.linalg.eigvalsh(_gram_from_kernel(kernel, points)).min() > -1e-12
    # y ||M(iy)|| is bounded by ||Gamma||^2 for a response and equals y^2 for z I
    h, port = _hermitian(rng, 4), _complex(rng, 4, 2)
    bound = np.linalg.norm(port, 2) ** 2
    for y in (1.0, 10.0, 1e3, 1e6):
        assert y * np.linalg.norm(_response(h, port, 1j * y), 2) <= bound + 1e-9
        assert y * np.linalg.norm(linear(1j * y), 2) == y**2
    assert 1e6**2 > bound


def test_form_selfadjoint_generator_with_a_negative_imaginary_part():
    """A generator self-adjoint for an indefinite form, with real spectrum, whose response is not Herglotz."""
    generator = np.array([[0.0, np.sqrt(2)], [-np.sqrt(2), 3.0]])
    form = np.diag([1.0, -1.0])
    np.testing.assert_allclose(form @ generator, (form @ generator).T, atol=1e-15)      # self-adjoint for the form
    assert not np.allclose(generator, generator.T)                                      # not Hilbert-self-adjoint
    np.testing.assert_allclose(np.sort(np.linalg.eigvals(generator).real), [1.0, 2.0], atol=1e-12)
    np.testing.assert_allclose(np.linalg.eigvals(generator).imag, 0.0, atol=1e-12)
    response = lambda z: np.linalg.inv(generator - z * np.eye(2))[0, 0]
    for z in (0.3 + 0.7j, -1.0 + 0.2j, 2.0 + 0.5j, 5.0 - 1.0j):
        np.testing.assert_allclose(response(z), 2 / (1 - z) - 1 / (2 - z), atol=1e-12)
    np.testing.assert_allclose(response(2.0 + 0.5j), -1.6 - 1.2j, atol=1e-12)           # imaginary part -6/5
    assert response(2.0 + 0.5j).imag < 0 < (2.0 + 0.5j).imag


def test_two_sources_of_positive_responses():
    rng = np.random.default_rng(102)
    for _ in range(30):
        # (a) the return A*A of a rectangular map
        rows, cols, p = rng.integers(2, 7), rng.integers(2, 6), rng.integers(1, 3)
        a, port = _complex(rng, rows, cols), _complex(rng, cols, p)
        z = _upper(rng)
        assert np.linalg.eigvalsh(_imag(_response(a.conj().T @ a, port, z))).min() > -1e-10
        # (b) a port into the retained block of a self-adjoint block operator
        n_p, n_q = rng.integers(2, 5), rng.integers(2, 6)
        h_p, h_q, c = _hermitian(rng, n_p), _hermitian(rng, n_q), _complex(rng, n_q, n_p)
        h = np.block([[h_p, c.conj().T], [c, h_q]])
        port_p = _complex(rng, n_p, p)
        embedded = np.vstack([port_p, np.zeros((n_q, p))])
        reduced = np.linalg.inv(h_p - z * np.eye(n_p) - c.conj().T @ _resolvent(h_q, z) @ c)
        full = _response(h, embedded, z)
        np.testing.assert_allclose(full, port_p.conj().T @ reduced @ port_p, atol=1e-9)
        assert np.linalg.eigvalsh(_imag(full)).min() > -1e-10


def test_cayley_transform_and_schur_kernel():
    rng = np.random.default_rng(103)
    for _ in range(40):
        n, p = rng.integers(3, 8), rng.integers(1, 4)
        h, port = _hermitian(rng, n), _complex(rng, n, p)
        identity = np.eye(p)
        schur = lambda a: (_response(h, port, a) - 1j * identity) @ np.linalg.inv(_response(h, port, a) + 1j * identity)
        z, w = _upper(rng), _upper(rng)
        m = _response(h, port, z)
        assert np.linalg.norm(np.linalg.inv(m + 1j * identity), 2) <= 1 + 1e-12
        s = schur(z)
        defect = identity - s.conj().T @ s
        expected = 4 * np.linalg.inv(m.conj().T - 1j * identity) @ _imag(m) @ np.linalg.inv(m + 1j * identity)
        np.testing.assert_allclose(defect, expected, atol=1e-9)
        assert np.linalg.norm(s, 2) <= 1 + 1e-12
        k_m = lambda a, b: (_response(h, port, a) - _response(h, port, b).conj().T) / (a - np.conj(b))
        k_s = lambda a, b: (identity - schur(a) @ schur(b).conj().T) / (-1j * (a - np.conj(b)))
        factor = 2 * np.linalg.inv(m + 1j * identity) @ k_m(z, w) @ np.linalg.inv(_response(h, port, w).conj().T - 1j * identity)
        np.testing.assert_allclose(k_s(z, w), factor, atol=1e-9)
        gram = _gram_from_kernel(k_s, [_upper(rng) for _ in range(4)])
        assert np.linalg.eigvalsh((gram + gram.conj().T) / 2).min() > -1e-9


def test_spectral_weights_and_time_domain():
    rng = np.random.default_rng(104)
    for _ in range(30):
        n, p = rng.integers(3, 8), rng.integers(1, 4)
        h, port = _hermitian(rng, n), _complex(rng, n, p)
        values, vectors = np.linalg.eigh(h)
        weights = [port.conj().T @ np.outer(v, v.conj()) @ port for v in vectors.T]
        np.testing.assert_allclose(sum(weights), port.conj().T @ port, atol=1e-10)
        assert all(np.linalg.eigvalsh(r).min() > -1e-10 for r in weights)
        z = _upper(rng)
        np.testing.assert_allclose(_response(h, port, z), sum(r / (lam - z) for r, lam in zip(weights, values)), atol=1e-9)
        correlation = lambda t: port.conj().T @ expm(-1j * t * h) @ port
        t = rng.uniform(-5, 5)
        np.testing.assert_allclose(correlation(t), sum(r * np.exp(-1j * t * lam) for r, lam in zip(weights, values)), atol=1e-9)
        times = rng.uniform(-4, 4, size=4)
        gram = np.block([[correlation(ti - tj) for tj in times] for ti in times])
        evolved = np.hstack([expm(1j * tj * h) @ port for tj in times])
        np.testing.assert_allclose(gram, evolved.conj().T @ evolved, atol=1e-9)
        assert np.linalg.eigvalsh((gram + gram.conj().T) / 2).min() > -1e-9


def _reached_spaces(h, port, rng):
    """Orthonormal bases of the resolvent span and of the time span."""
    n = len(h)
    resolvent = np.hstack([_resolvent(h, z) @ port for z in [_upper(rng) * s for s in (1, -1) for _ in range(n)]])
    time = np.hstack([expm(-1j * t * h) @ port for t in rng.uniform(-9, 9, size=2 * n)])
    return orth(resolvent, rcond=1e-9), orth(time, rcond=1e-9)


def test_reached_space_minimal_reduction_and_blindness():
    rng = np.random.default_rng(105)
    for _ in range(25):
        n_a, n_b, p = rng.integers(3, 6), rng.integers(2, 5), rng.integers(1, 3)
        # H is block diagonal in a hidden basis; the port lives in the first block
        mixing, _ = np.linalg.qr(_complex(rng, n_a + n_b, n_a + n_b))
        h_a, h_b = _hermitian(rng, n_a), _hermitian(rng, n_b)
        h = mixing @ np.block([[h_a, np.zeros((n_a, n_b))], [np.zeros((n_b, n_a)), h_b]]) @ mixing.conj().T
        port = mixing[:, :n_a] @ _complex(rng, n_a, p)
        res, time = _reached_spaces(h, port, rng)
        assert res.shape[1] == time.shape[1] <= n_a
        np.testing.assert_allclose(res @ res.conj().T, time @ time.conj().T, atol=1e-7)
        projection = res @ res.conj().T
        np.testing.assert_allclose(projection @ h, h @ projection, atol=1e-7)          # reduces H
        np.testing.assert_allclose(projection @ port, port, atol=1e-8)                 # contains Ran Gamma
        z = _upper(rng)
        h_min, port_min = res.conj().T @ h @ res, res.conj().T @ port
        np.testing.assert_allclose(_response(h, port, z), _response(h_min, port_min, z), atol=1e-8)
        # invisible vectors are exactly the complement; changing H there changes nothing at the port
        complement = null_space(res.conj().T)
        for t in rng.uniform(-6, 6, size=3):
            np.testing.assert_allclose(port.conj().T @ expm(-1j * t * h) @ complement, 0, atol=1e-8)
        moved = h + complement @ _hermitian(rng, complement.shape[1]) @ complement.conj().T
        np.testing.assert_allclose(_response(moved, port, z), _response(h, port, z), atol=1e-8)


def test_unitary_uniqueness_of_minimal_realizations():
    rng = np.random.default_rng(106)
    for _ in range(25):
        n, p = rng.integers(3, 7), rng.integers(1, 3)
        h1, port1 = _hermitian(rng, n), _complex(rng, n, p)
        assert _reached_spaces(h1, port1, rng)[0].shape[1] == n                      # minimal
        hidden, _ = np.linalg.qr(_complex(rng, n, n))
        h2, port2 = hidden @ h1 @ hidden.conj().T, hidden @ port1
        points = [_upper(rng) for _ in range(n)]
        points += [np.conj(z) for z in points]                   # every point together with its conjugate
        for z in points:
            np.testing.assert_allclose(_response(h1, port1, z), _response(h2, port2, z), atol=1e-9)
        # the unitary is determined by the responses: it sends R_1(z) Gamma_1 h to R_2(z) Gamma_2 h
        first = np.hstack([_resolvent(h1, z) @ port1 for z in points])
        second = np.hstack([_resolvent(h2, z) @ port2 for z in points])
        np.testing.assert_allclose(first.conj().T @ first, second.conj().T @ second, atol=1e-8)     # equal Gram matrices
        w = second @ np.linalg.pinv(first)
        np.testing.assert_allclose(w.conj().T @ w, np.eye(n), atol=1e-6)
        np.testing.assert_allclose(w @ h1, h2 @ w, atol=1e-6)
        np.testing.assert_allclose(w @ port1, port2, atol=1e-6)
        np.testing.assert_allclose(w, hidden, atol=1e-6)
        # the Gram formula at every pair of points that is not a conjugate pair
        m1 = lambda a: _response(h1, port1, a)
        gram = lambda a, b: (_resolvent(h1, a) @ port1).conj().T @ (_resolvent(h1, b) @ port1)
        for a in points[:3] + points[n:n + 3]:
            for b in points:
                if abs(b - np.conj(a)) > 1e-12:
                    np.testing.assert_allclose(gram(a, b), (m1(b) - m1(a).conj().T) / (b - np.conj(a)), atol=1e-7)
        # the conjugate pair w = conj(z): the formula is 0/0 there, and the Gram entry is the limit of the
        # quotient, which is the derivative of the response at conj(z).  Both come from the response alone.
        z = points[0]
        assert np.conj(z) in points
        entry = gram(z, np.conj(z))
        np.testing.assert_allclose(entry, port1.conj().T @ _resolvent(h1, np.conj(z)) @ _resolvent(h1, np.conj(z)) @ port1, atol=1e-9)
        errors = []
        for step in (1e-2, 1e-3, 1e-4):
            v = np.conj(z) + step * (0.6 - 0.8j)
            errors.append(np.linalg.norm((m1(v) - m1(z).conj().T) / (v - np.conj(z)) - entry))
        assert errors[0] > errors[1] > errors[2] and errors[2] < 1e-3 * max(1.0, np.linalg.norm(entry)) / z.imag
        m2 = lambda a: _response(h2, port2, a)
        np.testing.assert_allclose((_resolvent(h2, z) @ port2).conj().T @ (_resolvent(h2, np.conj(z)) @ port2), entry, atol=1e-8)
        # equal points in the intertwining step: R(zeta) R(z) Gamma at zeta = z is the limit of the difference quotient
        square = _resolvent(h1, z) @ _resolvent(h1, z) @ port1
        errors = []
        for step in (1e-2, 1e-3, 1e-4):
            zeta = z + step * (0.6 + 0.8j) * z.imag
            errors.append(np.linalg.norm((_resolvent(h1, zeta) - _resolvent(h1, z)) @ port1 / (zeta - z) - square))
        assert errors[0] > errors[1] > errors[2] and errors[2] < 1e-3 * np.linalg.norm(square) / z.imag
        np.testing.assert_allclose(w @ square, _resolvent(h2, z) @ _resolvent(h2, z) @ port2, atol=1e-5 * max(1.0, np.linalg.norm(square)))
        # a larger system with the same response has the same minimal part
        extra = rng.integers(1, 4)
        big = np.block([[h1, np.zeros((n, extra))], [np.zeros((extra, n)), _hermitian(rng, extra)]])
        big_port = np.vstack([port1, np.zeros((extra, p))])
        np.testing.assert_allclose(_response(big, big_port, z), m1(z), atol=1e-9)
        assert _reached_spaces(big, big_port, rng)[0].shape[1] == n


def _count(weights):
    return sum(np.linalg.matrix_rank(r, tol=1e-9) for r in weights) - np.linalg.matrix_rank(sum(weights), tol=1e-9)


def _weights(h, port):
    """Nonzero weights Gamma^* E_k Gamma.  Each eigenvector belongs to exactly one block of spectral_blocks."""
    out = []
    for _, block in spectral_blocks(h):
        piece = block.conj().T @ port
        weight = piece.conj().T @ piece
        if np.linalg.norm(weight) > 1e-10:
            out.append(weight)
    return out


def test_nearby_energies_are_counted_once():
    """Lines closer than a relative tolerance but farther apart than the spectral tolerance stay separate."""
    port = np.eye(2, dtype=complex)
    close = np.diag([1.0, 1.000001]).astype(complex)
    weights = _weights(close, port)
    assert [np.linalg.matrix_rank(r, tol=1e-9) for r in weights] == [1, 1]
    np.testing.assert_allclose(sum(weights), np.eye(2), atol=1e-12)                 # the weights sum to Gamma^* Gamma
    assert _count(weights) == 0 and _reached_by_spectral_projections(close, port).shape[1] == 2
    # a true degeneracy is one line of rank two, and the count is again zero
    weights = _weights(np.eye(2, dtype=complex), port)
    assert [np.linalg.matrix_rank(r, tol=1e-9) for r in weights] == [2] and _count(weights) == 0
    # one vector spread over three lines, two of them close: three rank-one weights and two indirect directions
    spread = np.diag([1.0, 1.000001, 2.0]).astype(complex)
    vector = np.ones((3, 1), dtype=complex) / np.sqrt(3)
    weights = _weights(spread, vector)
    assert len(weights) == 3 and _count(weights) == 2
    assert _reached_by_spectral_projections(spread, vector).shape[1] == 3
    np.testing.assert_allclose(sum(weights), vector.conj().T @ vector, atol=1e-12)


def test_finite_response_relations_and_enlarged_ports():
    rng = np.random.default_rng(107)
    for theta in (np.pi / 8, 0.2, 0.6):
        pair = PairModel(theta)
        shell = pair.F2
        h = shell.conj().T @ pair.H @ shell
        port = shell.conj().T @ pair.A
        weights = _weights(h, port)
        assert _count(weights) == 1
        # relations among the ranges of the weights
        ranges = [orth(r, rcond=1e-9) for r in weights]
        relations = sum(b.shape[1] for b in ranges) - np.linalg.matrix_rank(np.hstack(ranges), tol=1e-9)
        assert relations == 1
        np.testing.assert_allclose(orth(sum(weights), rcond=1e-9) @ orth(sum(weights), rcond=1e-9).conj().T,
                                   orth(np.hstack(ranges), rcond=1e-9) @ orth(np.hstack(ranges), rcond=1e-9).conj().T, atol=1e-8)
        # an invertible change of port coordinates preserves the count
        change = _complex(rng, 2, 2)
        assert _count(_weights(h, port @ change)) == 1
        # adding the endpoint preparation fills the indirect line
        endpoint = shell.conj().T @ PairModel(0.0).A
        assert _count(_weights(h, endpoint)) == 0
        both = np.hstack([port, endpoint])
        assert _count(_weights(h, both)) == 0 and np.linalg.matrix_rank(both, tol=1e-9) == 3
        # the port is injective exactly when the summed weight has trivial kernel
        assert np.linalg.matrix_rank(sum(weights), tol=1e-9) == 2 and np.linalg.matrix_rank(sum(_weights(h, both)), tol=1e-9) == 3
    # adding a port can also raise the count: the endpoint pi/4 plane plus one vector mixing two energies
    pair = PairModel(np.pi / 4)
    shell = pair.F2
    h = shell.conj().T @ pair.H @ shell
    port = shell.conj().T @ pair.A
    assert _count(_weights(h, port)) == 0
    mixed = shell.conj().T @ ((pair.e + pair.w0) / np.sqrt(2))
    assert _count(_weights(h, np.column_stack([port, mixed]))) == 1


def _reached_by_spectral_projections(h, port):
    """Orthonormal basis of the span of E_k Gamma over the spectral blocks of a Hermitian matrix."""
    pieces = [block @ (block.conj().T @ port) for _, block in spectral_blocks(h)]
    return orth(np.hstack(pieces), rcond=1e-9)


def _sector_model(rng, levels, nu, ell, weights=None):
    """g, then a block with L = +ell, then its mirror image with L = -ell.  H0 is the same on both blocks."""
    k = len(levels)
    h0 = np.diag(np.r_[0.0, levels, levels]).astype(complex)
    lz = np.diag(np.r_[0.0, ell * np.ones(k), -ell * np.ones(k)])
    if weights is None:
        vector = _complex(rng, k)
        vector /= np.linalg.norm(vector)
    else:
        vector = np.sqrt(np.asarray(weights, dtype=complex))
    g = np.eye(2 * k + 1)[:, 0].astype(complex)
    q_plus = np.r_[0, vector, np.zeros(k)]
    q_minus = np.r_[0, np.zeros(k), vector]
    return h0 - nu * lz, h0, g, q_plus, q_minus


def test_opposite_symmetry_sectors():
    rng = np.random.default_rng(108)
    nu, ell = 0.37, 2.0
    for _ in range(12):
        levels = np.sort(rng.uniform(0.5, 6.0, size=rng.integers(1, 5)))
        h_nu, h0, g, q_plus, q_minus = _sector_model(rng, levels, nu, ell)
        e = (q_plus + q_minus) / np.sqrt(2)
        d = (q_plus - q_minus) / (1j * np.sqrt(2))
        np.testing.assert_allclose(np.array([[np.vdot(a, b) for b in (g, e, d)] for a in (g, e, d)]), np.eye(3), atol=1e-12)
        port = np.column_stack([g, e])
        reached = _reached_by_spectral_projections(h_nu, port)
        p_min = reached @ reached.conj().T
        p_direct = port @ port.conj().T
        for t in np.r_[rng.uniform(0, 12, size=6), np.pi / (2 * nu * ell)]:
            c0 = np.vdot(q_plus, expm(-1j * t * h0) @ q_plus)
            np.testing.assert_allclose(c0, np.vdot(q_minus, expm(-1j * t * h0) @ q_minus), atol=1e-10)
            psi = expm(-1j * t * h_nu) @ e
            np.testing.assert_allclose(np.linalg.norm(p_min @ psi), 1.0, atol=1e-9)        # the orbit stays in the reached space
            np.testing.assert_allclose(np.vdot(e, psi), c0 * np.cos(nu * ell * t), atol=1e-9)
            np.testing.assert_allclose(np.vdot(d, psi), -c0 * np.sin(nu * ell * t), atol=1e-9)
            indirect = np.linalg.norm((p_min - p_direct) @ psi) ** 2
            np.testing.assert_allclose(indirect, 1 - abs(c0) ** 2 * np.cos(nu * ell * t) ** 2, atol=1e-9)
            outside = 1 - sum(abs(np.vdot(x, psi)) ** 2 for x in (g, e, d))
            np.testing.assert_allclose(outside, 1 - abs(c0) ** 2, atol=1e-9)
        assert np.isclose(np.linalg.norm((p_min - p_direct) @ expm(-1j * np.pi / (2 * nu * ell) * h_nu) @ e) ** 2, 1.0)
        # with distinct shifted lines the analysis vector d is reachable
        shifted = np.r_[levels - nu * ell, levels + nu * ell]
        if np.min(np.diff(np.sort(shifted))) > 1e-3:
            assert np.isclose(np.linalg.norm(p_min @ d), 1.0, atol=1e-8)
    # overlapping branches: weight one half at 0 and at 2 nu ell in each sector
    h_nu, h0, g, q_plus, q_minus = _sector_model(rng, np.array([0.0, 2 * nu * ell]), nu, ell, weights=[0.5, 0.5])
    e = (q_plus + q_minus) / np.sqrt(2)
    d = (q_plus - q_minus) / (1j * np.sqrt(2))
    port = np.column_stack([g, e])
    reached = _reached_by_spectral_projections(h_nu, port)
    assert reached.shape[1] == 4                                 # the ground and three excitation lines
    assert np.linalg.norm(reached @ reached.conj().T @ d) < 1 - 1e-3
    psi = expm(-1j * 1.3 * h_nu) @ e
    np.testing.assert_allclose(np.linalg.norm(reached.conj().T @ psi), 1.0, atol=1e-10)


def test_heat_smoothed_ports():
    rng = np.random.default_rng(109)
    count, p, tau = 300, 2, 0.02
    levels = np.arange(count) - 40.0
    coefficients = _complex(rng, count, p) * (1 + np.abs(levels))[:, None] ** 2        # grow polynomially
    damping = np.exp(-tau * levels**2 / 2)
    port = damping[:, None] * coefficients.conj()                # (Gamma h)_n = exp(-tau l_n^2/2) <gamma_n, h>
    summed = np.sum(np.exp(-tau * levels**2) * np.linalg.norm(coefficients, axis=1) ** 2)
    assert np.isclose(np.linalg.norm(port, "fro") ** 2, summed)
    h = np.diag(levels)
    for z in (0.4 + 0.3j, -7.0 - 1.1j):
        series = sum(np.exp(-tau * lam**2) / (lam - z) * np.outer(c, c.conj()) for lam, c in zip(levels, coefficients))
        np.testing.assert_allclose(_response(h.astype(complex), port, z), series, atol=1e-8 * abs(series).max())
    mass = sum(np.exp(-tau * lam**2) * np.outer(c, c.conj()) for lam, c in zip(levels, coefficients))
    np.testing.assert_allclose(port.conj().T @ port, mass, atol=1e-8 * abs(mass).max())
    assert np.isclose(np.trace(mass).real, summed)
    # the tail beyond the truncation is negligible: the last retained terms are already below rounding
    assert np.exp(-tau * levels[-1] ** 2) * np.linalg.norm(coefficients[-1]) ** 2 < 1e-12 * summed
    # a distributional trace bounded on the scale X_s: |gamma x| <= ||gamma||_s ||(1 + H^2)^(s/2) x||
    s = 3.0
    scale = (1 + levels**2) ** (s / 2)
    trace = _complex(rng, p, count) / (1 + np.abs(levels))[None, :] ** 0.6 * scale[None, :] / np.sqrt(count)
    norm_on_scale = np.linalg.norm(trace / scale[None, :], 2)
    smoothed = trace * damping[None, :]                          # gamma exp(-tau H^2/2), bounded on X
    assert np.linalg.norm(smoothed, 2) <= norm_on_scale * (scale * damping).max() + 1e-12
    assert np.linalg.norm(trace, 2) > 50 * np.linalg.norm(smoothed, 2)       # unsmoothed, it is far larger at this truncation
    heat_port = smoothed.conj().T
    assert np.linalg.eigvalsh(_imag(_response(h.astype(complex), heat_port, 0.3 + 0.2j))).min() > -1e-10
    # normalization: scalar port of unit mass, and the polar part of a general port
    scalar = port[:, :1]
    unit = scalar / np.linalg.norm(scalar)
    assert np.isclose((unit.conj().T @ unit)[0, 0].real, 1.0)
    weak = port @ np.diag([1.0, 1e-6])                           # one direction of very small weight
    modulus = sqrtm(weak.conj().T @ weak)
    polar = weak @ np.linalg.inv(modulus)
    np.testing.assert_allclose(polar.conj().T @ polar, np.eye(p), atol=1e-6)
    np.testing.assert_allclose(polar @ modulus, weak, atol=1e-10 * np.linalg.norm(weak))
