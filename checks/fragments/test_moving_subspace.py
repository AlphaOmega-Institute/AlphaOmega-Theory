"""Checks for the fragment on moving subspaces (manuscripts/consolidated/fragments/moving_subspace.tex).

The fragment carries TAUM's sections "Moving decompositions and projector curvature" and
"Characteristic classes of the projector bundle" into the consolidated manuscript.  Its
own labels are not in Universal_Model.tex yet, so this file cites the TAUM labels of the
statements it tests and the labels of the consolidated manuscript that the worked example
uses.  hbar = 1 and random numbers use fixed seeds.

A family of projections is built as Pi(x) = U(x) P U(x)^* with U(x) the exponential of a
polynomial in the parameters.  Its first derivatives are exact (Frechet derivative of the
matrix exponential).  Second derivatives, where a statement needs them, are central
differences of exact first derivatives.
"""
import itertools

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, expm_frechet, sqrtm

from pair.model import PairModel

LABELS = {
    "TAUM": ["lem:projection-tangent", "def:tangent-channel", "thm:tangent-realization",
             "def:Grassmann-connection", "thm:projector-curvature", "def:mixed-return", "cor:curvature-AV",
             "def:local-gauge-data", "thm:frame-covariance", "def:holonomy", "eq:fiberwise-Laplacian",
             "thm:chern-character", "eq:higher-chern-character", "cor:line-flux", "eq:ch2-return"],
    "UM": ["eq:prepared-schur", "eq:pair-invariance-amplitude", "eq:index-graded"],
}

STEP = 1e-5


def _complex(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _hermitian(rng, n):
    raw = _complex(rng, n, n)
    return (raw + raw.conj().T) / 2


def _comm(a, b):
    return a @ b - b @ a


class Family:
    """Pi(x) = U(x) P U(x)^* on C^n, of rank r, with p real parameters and exact first derivatives."""

    def __init__(self, rng, n, r, p, strength=0.5):
        self.n, self.r, self.p = n, r, p
        self.linear = [strength * _hermitian(rng, n) for _ in range(p)]
        self.quadratic = [0.3 * strength * _hermitian(rng, n) for _ in range(p)]
        basis, _ = np.linalg.qr(_complex(rng, n, n))
        self.xi0 = basis[:, :r]

    def _unitary(self, x):
        generator = 1j * sum(x[a] * self.linear[a] + x[a] ** 2 * self.quadratic[a] for a in range(self.p))
        derivatives = [expm_frechet(generator, 1j * (self.linear[a] + 2 * x[a] * self.quadratic[a]), compute_expm=False)
                       for a in range(self.p)]
        return expm(generator), derivatives

    def frame(self, x):
        u, du = self._unitary(np.asarray(x, dtype=float))
        return u @ self.xi0, [d @ self.xi0 for d in du]

    def proj(self, x):
        xi, _ = self.frame(x)
        return xi @ xi.conj().T

    def dproj(self, x):
        xi, dxi = self.frame(x)
        return xi @ xi.conj().T, [d @ xi.conj().T + xi @ d.conj().T for d in dxi]

    def curvature(self, x, a, b):
        """F_ab = Pi [d_a Pi, d_b Pi] Pi as an operator on C^n."""
        pi, dpi = self.dproj(x)
        return pi @ _comm(dpi[a], dpi[b]) @ pi

    def connection(self, x):
        xi, dxi = self.frame(x)
        return [xi.conj().T @ d for d in dxi]


def _shift(x, a, h):
    y = np.array(x, dtype=float)
    y[a] += h
    return y


def _diff(function, x, a):
    """Central difference of a matrix-valued function of the parameters."""
    return (function(_shift(x, a, STEP)) - function(_shift(x, a, -STEP))) / (2 * STEP)


def test_direct_rotation_of_nearby_projections():
    rng = np.random.default_rng(201)
    for _ in range(40):
        n, r = rng.integers(3, 8), rng.integers(1, 3)
        basis, _ = np.linalg.qr(_complex(rng, n, n))
        p = basis[:, :r] @ basis[:, :r].conj().T
        turn = expm(1j * rng.uniform(0.05, 0.45) * _hermitian(rng, n))
        q = turn @ p @ turn.conj().T
        if np.linalg.norm(p - q, 2) >= 0.98:
            continue
        defect = (p - q) @ (p - q)
        np.testing.assert_allclose(defect, p + q - p @ q - q @ p, atol=1e-12)
        np.testing.assert_allclose(p @ defect, defect @ p, atol=1e-12)
        eye = np.eye(n)
        raw = q @ p + (eye - q) @ (eye - p)
        np.testing.assert_allclose(raw.conj().T @ raw, eye - defect, atol=1e-12)
        np.testing.assert_allclose(raw @ raw.conj().T, eye - defect, atol=1e-12)
        rotation = raw @ np.linalg.inv(sqrtm(eye - defect))
        np.testing.assert_allclose(rotation.conj().T @ rotation, eye, atol=1e-10)
        np.testing.assert_allclose(rotation @ p @ rotation.conj().T, q, atol=1e-10)
    # projections of different rank are never closer than one
    for _ in range(20):
        n = rng.integers(3, 8)
        basis, _ = np.linalg.qr(_complex(rng, n, n))
        other, _ = np.linalg.qr(_complex(rng, n, n))
        one, two = basis[:, :1] @ basis[:, :1].conj().T, other[:, :2] @ other[:, :2].conj().T
        assert np.linalg.norm(one - two, 2) >= 1 - 1e-12


def test_tangent_map_of_a_moving_subspace():
    rng = np.random.default_rng(202)
    for _ in range(25):
        n, r, p = rng.integers(3, 7), rng.integers(1, 3), 3
        family = Family(rng, n, r, p)
        x = rng.uniform(-0.6, 0.6, size=p)
        pi, dpi = family.dproj(x)
        for a in range(p):                                         # the exact derivative agrees with a difference quotient
            np.testing.assert_allclose(dpi[a], _diff(family.proj, x, a), atol=1e-8)
        direction = rng.normal(size=p)
        y_pi = sum(direction[a] * dpi[a] for a in range(p))
        eye = np.eye(n)
        perp, symmetry = eye - pi, 2 * pi - eye
        np.testing.assert_allclose(y_pi, y_pi.conj().T, atol=1e-12)
        np.testing.assert_allclose(pi @ y_pi @ pi, 0, atol=1e-12)
        np.testing.assert_allclose(perp @ y_pi @ perp, 0, atol=1e-12)
        np.testing.assert_allclose(symmetry @ y_pi + y_pi @ symmetry, 0, atol=1e-12)
        tangent = perp @ y_pi @ pi                                 # B_Y, also the differential d_Y on K
        np.testing.assert_allclose(tangent.conj().T, pi @ y_pi @ perp, atol=1e-12)
        np.testing.assert_allclose(y_pi, tangent + tangent.conj().T, atol=1e-12)
        np.testing.assert_allclose(symmetry @ y_pi, tangent.conj().T - tangent, atol=1e-12)
        laplacian = tangent.conj().T @ tangent + tangent @ tangent.conj().T
        np.testing.assert_allclose(y_pi @ y_pi, laplacian, atol=1e-12)
        np.testing.assert_allclose((symmetry @ y_pi) @ (symmetry @ y_pi), -laplacian, atol=1e-12)
        # the two adjoints of Section 3, for the symmetry of this point
        krein = lambda t: symmetry @ t.conj().T @ symmetry
        np.testing.assert_allclose(krein(y_pi), -y_pi, atol=1e-12)
        np.testing.assert_allclose(krein(symmetry @ y_pi), symmetry @ y_pi, atol=1e-12)
        np.testing.assert_allclose((symmetry @ y_pi).conj().T, -symmetry @ y_pi, atol=1e-12)
        # block form in a basis adapted to the splitting, and paired nonzero spectra of the two returns
        xi, _ = family.frame(x)
        full, _ = np.linalg.qr(np.hstack([xi, _complex(rng, n, n - r)]))
        blocks = full.conj().T @ y_pi @ full
        np.testing.assert_allclose(blocks[:r, :r], 0, atol=1e-10)
        np.testing.assert_allclose(blocks[r:, r:], 0, atol=1e-10)
        lower = blocks[r:, :r]
        first = np.sort(np.linalg.eigvalsh(lower.conj().T @ lower))[::-1]
        second = np.sort(np.linalg.eigvalsh(lower @ lower.conj().T))[::-1]
        k = min(len(first), len(second))
        np.testing.assert_allclose(first[:k], second[:k], atol=1e-10)
        np.testing.assert_allclose(np.r_[first[k:], second[k:]], 0, atol=1e-10)


def _section(rng, family):
    """A smooth section s(x) = Pi(x) w(x) and its exact first derivatives."""
    n, p = family.n, family.p
    w0, w1, w2 = _complex(rng, n), _complex(rng, p, n), 0.5 * _complex(rng, p, n)

    def value(x):
        return family.proj(x) @ (w0 + x @ w1 + (x**2) @ w2)

    def derivative(x, a):
        pi, dpi = family.dproj(x)
        return dpi[a] @ (w0 + x @ w1 + (x**2) @ w2) + pi @ (w1[a] + 2 * x[a] * w2[a])

    return value, derivative


def test_curvature_is_the_antisymmetrized_mixed_return():
    rng = np.random.default_rng(203)
    for _ in range(12):
        n, r, p = rng.integers(3, 7), rng.integers(1, 3), 3
        family = Family(rng, n, r, p)
        x = rng.uniform(-0.5, 0.5, size=p)
        pi, dpi = family.dproj(x)
        perp = np.eye(n) - pi
        tangent = [perp @ d @ pi for d in dpi]
        value, derivative = _section(rng, family)
        s = value(x)
        np.testing.assert_allclose(pi @ s, s, atol=1e-12)
        covariant = lambda a: (lambda y: family.proj(y) @ derivative(y, a))        # x -> nabla_a s
        for a, b in itertools.permutations(range(p), 2):
            curvature = family.curvature(x, a, b)
            mixed = tangent[a].conj().T @ tangent[b] - tangent[b].conj().T @ tangent[a]
            np.testing.assert_allclose(curvature, mixed, atol=1e-12)
            np.testing.assert_allclose(curvature.conj().T, -curvature, atol=1e-12)
            # definition: [nabla_a, nabla_b] s, with the outer derivative taken numerically
            second = pi @ _diff(covariant(b), x, a) - pi @ _diff(covariant(a), x, b)
            np.testing.assert_allclose(second, curvature @ s, atol=1e-7)
        # vector fields with a bracket: Y = d_0, Z = x_0 d_1, [Y, Z] = d_1
        nabla_z = lambda y: y[0] * family.proj(y) @ derivative(y, 1)
        nabla_y = covariant(0)
        left = pi @ _diff(nabla_z, x, 0) - x[0] * pi @ _diff(nabla_y, x, 1) - pi @ derivative(x, 1)
        np.testing.assert_allclose(left, x[0] * family.curvature(x, 0, 1) @ s, atol=1e-7)
        # compatibility with the inner product
        other, other_derivative = _section(rng, family)
        for a in range(p):
            pairing = lambda y: np.array([[np.vdot(value(y), other(y))]])
            expected = np.vdot(pi @ derivative(x, a), other(x)) + np.vdot(s, pi @ other_derivative(x, a))
            np.testing.assert_allclose(_diff(pairing, x, a)[0, 0], expected, atol=1e-7)
        # the mixed returns of a unit vector form a positive Hermitian matrix in the two directions;
        # its antisymmetric part is one half of the curvature and its diagonal is a single positive return
        xi, _ = family.frame(x)
        psi = xi @ (lambda c: c / np.linalg.norm(c))(_complex(rng, r))
        tensor = np.array([[np.vdot(tangent[a] @ psi, tangent[b] @ psi) for b in range(p)] for a in range(p)])
        np.testing.assert_allclose(tensor, tensor.conj().T, atol=1e-12)
        assert np.linalg.eigvalsh(tensor).min() > -1e-12
        for a, b in itertools.combinations(range(p), 2):
            np.testing.assert_allclose(tensor[a, b] - tensor[b, a], np.vdot(psi, family.curvature(x, a, b) @ psi), atol=1e-12)


def test_frames_and_change_of_frame():
    rng = np.random.default_rng(204)
    for _ in range(12):
        n, r, p = rng.integers(3, 7), rng.integers(1, 4), 2
        r = min(r, n - 1)
        family = Family(rng, n, r, p)
        x = rng.uniform(-0.5, 0.5, size=p)
        xi, dxi = family.frame(x)
        pi = xi @ xi.conj().T
        connection = family.connection(x)
        for a in range(p):
            np.testing.assert_allclose(connection[a].conj().T, -connection[a], atol=1e-12)
        # nabla_a (Xi c) = Xi (d_a c + A_a c)
        c0, c1 = _complex(rng, r), _complex(rng, p, r)
        for a in range(p):
            ambient = dxi[a] @ (c0 + x @ c1) + xi @ c1[a]
            np.testing.assert_allclose(pi @ ambient, xi @ (c1[a] + connection[a] @ (c0 + x @ c1)), atol=1e-12)
        # the frame curvature represents the projector curvature
        frame_curvature = (_diff(lambda y: family.connection(y)[1], x, 0) - _diff(lambda y: family.connection(y)[0], x, 1)
                           + _comm(connection[0], connection[1]))
        np.testing.assert_allclose(frame_curvature, xi.conj().T @ family.curvature(x, 0, 1) @ xi, atol=1e-7)
        # change of frame by a position-dependent unitary matrix
        gens = [0.7 * _hermitian(rng, r) for _ in range(p)]

        def change(y):
            return expm(1j * sum(y[a] * gens[a] for a in range(p)))

        def changed_connection(y, a):
            frame, dframe = family.frame(y)
            q = change(y)
            dq = expm_frechet(1j * sum(y[b] * gens[b] for b in range(p)), 1j * gens[a], compute_expm=False)
            return (frame @ q).conj().T @ (dframe[a] @ q + frame @ dq), q, dq

        new = []
        for a in range(p):
            matrix, q, dq = changed_connection(x, a)
            np.testing.assert_allclose(matrix, q.conj().T @ connection[a] @ q + q.conj().T @ dq, atol=1e-11)
            new.append(matrix)
        new_curvature = (_diff(lambda y: changed_connection(y, 1)[0], x, 0) - _diff(lambda y: changed_connection(y, 0)[0], x, 1)
                         + _comm(new[0], new[1]))
        q = change(x)
        np.testing.assert_allclose(new_curvature, q.conj().T @ frame_curvature @ q, atol=1e-7)
        # the local frame built from one orthonormal basis and the projections near it
        base = np.zeros(p)
        xi_base, _ = family.frame(base)
        for _ in range(4):
            y = rng.uniform(-0.3, 0.3, size=p)
            pi_y = family.proj(y)
            gram = xi_base.conj().T @ pi_y @ xi_base
            assert np.linalg.eigvalsh(gram).min() > 1e-3
            local = pi_y @ xi_base @ np.linalg.inv(sqrtm(gram))
            np.testing.assert_allclose(local.conj().T @ local, np.eye(r), atol=1e-10)
            np.testing.assert_allclose(local @ local.conj().T, pi_y, atol=1e-10)


def _transport(family, path, velocity, t_end=1.0):
    """Propagator of ds/dt = [dPi/dt, Pi] s along a path in the parameter space."""
    n = family.n

    def rhs(t, flat):
        pi, dpi = family.dproj(path(t))
        rate = sum(v * d for v, d in zip(velocity(t), dpi))
        return (_comm(rate, pi) @ flat.reshape(n, n)).ravel()

    solution = solve_ivp(rhs, (0.0, t_end), np.eye(n, dtype=complex).ravel(), method="DOP853", rtol=1e-11, atol=1e-13)
    return solution.y[:, -1].reshape(n, n)


def test_transport_and_holonomy():
    rng = np.random.default_rng(205)
    for _ in range(5):
        n, r, p = rng.integers(3, 6), rng.integers(1, 3), 2
        family = Family(rng, n, r, p)
        centre, radius = rng.uniform(-0.2, 0.2, size=p), 0.35
        path = lambda t: centre + radius * np.array([np.cos(2 * np.pi * t) - 1, np.sin(2 * np.pi * t)])
        velocity = lambda t: 2 * np.pi * radius * np.array([-np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)])
        # the generator identity behind the intertwining
        pi, dpi = family.dproj(path(0.3))
        rate = sum(v * d for v, d in zip(velocity(0.3), dpi))
        np.testing.assert_allclose(_comm(_comm(rate, pi), pi), rate, atol=1e-12)
        np.testing.assert_allclose(_comm(rate, pi).conj().T, -_comm(rate, pi), atol=1e-12)
        # unitary propagator that carries the starting subspace onto the moving one
        for t_end in (0.37, 1.0):
            propagator = _transport(family, path, velocity, t_end)
            np.testing.assert_allclose(propagator.conj().T @ propagator, np.eye(n), atol=1e-9)
            np.testing.assert_allclose(propagator @ family.proj(path(0.0)) @ propagator.conj().T, family.proj(path(t_end)), atol=1e-9)
        xi0, _ = family.frame(path(0.0))
        holonomy = xi0.conj().T @ propagator @ xi0
        np.testing.assert_allclose(holonomy.conj().T @ holonomy, np.eye(r), atol=1e-9)
        np.testing.assert_allclose(propagator @ xi0, xi0 @ holonomy, atol=1e-9)        # the fibre returns to itself
        # the same matrix from the transport equation in the frame
        def frame_rhs(t, flat):
            matrices = family.connection(path(t))
            return (-sum(v * m for v, m in zip(velocity(t), matrices)) @ flat.reshape(r, r)).ravel()

        in_frame = solve_ivp(frame_rhs, (0.0, 1.0), np.eye(r, dtype=complex).ravel(), method="DOP853",
                             rtol=1e-11, atol=1e-13).y[:, -1].reshape(r, r)
        np.testing.assert_allclose(in_frame, holonomy, atol=1e-8)
        # transport is the limit of successive projections along the path
        steps = 3000
        carried = xi0
        for k in range(1, steps + 1):
            carried = family.proj(path(k / steps)) @ carried
        np.testing.assert_allclose(xi0.conj().T @ carried, holonomy, atol=5e-3)
        # a change of frame conjugates the holonomy by its value at the base point
        q0 = expm(1j * _hermitian(rng, r))
        np.testing.assert_allclose((xi0 @ q0).conj().T @ propagator @ (xi0 @ q0), q0.conj().T @ holonomy @ q0, atol=1e-10)
    # rank one: the holonomy of the boundary of a small square is exp(-integral of the curvature), not a root of unity
    family = Family(rng, 4, 1, 2)
    side = 0.4
    corners = [np.array(c) for c in ((0.0, 0.0), (side, 0.0), (side, side), (0.0, side), (0.0, 0.0))]
    total = np.eye(4, dtype=complex)
    for start, end in zip(corners[:-1], corners[1:]):
        path = lambda t, start=start, end=end: start + t * (end - start)
        velocity = lambda t, start=start, end=end: end - start
        total = _transport(family, path, velocity) @ total
    xi0, _ = family.frame(corners[0])
    phase = (xi0.conj().T @ total @ xi0)[0, 0]
    nodes, weights = np.polynomial.legendre.leggauss(14)
    nodes, weights = side * (nodes + 1) / 2, side * weights / 2
    flux = sum(wa * wb * np.trace(family.curvature([a, b], 0, 1)) for a, wa in zip(nodes, weights) for b, wb in zip(nodes, weights))
    np.testing.assert_allclose(phase, np.exp(-flux), atol=1e-8)
    assert abs(flux.real) < 1e-10 and 1e-3 < abs(flux.imag) < 1.0


def _wedge_trace(dpi, pi, indices):
    """Coefficient of Tr(Pi (dPi)^(wedge m)) on the given ordered tuple of m coordinate indices."""
    total = 0.0
    for permutation in itertools.permutations(range(len(indices))):
        sign = np.linalg.det(np.eye(len(indices))[list(permutation)])
        product = pi
        for position in permutation:
            product = product @ dpi[indices[position]]
        total += sign * np.trace(product)
    return total


def test_characteristic_forms_are_closed():
    rng = np.random.default_rng(206)
    # products of an odd number of derivatives are odd for the symmetry and have no trace
    family = Family(rng, 6, 2, 5, strength=0.35)
    x = rng.uniform(-0.4, 0.4, size=5)
    pi, dpi = family.dproj(x)
    symmetry = 2 * pi - np.eye(6)
    for a, b in itertools.permutations(range(5), 2):
        np.testing.assert_allclose(_comm(dpi[a] @ dpi[b], pi), 0, atol=1e-12)        # d Pi wedge d Pi commutes with Pi
    for triple in itertools.product(range(5), repeat=3):
        odd = dpi[triple[0]] @ dpi[triple[1]] @ dpi[triple[2]]
        np.testing.assert_allclose(symmetry @ odd @ symmetry, -odd, atol=1e-12)
        assert abs(np.trace(odd)) < 1e-12
    # degree two: d Tr(Pi dPi dPi) = 0 on every triple of coordinates
    two_form = lambda a, b: (lambda y: np.array([[np.trace(family.curvature(y, a, b))]]))
    for a, b, c in itertools.combinations(range(5), 3):
        terms = [_diff(two_form(b, c), x, a)[0, 0], _diff(two_form(c, a), x, b)[0, 0], _diff(two_form(a, b), x, c)[0, 0]]
        assert abs(sum(terms)) < 1e-7 * (1 + max(abs(t) for t in terms))
    # degree four: the powers of the curvature are Pi (dPi)^4, and d Tr(Pi (dPi)^4) = 0 in five parameters
    quad = (0, 1, 2, 3)
    curvature = {(a, b): family.curvature(x, a, b) for a, b in itertools.permutations(quad, 2)}
    squared = sum(np.linalg.det(np.eye(4)[list(perm)]) * np.trace(curvature[(quad[perm[0]], quad[perm[1]])] @ curvature[(quad[perm[2]], quad[perm[3]])])
                  for perm in itertools.permutations(range(4))) / 4
    np.testing.assert_allclose(squared, _wedge_trace(dpi, pi, quad), atol=1e-11)
    four_form = lambda indices: (lambda y: np.array([[_wedge_trace(family.dproj(y)[1], family.proj(y), indices)]]))
    terms = []
    for a in range(5):
        rest = tuple(i for i in range(5) if i != a)
        terms.append((-1) ** a * _diff(four_form(rest), x, a)[0, 0])
    assert max(abs(t) for t in terms) > 1e-4                                         # the single terms do not vanish
    assert abs(sum(terms)) < 1e-6 * (1 + max(abs(t) for t in terms))


SIGMA = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], dtype=complex)]


def _two_band(m):
    def h(kx, ky):
        return np.sin(kx) * SIGMA[0] + np.sin(ky) * SIGMA[1] + (m + np.cos(kx) + np.cos(ky)) * SIGMA[2]
    return h


def _lower_frame(matrix, rank):
    values, vectors = np.linalg.eigh(matrix)
    assert values[rank] - values[rank - 1] > 0.2                                     # the selected bands are separated
    return vectors[:, :rank]


def _lattice_flux(frames):
    """(1/2 pi) times the sum of the phases of the edge determinants around every cell of a periodic grid."""
    n1, n2 = len(frames), len(frames[0])
    edge = lambda v, w: np.linalg.det(frames[w[0] % n1][w[1] % n2].conj().T @ frames[v[0] % n1][v[1] % n2])
    total = 0.0
    for i in range(n1):
        for j in range(n2):
            corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1), (i, j)]
            total += np.angle(np.prod([edge(v, w) for v, w in zip(corners[:-1], corners[1:])]))
    return total / (2 * np.pi)


def _torus_flux(hamiltonian, rank, grid):
    """Lattice value and quadrature value of (i / 2 pi) times the integral of the trace of the curvature."""
    points = 2 * np.pi * np.arange(grid) / grid
    frames = [[_lower_frame(hamiltonian(kx, ky), rank) for ky in points] for kx in points]
    lattice = _lattice_flux(frames)
    proj = lambda kx, ky: (lambda f: f @ f.conj().T)(_lower_frame(hamiltonian(kx, ky), rank))
    step, quadrature = 1e-4, 0.0
    for kx in points:
        for ky in points:
            pi = proj(kx, ky)
            dx = (proj(kx + step, ky) - proj(kx - step, ky)) / (2 * step)
            dy = (proj(kx, ky + step) - proj(kx, ky - step)) / (2 * step)
            quadrature += np.trace(pi @ _comm(dx, dy))
    quadrature = (1j / (2 * np.pi)) * quadrature * (2 * np.pi / grid) ** 2
    return lattice, quadrature


def test_integer_flux():
    rng = np.random.default_rng(207)
    # the sphere of lines in C^2: curvature (i/2) sin(theta), flux -1, and +1 for the orthogonal line
    nodes, weights = np.polynomial.legendre.leggauss(20)
    thetas, weights = np.pi * (nodes + 1) / 2, np.pi * weights / 2
    for sign, expected in ((1, -1.0), (-1, 1.0)):
        flux = 0.0
        for theta, weight in zip(thetas, weights):
            phi = rng.uniform(0, 2 * np.pi)
            direction = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            d_theta = np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)])
            d_phi = np.array([-np.sin(theta) * np.sin(phi), np.sin(theta) * np.cos(phi), 0.0])
            dot = lambda vector: sum(c * s for c, s in zip(vector, SIGMA))
            pi = (np.eye(2) + sign * dot(direction)) / 2
            density = np.trace(pi @ _comm(sign * dot(d_theta) / 2, sign * dot(d_phi) / 2))
            np.testing.assert_allclose(density, sign * 0.5j * np.sin(theta), atol=1e-12)
            flux += weight * 2 * np.pi * density
        np.testing.assert_allclose(1j * flux / (2 * np.pi), expected, atol=1e-12)
    # rank two over a torus: the sum of edge determinants used in the proof is an integer and equals the integral
    mix, _ = np.linalg.qr(_complex(rng, 4, 4))
    coupling = 0.15 * _hermitian(rng, 4)
    results = {}
    for masses in ((1.0, 1.0), (1.0, -1.0), (1.0, 3.0), (3.0, 3.0)):
        first, second = _two_band(masses[0]), _two_band(masses[1])

        def hamiltonian(kx, ky, first=first, second=second):
            block = np.zeros((4, 4), dtype=complex)
            block[:2, :2], block[2:, 2:] = first(kx, ky), second(kx, ky)
            return mix @ (block + coupling) @ mix.conj().T

        lattice, quadrature = _torus_flux(hamiltonian, 2, 36)
        assert abs(lattice - round(lattice)) < 1e-9
        np.testing.assert_allclose(quadrature, lattice, atol=2e-4)
        results[masses] = round(lattice)
    assert abs(results[(1.0, 1.0)]) == 2 and abs(results[(1.0, 3.0)]) == 1
    assert results[(1.0, -1.0)] == 0 and results[(3.0, 3.0)] == 0


def test_second_chern_number_on_the_four_sphere():
    """Rank two over S^4: the first Chern form vanishes and ch_2 integrates to an integer of modulus one."""
    rng = np.random.default_rng(208)
    s0, (s1, s2, s3) = np.eye(2), SIGMA
    gamma = [np.kron(s1, s1), np.kron(s1, s2), np.kron(s1, s3), np.kron(s2, s0), np.kron(s3, s0)]
    for a, b in itertools.product(range(5), repeat=2):
        np.testing.assert_allclose(gamma[a] @ gamma[b] + gamma[b] @ gamma[a], 2 * (a == b) * np.eye(4), atol=1e-12)
    volume = 8 * np.pi**2 / 3
    values = []
    for _ in range(6):
        frame, _ = np.linalg.qr(rng.normal(size=(5, 5)))
        if np.linalg.det(frame) < 0:
            frame[:, 4] *= -1                                       # (n, e_1, ..., e_4) positively oriented
        point, tangents = frame[:, 0], frame[:, 1:].T
        dot = lambda vector: sum(c * g for c, g in zip(vector, gamma))
        pi = (np.eye(4) + dot(point)) / 2
        np.testing.assert_allclose(pi @ pi, pi, atol=1e-12)
        assert round(np.trace(pi).real) == 2
        dpi = [dot(tangent) / 2 for tangent in tangents]
        curvature = {(a, b): pi @ _comm(dpi[a], dpi[b]) @ pi for a, b in itertools.permutations(range(4), 2)}
        for key in curvature:
            assert abs(np.trace(curvature[key])) < 1e-12              # first Chern form
        density = 2 * (np.trace(curvature[(0, 1)] @ curvature[(2, 3)]) - np.trace(curvature[(0, 2)] @ curvature[(1, 3)])
                       + np.trace(curvature[(0, 3)] @ curvature[(1, 2)]))
        np.testing.assert_allclose(density, _wedge_trace(dpi, pi, (0, 1, 2, 3)), atol=1e-12)
        values.append(0.5 * (1j / (2 * np.pi)) ** 2 * density * volume)
    np.testing.assert_allclose(values, values[0], atol=1e-12)        # the density is the same at every point
    np.testing.assert_allclose(abs(values[0]), 1.0, atol=1e-12)
    assert abs(values[0].imag) < 1e-12


def _pair_plane(theta, s=0.4):
    pair = PairModel(theta)
    shell = pair.F2
    h = shell.conj().T @ pair.exchange_hamiltonian(s) @ shell
    plane = shell.conj().T @ pair.A
    vectors = {name: shell.conj().T @ getattr(pair, name) for name in ("g", "e", "m", "b", "dvec")}
    return h, plane @ plane.conj().T, vectors, pair.j - s, pair.u, pair.v


def test_turning_and_waiting_in_the_pair_model():
    rng = np.random.default_rng(209)
    eye = np.eye(6)
    for theta in (0.2, np.pi / 8, 0.6):
        h, p, vec, delta, u, v = _pair_plane(theta)
        b, d = vec["b"], vec["dvec"]
        np.testing.assert_allclose(p, np.outer(vec["g"], vec["g"].conj()) + np.outer(b, b.conj()), atol=1e-12)
        step = 1e-5
        d_theta = (_pair_plane(theta + step)[1] - _pair_plane(theta - step)[1]) / (2 * step)
        d_time = -1j * _comm(h, p)
        turning = (eye - p) @ d_theta @ p
        waiting = (eye - p) @ d_time @ p
        np.testing.assert_allclose(turning, 2 * np.outer(d, b.conj()), atol=1e-8)
        coupling = (eye - p) @ h @ p                                 # C of the Schur formula
        np.testing.assert_allclose(coupling, -delta * u * v * np.outer(d, b.conj()), atol=1e-12)
        np.testing.assert_allclose(waiting, -1j * coupling, atol=1e-12)
        curvature = turning.conj().T @ waiting - waiting.conj().T @ turning
        np.testing.assert_allclose(curvature, 4j * delta * u * v * np.outer(b, b.conj()), atol=1e-7)
        np.testing.assert_allclose(p @ _comm(d_theta, d_time) @ p, curvature, atol=1e-7)
        np.testing.assert_allclose(np.trace(curvature), 2j * delta * np.sin(4 * theta), atol=1e-7)
        # at a later time everything is conjugated by the evolution
        t = rng.uniform(0.3, 5.0)
        evolution = expm(-1j * t * h)
        later = lambda q: evolution @ q @ evolution.conj().T
        d_theta_later = (later(_pair_plane(theta + step)[1]) - later(_pair_plane(theta - step)[1])) / (2 * step)
        np.testing.assert_allclose(later(p) @ _comm(d_theta_later, -1j * _comm(h, later(p))) @ later(p), later(curvature), atol=1e-7)
        # one period: exact return, holonomy exp(2 pi i u^2), and the factor from the mean energy
        period = 2 * np.pi / delta
        np.testing.assert_allclose(expm(-1j * period * h) @ b, b, atol=1e-10)

        def rhs(time, flat):
            pi = expm(-1j * time * h) @ p @ expm(1j * time * h)
            return (_comm(-1j * _comm(h, pi), pi) @ flat.reshape(6, 6)).ravel()

        transport = solve_ivp(rhs, (0.0, period), eye.astype(complex).ravel(), method="DOP853",
                              rtol=1e-11, atol=1e-13).y[:, -1].reshape(6, 6)
        np.testing.assert_allclose(transport @ b, np.exp(2j * np.pi * u**2) * b, atol=1e-8)
        np.testing.assert_allclose(transport @ vec["g"], vec["g"], atol=1e-8)
        mean_energy = np.vdot(b, h @ b).real
        np.testing.assert_allclose(mean_energy, delta * u**2, atol=1e-12)
        np.testing.assert_allclose(np.exp(-1j * mean_energy * period) * np.exp(2j * np.pi * u**2), 1.0, atol=1e-10)
        if np.isclose(theta, np.pi / 8, rtol=0.0, atol=1e-12):
            np.testing.assert_allclose(transport @ b, -b, atol=1e-8)
    # flux over all settings and one period: quadrature of the trace, and the edge determinants of the proof
    delta = _pair_plane(0.1)[3]
    period = 2 * np.pi / delta
    nodes, weights = np.polynomial.legendre.leggauss(12)
    thetas, weights = np.pi / 4 * (nodes + 1) / 2, np.pi / 4 * weights / 2
    flux = 0.0
    for theta, weight in zip(thetas, weights):
        h, p, vec, delta, u, v = _pair_plane(theta)
        step = 1e-5
        d_theta = (_pair_plane(theta + step)[1] - _pair_plane(theta - step)[1]) / (2 * step)
        flux += weight * period * np.trace(p @ _comm(d_theta, -1j * _comm(h, p)) @ p)
    np.testing.assert_allclose(1j * flux / (2 * np.pi), -1.0, atol=1e-7)
    n_theta, n_time = 24, 48
    frames = []
    for i in range(n_theta + 1):
        h, p, vec, delta, u, v = _pair_plane(np.pi / 4 * i / n_theta)
        frames.append([expm(-1j * (period * k / n_time) * h) @ np.column_stack([vec["g"], vec["b"]]) for k in range(n_time)])
    total = 0.0
    edge = lambda v, w: np.linalg.det(frames[w[0]][w[1] % n_time].conj().T @ frames[v[0]][v[1] % n_time])
    for i in range(n_theta):
        for k in range(n_time):
            corners = [(i, k), (i + 1, k), (i + 1, k + 1), (i, k + 1), (i, k)]
            total += np.angle(np.prod([edge(v, w) for v, w in zip(corners[:-1], corners[1:])]))
    np.testing.assert_allclose(total / (2 * np.pi), -1.0, atol=1e-9)
    # the planes at the two end settings do not depend on the time waited, so the parameters close up to a sphere
    for theta in (0.0, np.pi / 4):
        h, p, *_ = _pair_plane(theta)
        np.testing.assert_allclose(_comm(h, p), 0, atol=1e-12)
    # for any Hamiltonian and projection the tangent map of waiting is -i times the coupling
    for _ in range(10):
        n, r = 6, rng.integers(1, 4)
        h = _hermitian(rng, n)
        basis, _ = np.linalg.qr(_complex(rng, n, n))
        p = basis[:, :r] @ basis[:, :r].conj().T
        moved = lambda t: expm(-1j * t * h) @ p @ expm(1j * t * h)
        rate = (moved(1e-5) - moved(-1e-5)) / 2e-5
        np.testing.assert_allclose((np.eye(n) - p) @ rate @ p, -1j * (np.eye(n) - p) @ h @ p, atol=1e-8)
