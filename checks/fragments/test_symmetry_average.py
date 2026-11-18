"""Checks for the fragment on symmetry averaging (manuscripts/consolidated/fragments/symmetry_average.tex).

The fragment carries TAUM's section "Compact symmetry averaging and invariant compression"
into the consolidated manuscript.  Its own labels are not in Universal_Model.tex yet, so
this file cites the TAUM labels of the statements it tests and the labels of the
consolidated manuscript that the fragment uses.  hbar = 1 and random numbers use fixed seeds.

Averages over a circle group and over SU(2) are computed by quadrature rules that are exact
for the matrix elements that occur.
"""
import itertools

import numpy as np
from scipy.linalg import expm, null_space

from lib.fock import dag
from lib.linalg import minimal_space
from pair.model import PairModel

LABELS = {
    "TAUM": ["thm:Haar-projector", "cor:finite-group-average", "eq:Haar-projector", "eq:Haar-properties",
             "eq:J-compatible-representation", "eq:Haar-J-commute", "eq:finite-group-projector"],
    "UM": ["prop:preparation-split", "prop:passive-two-spans", "prop:passive-opposite-sectors",
           "prop:shell-accuracy", "eq:pair-two-setting-ranges", "subsec:moving-pair"],
}


def _complex(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _hermitian(rng, n):
    raw = _complex(rng, n, n)
    return (raw + raw.conj().T) / 2


def _unitary(rng, n):
    return np.linalg.qr(_complex(rng, n, n))[0]


def _average(operators):
    return sum(operators) / len(operators)


def _fixed_projection(operators):
    """Projection onto the vectors fixed by every operator, from the joint null space of the operators minus one."""
    n = operators[0].shape[0]
    basis = null_space(np.vstack([op - np.eye(n) for op in operators]), rcond=1e-10)
    return basis @ basis.conj().T


def _is_projection(pi):
    np.testing.assert_allclose(pi @ pi, pi, atol=1e-10)
    np.testing.assert_allclose(pi, pi.conj().T, atol=1e-10)


def _permutations(n):
    return [np.eye(n)[list(p)] for p in itertools.permutations(range(n))]


def _circle(charges, mix):
    """phi -> exp(-i phi N) with N = mix diag(charges) mix^*, and a rule that averages it exactly."""
    generator = mix @ np.diag(charges) @ mix.conj().T
    nodes = 2 * int(max(abs(c) for c in charges)) + 3
    kappa = lambda phi: mix @ np.diag(np.exp(-1j * phi * np.asarray(charges))) @ mix.conj().T
    return generator, kappa, [2 * np.pi * k / nodes for k in range(nodes)]


def _su2(alpha, beta, gamma):
    sz, sy = np.diag([1.0, -1.0]), np.array([[0, -1j], [1j, 0]])
    return expm(-0.5j * alpha * sz) @ expm(-0.5j * beta * sy) @ expm(-0.5j * gamma * sz)


def _su2_average(power, order):
    """Haar average of u tensored with itself `power` times, exact for that power."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total, norm = 0.0, 0.0
    for x, w in zip(nodes, weights):
        for alpha in 2 * np.pi * np.arange(2 * order) / (2 * order):
            for gamma in 4 * np.pi * np.arange(4 * order) / (4 * order):
                u = _su2(alpha, np.arccos(x), gamma)
                rep = u
                for _ in range(power - 1):
                    rep = np.kron(rep, u)
                total, norm = total + w * rep, norm + w
    return total / norm


def test_average_is_the_projection_onto_invariant_vectors():
    rng = np.random.default_rng(301)
    # finite groups: permutations of three and four points, on the points and on ordered pairs, in a random basis
    for n in (3, 4):
        for tensor in (False, True):
            group = [np.kron(p, p) if tensor else p for p in _permutations(n)]
            mix = _unitary(rng, group[0].shape[0])
            group = [mix @ g @ mix.conj().T for g in group]
            pi = _average(group)
            _is_projection(pi)
            np.testing.assert_allclose(pi, _fixed_projection(group), atol=1e-10)
            assert round(np.trace(pi).real) == (2 if tensor else 1)        # orbits on pairs, and on points
            for g in group:
                np.testing.assert_allclose(g @ pi, pi, atol=1e-10)
                np.testing.assert_allclose(pi @ g, pi, atol=1e-10)
            # a closed invariant subspace is mapped into itself
            start = _complex(rng, group[0].shape[0])
            orbit = np.column_stack([g @ start for g in group])
            basis = np.linalg.svd(orbit, full_matrices=False)[0][:, :np.linalg.matrix_rank(orbit, tol=1e-9)]
            inside = basis @ basis.conj().T
            assert basis.shape[1] < group[0].shape[0] or not tensor     # on pairs the orbit spans a proper invariant subspace
            for g in group:
                np.testing.assert_allclose(inside @ g @ inside, g @ inside, atol=1e-9)
            np.testing.assert_allclose(inside @ pi @ inside, pi @ inside, atol=1e-9)
    # a circle group
    charges = [0, 0, 1, -2, 3, 5, -5, 0]
    generator, kappa, angles = _circle(charges, _unitary(rng, len(charges)))
    pi = _average([kappa(phi) for phi in angles])
    _is_projection(pi)
    assert round(np.trace(pi).real) == 3
    np.testing.assert_allclose(generator @ pi, 0, atol=1e-10)             # the spectral projection of the generator at zero
    for phi in rng.uniform(0, 2 * np.pi, size=4):
        np.testing.assert_allclose(kappa(phi) @ pi, pi, atol=1e-10)
        np.testing.assert_allclose(pi @ kappa(phi), pi, atol=1e-10)
    # twisted averages give the other eigenspaces of the generator
    for charge in (1, -2, 5, 4):
        twisted = _average([np.exp(1j * charge * phi) * kappa(phi) for phi in angles])
        _is_projection(twisted)
        assert round(np.trace(twisted).real) == charges.count(charge)
        np.testing.assert_allclose(generator @ twisted, charge * twisted, atol=1e-9)
    # a compact group that is not abelian: SU(2) on two, three and four spins
    for power, order, rank in ((2, 3, 1), (3, 4, 0), (4, 5, 2)):
        pi = _su2_average(power, order)
        _is_projection(pi)
        assert round(np.trace(pi).real) == rank
        for _ in range(3):
            u = _su2(*rng.uniform(0, 2 * np.pi, size=3))
            rep = u
            for _ in range(power - 1):
                rep = np.kron(rep, u)
            np.testing.assert_allclose(rep @ pi, pi, atol=1e-10)
            np.testing.assert_allclose(pi @ rep, pi, atol=1e-10)
    singlet = np.array([0, 1, -1, 0]) / np.sqrt(2)
    np.testing.assert_allclose(_su2_average(2, 3), np.outer(singlet, singlet), atol=1e-10)


def test_subgroups_and_refinement():
    rng = np.random.default_rng(302)
    charges = list(range(-6, 7)) + [0, 4, 8, -8]
    generator, kappa, angles = _circle(charges, _unitary(rng, len(charges)))
    whole = _average([kappa(phi) for phi in angles])
    previous = np.eye(len(charges))
    for n in (1, 2, 4, 8, 16):
        cyclic = _average([kappa(2 * np.pi * k / n) for k in range(n)])
        _is_projection(cyclic)
        assert round(np.trace(cyclic).real) == sum(c % n == 0 for c in charges)
        np.testing.assert_allclose(whole @ cyclic, whole, atol=1e-10)      # the smaller group fixes at least as much
        np.testing.assert_allclose(cyclic @ whole, whole, atol=1e-10)
        np.testing.assert_allclose(previous @ cyclic, cyclic, atol=1e-10)  # the chain of projections decreases
        previous = cyclic
        assert np.allclose(cyclic, whole) == (n > max(abs(c) for c in charges))
    # strong convergence without norm convergence, on a long truncation of a space with all charges
    count = 4097
    values = np.arange(count) - count // 2
    vector = 1.0 / (1.0 + np.abs(values)) ** 1.5
    vector /= np.linalg.norm(vector)
    errors = []
    for k in range(1, 12):
        n = 2**k
        fixed = (values % n == 0) & (values != 0)                        # charges fixed by Z_n and not by the circle
        assert fixed.any()                                               # the difference of projections has norm one
        errors.append(np.linalg.norm(vector[fixed]))
    assert all(a > b for a, b in zip(errors[:-1], errors[1:])) and errors[-1] < 1e-4


def _factorial(k):
    out = 1
    for i in range(2, k + 1):
        out *= i
    return out


def _cubic_rotations():
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            matrix = np.zeros((3, 3))
            for row, (col, sign) in enumerate(zip(perm, signs)):
                matrix[row, col] = sign
            if np.isclose(np.linalg.det(matrix), 1.0):
                out.append(matrix)
    return out


def _monomials(degree):
    return [(a, b, degree - a - b) for a in range(degree + 1) for b in range(degree + 1 - a)]


def _substitution(matrix, degree):
    """Action of a signed permutation matrix on homogeneous polynomials, p(x) -> p(matrix^T x), in the monomial basis."""
    basis = _monomials(degree)
    index = {m: i for i, m in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)))
    for mono in basis:
        exponents, sign = [0, 0, 0], 1.0
        for i in range(3):                                   # x_i is replaced by sum_j matrix[j, i] x_j
            j = int(np.nonzero(matrix[:, i])[0][0])
            exponents[j] += mono[i]
            sign *= matrix[j, i] ** mono[i]
        out[index[tuple(exponents)], index[mono]] = sign
    return out


def _rotation_generators(degree):
    """The derivations x_i d_j - x_j d_i on homogeneous polynomials, in the monomial basis."""
    basis = _monomials(degree)
    index = {m: i for i, m in enumerate(basis)}
    out = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        matrix = np.zeros((len(basis), len(basis)))
        for mono in basis:
            for a, b, sign in ((i, j, 1.0), (j, i, -1.0)):   # sign * x_a d_b
                if mono[b] > 0:
                    target = list(mono)
                    target[b] -= 1
                    target[a] += 1
                    matrix[index[tuple(target)], index[mono]] += sign * mono[b]
        out.append(matrix)
    return out


def test_cubic_rotations_fix_more_than_all_rotations():
    group = _cubic_rotations()
    assert len(group) == 24
    products = {tuple(np.rint(a @ b).astype(int).ravel()) for a in group for b in group}
    assert products == {tuple(g.astype(int).ravel()) for g in group}        # closed under multiplication
    cubic_counts, rotation_counts = [], []
    for degree in range(7):
        basis = _monomials(degree)
        scale = np.array([np.sqrt(float(np.prod([_factorial(k) for k in m]))) for m in basis])
        action = [_substitution(g, degree) for g in group]
        unitary = [a * scale[:, None] / scale[None, :] for a in action]     # orthonormal for the Fischer inner product
        for u in unitary:
            np.testing.assert_allclose(u.T @ u, np.eye(len(basis)), atol=1e-12)
        np.testing.assert_allclose(_substitution(group[3] @ group[7], degree), action[3] @ action[7], atol=1e-12)
        pi = _average(unitary)
        _is_projection(pi)
        cubic_counts.append(round(np.trace(pi)))
        generators = _rotation_generators(degree)
        fixed = null_space(np.vstack(generators), rcond=1e-10)
        rotation_counts.append(fixed.shape[1])
        if degree:
            # a quarter turn generated by one derivation is one of the 24 substitutions
            quarter = expm(np.pi / 2 * generators[0])
            assert any(np.allclose(quarter, a, atol=1e-9) for a in action)
        # everything fixed by all rotations is fixed by the cube group
        for column in fixed.T:
            np.testing.assert_allclose(_average(action) @ column, column, atol=1e-10)
        if degree == 4:
            index = {m: i for i, m in enumerate(basis)}
            quartic, mixed, radial = np.zeros(len(basis)), np.zeros(len(basis)), np.zeros(len(basis))
            for m in ((4, 0, 0), (0, 4, 0), (0, 0, 4)):
                quartic[index[m]] = 1.0
                radial[index[m]] = 1.0
            for m in ((2, 2, 0), (0, 2, 2), (2, 0, 2)):
                mixed[index[m]] = 1.0
                radial[index[m]] = 2.0
            mean = _average(action)
            for vector in (quartic, mixed):
                np.testing.assert_allclose(mean @ vector, vector, atol=1e-12)
            assert np.linalg.matrix_rank(np.column_stack([quartic, mixed])) == 2 == round(np.trace(mean))
            assert fixed.shape[1] == 1
            np.testing.assert_allclose(abs(np.vdot(fixed[:, 0], radial)) / np.linalg.norm(radial), 1.0, atol=1e-10)
    assert cubic_counts == [1, 0, 1, 0, 2, 0, 3]
    assert rotation_counts == [1, 0, 1, 0, 1, 0, 1]


def test_sectors_and_covariant_maps():
    rng = np.random.default_rng(303)
    perms = _permutations(3)
    for _ in range(8):
        mix_in, mix_out = _unitary(rng, 6), _unitary(rng, 9)
        rep_in = [mix_in @ np.kron(np.eye(2), p) @ mix_in.conj().T for p in perms]
        rep_out = [mix_out @ np.kron(p, p) @ mix_out.conj().T for p in perms]
        both = [np.block([[a, np.zeros((6, 9))], [np.zeros((9, 6)), b]]) for a, b in zip(rep_in, rep_out)]
        symmetry = np.diag(np.r_[np.ones(6), -np.ones(9)])
        pi, pi_in, pi_out = _average(both), _average(rep_in), _average(rep_out)
        np.testing.assert_allclose(pi @ symmetry, symmetry @ pi, atol=1e-10)
        np.testing.assert_allclose(pi, np.block([[pi_in, np.zeros((6, 9))], [np.zeros((9, 6)), pi_out]]), atol=1e-10)
        np.testing.assert_allclose(symmetry @ pi.conj().T @ symmetry, pi, atol=1e-10)      # self-adjoint for the form
        # a covariant map, obtained by averaging an arbitrary one, with a kernel
        raw = _complex(rng, 9, 6)
        covariant = _average([b @ raw @ a.conj().T for a, b in zip(rep_in, rep_out)])
        if rng.random() < 0.5:
            covariant = covariant @ (np.eye(6) - pi_in)       # still covariant, with the invariant input vectors in its kernel
        for a, b in zip(rep_in, rep_out):
            np.testing.assert_allclose(covariant @ a, b @ covariant, atol=1e-9)
            np.testing.assert_allclose(covariant.conj().T @ b, a @ covariant.conj().T, atol=1e-9)
        np.testing.assert_allclose(covariant @ pi_in, pi_out @ covariant, atol=1e-9)
        np.testing.assert_allclose(covariant.conj().T @ pi_out, pi_in @ covariant.conj().T, atol=1e-9)
        left, singular, right = np.linalg.svd(covariant)
        rank = int(np.sum(singular > 1e-9))
        support_in = right[:rank].conj().T @ right[:rank]                    # onto the orthogonal complement of the kernel
        support_out = left[:, :rank] @ left[:, :rank].conj().T                # onto the range
        for a, b in zip(rep_in, rep_out):
            np.testing.assert_allclose(support_in @ a, a @ support_in, atol=1e-8)
            np.testing.assert_allclose(support_out @ b, b @ support_out, atol=1e-8)
        np.testing.assert_allclose(support_in @ pi_in, pi_in @ support_in, atol=1e-8)
        np.testing.assert_allclose(support_out @ pi_out, pi_out @ support_out, atol=1e-8)
        # a map that is not covariant need not respect the invariant vectors
        assert np.linalg.norm(raw @ pi_in - pi_out @ raw) > 1e-3


def test_invariant_port_and_compressed_response():
    rng = np.random.default_rng(304)
    perms = _permutations(4)
    for _ in range(6):
        n = 16
        mix = _unitary(rng, n)
        group = [mix @ np.kron(p, p) @ mix.conj().T for p in perms]
        pi = _average(group)
        rank = round(np.trace(pi).real)
        assert rank == 2
        raw = _hermitian(rng, n)
        h = _average([g @ raw @ g.conj().T for g in group])                       # commutes with the representation
        for g in group:
            np.testing.assert_allclose(g @ h, h @ g, atol=1e-9)
        np.testing.assert_allclose(pi @ h, h @ pi, atol=1e-9)
        port = pi @ _complex(rng, n, 2)
        reached = minimal_space(h, port)
        np.testing.assert_allclose(pi @ reached, reached, atol=1e-8)                # the reached space lies in the invariant space
        basis = null_space(np.eye(n) - pi, rcond=1e-9)
        h_fixed, port_fixed = basis.conj().T @ h @ basis, basis.conj().T @ port
        outside = np.eye(n) - pi
        broken = h + outside @ _hermitian(rng, n) @ outside                         # changes H off the invariant space and breaks the symmetry
        assert max(np.linalg.norm(g @ broken - broken @ g) for g in group) > 1e-3
        for z in (0.3 + 0.4j, -1.1 - 0.2j):
            response = port.conj().T @ np.linalg.inv(h - z * np.eye(n)) @ port
            np.testing.assert_allclose(response, port_fixed.conj().T @ np.linalg.inv(h_fixed - z * np.eye(rank)) @ port_fixed, atol=1e-9)
            np.testing.assert_allclose(response, port.conj().T @ np.linalg.inv(broken - z * np.eye(n)) @ port, atol=1e-9)
        for t in rng.uniform(-5, 5, size=3):
            np.testing.assert_allclose(outside @ expm(-1j * t * h) @ port, 0, atol=1e-9)
        # a port that is not invariant reaches outside
        generic = _complex(rng, n, 1)
        assert np.linalg.norm(outside @ minimal_space(h, generic)) > 1e-3
    # parity: the group of two elements, and its even space
    parity = np.diag([1.0, -1.0, 1.0, -1.0, -1.0, 1.0])
    mix = _unitary(rng, 6)
    parity = mix @ parity @ mix.conj().T
    even = _average([np.eye(6), parity])
    _is_projection(even)
    np.testing.assert_allclose(parity @ even, even, atol=1e-12)
    assert round(np.trace(even).real) == 3


def _two_particle(pair, one_body):
    """Action on the two-particle space of a unitary change of the four modes c1, c2, d1, d2."""
    create = [dag(x) for x in (pair.c[0], pair.c[1], pair.d[0], pair.d[1])]
    new = [sum(one_body[b, a] * create[b] for b in range(4)) for a in range(4)]
    out = np.zeros((16, 16), dtype=complex)
    for a, b in itertools.combinations(range(4), 2):
        out += np.outer(new[a] @ new[b] @ pair.vac, (create[a] @ create[b] @ pair.vac).conj())
    return pair.F2.conj().T @ out @ pair.F2


def _orbital_rotation(phi):
    c, s = np.cos(phi), np.sin(phi)
    out = np.zeros((4, 4))
    for i in (0, 1):
        out[i, i], out[i + 2, i], out[i, i + 2], out[i + 2, i + 2] = c, s, -s, c
    return out


def test_two_symmetries_of_the_pair_model():
    theta = 0.3
    pair = PairModel(theta)
    shell = pair.F2
    vec = {name: shell.conj().T @ getattr(pair, name) for name in ("g", "e", "m", "b")}
    silent = shell.conj().T @ pair.W
    line = lambda v: np.outer(v, v.conj())
    plane = lambda angle: (lambda a: a @ a.conj().T)(shell.conj().T @ PairModel(angle).A)
    hamiltonian = lambda s: shell.conj().T @ pair.exchange_hamiltonian(s) @ shell
    eye = np.eye(6)
    # exchange of the two orbitals
    one_body = np.zeros((4, 4))
    one_body[2, 0] = one_body[0, 2] = one_body[3, 1] = one_body[1, 3] = 1.0
    exchange = _two_particle(pair, one_body)
    np.testing.assert_allclose(exchange @ exchange, eye, atol=1e-12)
    for s in (0.0, 0.4, 1.0):
        np.testing.assert_allclose(exchange @ hamiltonian(s), hamiltonian(s) @ exchange, atol=1e-12)
    for name, sign in (("g", 1), ("m", 1), ("e", -1)):
        np.testing.assert_allclose(exchange @ vec[name], sign * vec[name], atol=1e-12)
    np.testing.assert_allclose(exchange @ silent, -silent, atol=1e-12)
    pi_exchange = (eye + exchange) / 2
    np.testing.assert_allclose(pi_exchange, line(vec["g"]) + line(vec["m"]), atol=1e-12)
    np.testing.assert_allclose(pi_exchange, plane(np.pi / 4), atol=1e-12)
    for s in (0.0, 0.4):                                      # the orbit of the invariant preparation stays in the invariant plane
        reached = minimal_space(hamiltonian(s), shell.conj().T @ PairModel(np.pi / 4).A)
        np.testing.assert_allclose(reached @ reached.conj().T, pi_exchange, atol=1e-9)
        assert minimal_space(hamiltonian(s), shell.conj().T @ pair.A).shape[1] == 3
    np.testing.assert_allclose(pi_exchange @ plane(0.0), plane(0.0) @ pi_exchange, atol=1e-12)
    assert np.linalg.norm(pi_exchange - plane(0.0)) > 0.5
    for angle in (0.1, theta, 0.7):
        assert np.linalg.norm(pi_exchange @ plane(angle) - plane(angle) @ pi_exchange) > 1e-2
    # rotation of the two orbitals
    kappa = lambda phi: _two_particle(pair, _orbital_rotation(phi))
    np.testing.assert_allclose(kappa(0.3) @ kappa(0.5), kappa(0.8), atol=1e-12)
    np.testing.assert_allclose(kappa(theta) @ (shell.conj().T @ PairModel(0.0).A), shell.conj().T @ pair.A, atol=1e-12)
    np.testing.assert_allclose(kappa(theta) @ vec["e"], vec["b"], atol=1e-12)
    np.testing.assert_allclose(kappa(0.4) @ vec["g"], vec["g"], atol=1e-12)
    np.testing.assert_allclose(kappa(0.4) @ silent, silent, atol=1e-12)
    np.testing.assert_allclose(kappa(0.4) @ vec["e"], np.cos(0.8) * vec["e"] + np.sin(0.8) * vec["m"], atol=1e-12)
    cyclic = lambda n: _average([kappa(2 * np.pi * k / n) for k in range(n)])
    pi_rotation = cyclic(64)
    np.testing.assert_allclose(pi_rotation, line(vec["g"]) + silent @ silent.conj().T, atol=1e-12)
    assert round(np.trace(pi_rotation).real) == 4
    np.testing.assert_allclose(pi_rotation @ plane(theta), line(vec["g"]), atol=1e-12)
    for s, commutes in ((0.0, False), (0.4, False), (1.0, True)):
        assert np.allclose(kappa(0.3) @ hamiltonian(s), hamiltonian(s) @ kappa(0.3)) == commutes
    np.testing.assert_allclose(hamiltonian(1.0), -2 * line(vec["g"]), atol=1e-12)
    # finite subgroups of the rotations
    np.testing.assert_allclose(kappa(np.pi), eye, atol=1e-12)
    np.testing.assert_allclose(cyclic(2), eye, atol=1e-12)
    for n in (3, 4, 5, 7):
        np.testing.assert_allclose(cyclic(n), pi_rotation, atol=1e-12)
    # the two units of charge of the pair (e, m)
    for charge, rank in ((0, 4), (2, 1), (-2, 1), (1, 0), (-1, 0)):
        twisted = _average([np.exp(1j * charge * 2 * np.pi * k / 64) * kappa(2 * np.pi * k / 64) for k in range(64)])
        assert round(np.trace(twisted).real) == rank
    # the exchange is a reflection; together with the rotations it fixes only the ground
    np.testing.assert_allclose(exchange @ kappa(0.3) @ exchange, kappa(-0.3), atol=1e-12)
    orthogonal = [kappa(2 * np.pi * k / 16) for k in range(16)] + [exchange @ kappa(2 * np.pi * k / 16) for k in range(16)]
    np.testing.assert_allclose(_average(orthogonal), line(vec["g"]), atol=1e-12)
    np.testing.assert_allclose(pi_exchange @ pi_rotation, line(vec["g"]), atol=1e-12)
