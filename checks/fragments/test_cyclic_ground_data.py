"""Checks for the clauses added to Section 5: pointed representations and vector states,
the GNS intertwiner, the space generated from a ground and its complement, the
independence of uniqueness, cyclicity and a gap, and purity and faithfulness."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import null_space

from lib.fock import dag
from lib.linalg import rank

LABELS = {
    "LG": ["def:pointed-representation", "def:vector-state", "prop:vector-state-positive",
           "def:ground-datum", "rem:ground-properties-independent", "prop:cyclic-is-gns",
           "def:pure-faithful-state", "def:physical-cyclic-subspace",
           "prop:physical-centered-cyclic", "def:dark-complement", "prop:no-dark-sector",
           "q:luminous-purity", "choice:determination-algebra",
           "prop:covariance-content-inner-product"],
    "UM": ["sec:determination", "prop:cyclic-determination-content",
           "eq:determination-covariance", "part:problem-origin", "sec:ground-correspondence",
           "eq:ground-return-energy"],
}

RNG = np.random.default_rng(5101)


def random_complex(*shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def unit(vector):
    return vector / np.linalg.norm(vector)


def projection(vector):
    return np.outer(vector, vector.conj())


def span_projection(vectors, tol=1e-10):
    """Orthogonal projection onto the span of the columns."""
    left, singular, _ = np.linalg.svd(vectors, full_matrices=False)
    basis = left[:, singular > tol]
    return basis @ dag(basis)


def generated_space(basis_of_algebra, y):
    """Projection onto the closed span of X y, X running over a basis of the algebra."""
    return span_projection(np.column_stack([X @ y for X in basis_of_algebra]))


def commutant_dimension(basis_of_algebra):
    """Dimension of the set of T commuting with every element of the algebra."""
    n = basis_of_algebra[0].shape[0]
    rows = [np.kron(X, np.eye(n)) - np.kron(np.eye(n), X.T) for X in basis_of_algebra]
    return null_space(np.vstack(rows)).shape[1]


def is_separating(basis_of_algebra, y):
    """y is separating when X -> X y is injective on the algebra."""
    images = np.column_stack([X @ y for X in basis_of_algebra])
    return rank(images) == len(basis_of_algebra)


def matrix_units(n):
    units = []
    for i in range(n):
        for j in range(n):
            E = np.zeros((n, n), complex)
            E[i, j] = 1
            units.append(E)
    return units


def diagonal_units(n):
    return [np.diag(np.eye(n)[k]).astype(complex) for k in range(n)]


def cell_units(cells, n):
    units = []
    for cell in cells:
        d = np.zeros(n)
        d[list(cell)] = 1
        units.append(np.diag(d).astype(complex))
    return units


def path_hamiltonian(n):
    H = np.zeros((n, n))
    for k in range(n - 1):
        H[k, k] += 0.5
        H[k + 1, k + 1] += 0.5
        H[k, k + 1] -= 0.5
        H[k + 1, k] -= 0.5
    return H


FINE_GROUND = np.ones(5) / np.sqrt(5)
CELLS = ((0,), (1, 2, 3, 4))


def test_pointed_representation_vector_state_and_abstract_faithfulness():
    # The algebra M_2 + C, represented on C^2 (x) C^2 by X + c -> X (x) I. The summand C is
    # in the kernel of the representation.
    def pi(X, c):
        return np.kron(X, np.eye(2))

    assert_allclose(pi(np.eye(2), 1.0), np.eye(4))
    X1, X2 = random_complex(2, 2), random_complex(2, 2)
    assert_allclose(pi(X1 @ X2, 0.3 * 0.7), pi(X1, 0.3) @ pi(X2, 0.7), atol=1e-12)
    assert_allclose(pi(dag(X1), 0.0), dag(pi(X1, 0.0)), atol=1e-12)

    y = np.array([1, 0, 0, 1]) / np.sqrt(2)
    image_basis = [np.kron(E, np.eye(2)) for E in matrix_units(2)]

    def state(X, c):
        return np.vdot(y, pi(X, c) @ y)

    for _ in range(5):
        X, c = random_complex(2, 2), complex(*RNG.normal(size=2))
        Y, d = random_complex(2, 2), complex(*RNG.normal(size=2))
        # linear, positive and normalized
        assert_allclose(state(X + 2j * Y, c + 2j * d), state(X, c) + 2j * state(Y, d), atol=1e-12)
        value = state(dag(X) @ X, np.conj(c) * c)
        assert abs(value.imag) < 1e-12 and value.real >= -1e-12
        assert_allclose(value, np.linalg.norm(pi(X, c) @ y) ** 2, atol=1e-12)
        # the vector state is the trace state, and a phase of y changes nothing
        assert_allclose(state(X, c), np.trace(X) / 2, atol=1e-12)
        assert_allclose(np.vdot(np.exp(0.8j) * y, pi(X, c) @ (np.exp(0.8j) * y)), state(X, c), atol=1e-12)
    assert_allclose(state(np.eye(2), 1.0), 1, atol=1e-12)

    # y is cyclic and separating for the image, so the state of the image is faithful.
    assert_allclose(generated_space(image_basis, y), np.eye(4), atol=1e-12)
    assert is_separating(image_basis, y)
    # Pulled back to M_2 + C it is not faithful: 0 + 1 is nonzero and has value zero.
    assert_allclose(state(np.zeros((2, 2)), 1.0), 0, atol=1e-12)
    # The state is not pure: the commutant I (x) M_2 is four-dimensional, and the trace state
    # is the equal mixture of the two states X -> X_00 and X -> X_11.
    assert commutant_dimension(image_basis) == 4
    X = random_complex(2, 2)
    assert_allclose(state(X, 0.0), (X[0, 0] + X[1, 1]) / 2, atol=1e-12)

    # the space generated from other vectors
    product = np.array([1, 0, 0, 0], complex)
    P = generated_space(image_basis, product)
    assert rank(P) == 2
    assert_allclose(P @ product, product, atol=1e-12)
    assert_allclose(generated_space(image_basis, -3.7j * product), P, atol=1e-12)
    assert rank(np.column_stack([B @ np.zeros(4) for B in image_basis])) == 0


def gns(basis_of_algebra, y, tol=1e-10):
    """GNS data of the vector state on the algebra spanned by `basis_of_algebra`.

    Returns coordinates of the classes [X_i] in an orthonormal basis of the GNS space,
    and the matrices of left multiplication by each basis element."""
    m = len(basis_of_algebra)
    gram = np.array([[np.vdot(y, dag(Xi) @ Xj @ y) for Xj in basis_of_algebra] for Xi in basis_of_algebra])
    values, vectors = np.linalg.eigh(gram)
    keep = values > tol
    # classes [X_i] -> coordinates c_i with <c_i, c_j> = gram_ij
    coordinates = (np.diag(np.sqrt(values[keep])) @ dag(vectors[:, keep]))
    # products X_a X_j expanded in the basis (the algebras below are closed under products)
    flat = np.column_stack([B.reshape(-1) for B in basis_of_algebra])
    left = []
    for Xa in basis_of_algebra:
        products = np.column_stack([np.linalg.lstsq(flat, (Xa @ Xj).reshape(-1), rcond=None)[0] for Xj in basis_of_algebra])
        assert_allclose(flat @ products, np.column_stack([(Xa @ Xj).reshape(-1) for Xj in basis_of_algebra]), atol=1e-10)
        # [X_a X_j] = sum_k products[k, j] [X_k]
        left.append(coordinates @ products @ np.linalg.pinv(coordinates))
    assert coordinates.shape == (int(keep.sum()), m)
    return coordinates, left


def test_gns_intertwiner_unitary_and_unique():
    for algebra, expected_dimension in ((diagonal_units(5), 5), (cell_units(CELLS, 5), 2), (matrix_units(5), 5)):
        y = FINE_GROUND
        coordinates, left = gns(algebra, y)
        assert coordinates.shape[0] == expected_dimension
        # U [X] = X y is an isometry from the GNS space onto K
        vectors = np.column_stack([X @ y for X in algebra])
        U = vectors @ np.linalg.pinv(coordinates)
        assert_allclose(U @ coordinates, vectors, atol=1e-10)
        assert_allclose(dag(U) @ U, np.eye(expected_dimension), atol=1e-10)
        assert_allclose(U @ dag(U), generated_space(algebra, y), atol=1e-10)
        # U carries the class of the identity to y and intertwines the two actions
        flat = np.column_stack([B.reshape(-1) for B in algebra])
        unit_class = coordinates @ np.linalg.lstsq(flat, np.eye(5).reshape(-1).astype(complex), rcond=None)[0]
        assert_allclose(U @ unit_class, y, atol=1e-10)
        for Xa, La in zip(algebra, left):
            assert_allclose(U @ La, Xa @ U, atol=1e-9)

    # Uniqueness: an operator that commutes with the algebra and fixes the cyclic vector is I.
    # For the site algebra the commutant is diagonal, and diag(w) y = y has the one solution
    # w = (1, ..., 1) because every entry of y is nonzero.
    columns = np.column_stack([E @ FINE_GROUND for E in diagonal_units(5)])
    assert rank(columns) == 5
    assert_allclose(np.linalg.solve(columns, FINE_GROUND), np.ones(5), atol=1e-12)
    # For M_2 (x) I with the maximally entangled vector the commutant is I (x) M_2, and
    # (I (x) w) y = y is a linear system in w with the one solution w = I.
    y = np.array([1, 0, 0, 1]) / np.sqrt(2)
    columns = np.column_stack([np.kron(np.eye(2), E) @ y for E in matrix_units(2)])
    assert rank(columns) == 4
    w = np.linalg.solve(columns, y.astype(complex))
    assert_allclose(w.reshape(2, 2), np.eye(2), atol=1e-12)


def test_variance_generated_space_and_its_complement():
    n = 6
    y = unit(random_complex(n))
    Q = np.eye(n) - projection(y)
    for _ in range(10):
        X = random_complex(n, n)
        omega = lambda T: np.vdot(y, T @ y)
        variance = omega(dag(X) @ X) - abs(omega(X)) ** 2
        assert_allclose(variance, np.linalg.norm(Q @ X @ y) ** 2, atol=1e-10)
    # zero variance exactly on the ray
    X = projection(y) * 2.5 + Q @ random_complex(n, n) @ Q
    assert_allclose(np.vdot(y, dag(X) @ X @ y) - abs(np.vdot(y, X @ y)) ** 2, 0, atol=1e-10)

    # The unequal-cell path of Section 4: the cell algebra generates D = Ran A from the fine
    # ground, and the complement of the generated space is the defect ker A*.
    A = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]]) / 2
    P_D = span_projection(A)
    cells = cell_units(CELLS, 5)
    P_K = generated_space(cells, FINE_GROUND)
    assert_allclose(P_K, P_D, atol=1e-12)
    defect = null_space(A.T)
    assert defect.shape[1] == 3
    assert_allclose(np.eye(5) - P_K, defect @ dag(defect), atol=1e-12)
    for X in cells:
        assert_allclose(X @ P_K, P_K @ X, atol=1e-12)
        assert_allclose(dag(defect) @ X @ FINE_GROUND, 0, atol=1e-12)
    # centered vectors lie in K and are orthogonal to y
    for X in cells:
        centered = X @ FINE_GROUND - np.vdot(FINE_GROUND, X @ FINE_GROUND) * FINE_GROUND
        assert_allclose(np.vdot(FINE_GROUND, centered), 0, atol=1e-12)
        assert_allclose(P_K @ centered, centered, atol=1e-12)
    # with the site algebra the ground is cyclic and the defect lies in the generated space
    assert_allclose(generated_space(diagonal_units(5), FINE_GROUND), np.eye(5), atol=1e-12)

    # the same statements on the input side: the two-site path with its site algebra
    x = np.ones(2) / np.sqrt(2)
    sites = diagonal_units(2)
    assert_allclose(generated_space(sites, x), np.eye(2), atol=1e-12)
    centered = np.column_stack([X @ x - np.vdot(x, X @ x) * x for X in sites])
    assert_allclose(span_projection(centered), np.eye(2) - projection(x), atol=1e-12)


def lowest_projection(H, tol=1e-9):
    values, vectors = np.linalg.eigh(H)
    basis = vectors[:, values < values[0] + tol]
    return values[0], basis @ dag(basis), values


def test_ground_projection_and_three_independent_conditions():
    gap5 = 1 - np.cos(np.pi / 5)
    # unique and gapped, not cyclic for the cell algebra
    E0, P0, values = lowest_projection(path_hamiltonian(5))
    assert_allclose(E0, 0, atol=1e-12)
    assert rank(P0) == 1
    assert_allclose(P0, projection(FINE_GROUND), atol=1e-12)
    assert_allclose(values[1], gap5, atol=1e-12)
    assert rank(generated_space(cell_units(CELLS, 5), FINE_GROUND)) == 2

    # cyclic and gapped, not unique: two separate five-site paths
    H = np.zeros((10, 10))
    H[:5, :5] = path_hamiltonian(5)
    H[5:, 5:] = path_hamiltonian(5)
    E0, P0, values = lowest_projection(H)
    assert rank(P0) == 2
    assert_allclose(values[2], gap5, atol=1e-12)
    sites = diagonal_units(10)
    for weights in ((1, 1), (1, 2), (3, -1j)):
        y = unit(np.concatenate((weights[0] * np.ones(5), weights[1] * np.ones(5))))
        assert_allclose(H @ y, E0 * y, atol=1e-12)
        assert_allclose(generated_space(sites, y), np.eye(10), atol=1e-12)
    # the gap inequality H >= E0 + Delta (I - P0) holds although the ground is not unique
    assert np.linalg.eigvalsh(H - E0 * np.eye(10) - gap5 * (np.eye(10) - P0)).min() > -1e-12

    # unique and cyclic, not gapped: H e_0 = 0, H e_n = e_n / n, and the shift
    for N in (50, 200):
        H = np.diag(np.concatenate(([0.0], 1.0 / np.arange(1, N + 1))))
        S = np.diag(np.ones(N), -1)                 # S e_n = e_{n+1}, truncated
        e0 = np.eye(N + 1)[0]
        E0, P0, values = lowest_projection(H)
        assert rank(P0) == 1
        assert_allclose(P0, projection(e0), atol=1e-12)
        # the bottom of the spectrum on the complement of e_0 is 1/N and tends to zero
        assert_allclose(values[1], 1.0 / N, atol=1e-12)
        orbit = np.column_stack([np.linalg.matrix_power(S, k) @ e0 for k in range(N + 1)])
        assert_allclose(orbit, np.eye(N + 1), atol=1e-12)


def test_purity_and_faithfulness_on_the_five_site_path():
    y = FINE_GROUND
    omega = lambda X: np.vdot(y, X @ y)
    # site algebra: cyclic and separating, so faithful; not pure
    sites = diagonal_units(5)
    assert_allclose(generated_space(sites, y), np.eye(5), atol=1e-12)
    assert is_separating(sites, y)
    assert commutant_dimension(sites) == 5
    X = np.diag(random_complex(5))
    assert_allclose(omega(X), np.mean(np.diag(X)), atol=1e-12)
    # full matrix algebra: pure, not faithful
    full = matrix_units(5)
    assert commutant_dimension(full) == 1
    assert not is_separating(full, y)
    v = unit(random_complex(5))
    v -= np.vdot(y, v) * y
    killer = np.outer(np.eye(5)[0], v.conj())
    assert np.linalg.norm(killer) > 0.1
    assert_allclose(killer @ y, 0, atol=1e-12)
    # cell algebra: not cyclic, faithful on the cell algebra, a mixture with weights 1/5 and 4/5
    cells = cell_units(CELLS, 5)
    assert rank(generated_space(cells, y)) == 2
    assert is_separating(cells, y)
    c1, c2 = random_complex(2)
    assert_allclose(omega(c1 * cells[0] + c2 * cells[1]), c1 / 5 + 4 * c2 / 5, atol=1e-12)
    assert commutant_dimension(cells) == 1 + 16


def test_pure_and_faithful_only_for_the_scalars():
    # M_k (x) I_m on C^k (x) C^m with a vector of full Schmidt rank min(k, m)
    for k in (1, 2, 3):
        for m in (1, 2, 3):
            algebra = [np.kron(E, np.eye(m)) for E in matrix_units(k)]
            y = unit(random_complex(k * m))
            cyclic = rank(generated_space(algebra, y)) == k * m
            separating = is_separating(algebra, y)
            pure = commutant_dimension(algebra) == 1
            assert cyclic == (m <= k)
            assert separating == (k <= m)
            assert pure == (m == 1)
            if cyclic and pure and separating:
                assert k == m == 1
            # the GNS space of a faithful state has the dimension of the algebra
            coordinates, _ = gns(algebra, y)
            assert coordinates.shape[0] == (k * k if separating else rank(generated_space(algebra, y)))
    # In an irreducible algebra of matrices, the projection onto a vector orthogonal to y
    # is a nonzero element that annihilates y: transitivity in finite dimensions.
    y = unit(random_complex(4))
    eta = unit(random_complex(4))
    eta = unit(eta - np.vdot(y, eta) * y)
    X = projection(eta)
    assert_allclose(X @ y, 0, atol=1e-12)
    assert_allclose(X @ eta, eta, atol=1e-12)
