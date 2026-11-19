"""Checks for the clauses added to Sections 2 to 4: transport of observables by the
polar map, block forms of the graph maps, and the plane of the two reference vectors."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import block_diag, expm
from scipy.special import zeta

from lib.fock import dag
from lib.linalg import hermitian_power, rank, spectral_blocks

LABELS = {
    "LG": ["prop:pure-projector-correspondence", "thm:alphaA-star-isomorphism",
           "cor:paired-spectral-measures", "thm:WA-isometry", "prop:phase-graph-projection",
           "def:paired-observable-embedding", "thm:paired-algebra-representation",
           "thm:paired-algebra-graph-intertwining", "prop:mixed-state-correspondence",
           "eq:mixed-sector-compressions", "rem:direct-sum-not-partial-trace",
           "prop:phase-graph-factorization", "lem:ground-channel-bounded",
           "thm:ground-channel-polar", "prop:u0-uniqueness", "def:ground-graph-lines",
           "thm:ground-graphs-neutral", "prop:normalized-ground-pair",
           "prop:nonduality-without-identity", "def:ground-polar-operators",
           "thm:ground-polar-algebra", "prop:ground-graph-eigenspaces",
           "eq:Xi-Krein-cross-pairing", "prop:projectivized-neutral-locus",
           "id:companion-branch-status", "thm:ground-compression-full-polar",
           "prop:ground-bridge-equivalences", "prop:three-part-phenomenal-decomposition"],
    "UM": ["prop:matched-born-probabilities", "eq:polar-active-supports",
           "eq:polar-graph-projector", "sec:neutral-graphs-phase",
           "eq:restricted-ground-orientation", "prop:ground-polar-compatibility",
           "eq:actual-ground-partner"],
}

RNG = np.random.default_rng(20261119)
N_IN, N_OUT = 3, 5


def random_complex(*shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def unit(vector):
    return vector / np.linalg.norm(vector)


def projection(vector):
    return np.outer(vector, vector.conj())


def polar_parts(A):
    """The partial isometry that vanishes on ker A, and the modulus |A|."""
    left, singular, right = np.linalg.svd(A)
    # the examples have singular values either above 0.1 or zero to rounding
    assert np.sum(singular < 1e-10) == np.sum(singular < 0.1)
    r = int(np.sum(singular > 0.1))
    u = left[:, :r] @ right[:r]
    modulus = dag(right[:r]) @ np.diag(singular[:r]) @ right[:r]
    assert_allclose(u @ modulus, A, atol=1e-12)
    return u, modulus


def off_diagonal(lower):
    """Block operator on the direct sum with `lower` below the diagonal and its adjoint above."""
    n_out, n_in = lower.shape
    return np.block([[np.zeros((n_in, n_in)), dag(lower)], [lower, np.zeros((n_out, n_out))]])


def graph_map(V, phase):
    return np.vstack((np.eye(V.shape[1]), np.exp(1j * phase) * V)) / np.sqrt(2)


# An isometry from the input space onto a three-dimensional subspace D of the output space.
V = np.linalg.qr(random_complex(N_OUT, N_IN))[0]
PI_D = V @ dag(V)
ACTIVE = block_diag(np.eye(N_IN), PI_D)
GRADING = block_diag(np.eye(N_IN), -np.eye(N_OUT))
REFLECTION = off_diagonal(V)
COMPLEX_STRUCTURE = GRADING @ REFLECTION
PHASES = (0.0, 0.72, np.pi / 2, np.pi, 4.1)


def paired(X):
    return block_diag(X, V @ X @ dag(V))


def test_transport_of_projections_spectra_and_functions():
    assert_allclose(dag(V) @ V, np.eye(N_IN), atol=1e-12)
    psi, chi = unit(random_complex(N_IN)), unit(random_complex(N_IN))
    # rank-one projections and their trace overlaps
    assert_allclose(V @ projection(psi) @ dag(V), projection(V @ psi), atol=1e-12)
    overlap = abs(np.vdot(chi, psi)) ** 2
    assert_allclose(np.trace(projection(psi) @ projection(chi)).real, overlap, atol=1e-12)
    assert_allclose(np.trace(projection(V @ psi) @ projection(V @ chi)).real, overlap, atol=1e-12)

    # a self-adjoint observable with a repeated eigenvalue
    frame = np.linalg.qr(random_complex(N_IN, N_IN))[0]
    levels = np.array([0.3, 0.3, 1.2])
    X = frame @ np.diag(levels) @ dag(frame)
    X_D = V @ X @ dag(V)
    # on D the transported operator has the spectrum of X; outside D it vanishes
    assert_allclose(np.sort(np.linalg.eigvalsh(dag(V) @ X_D @ V)), np.sort(levels), atol=1e-12)
    assert_allclose(X_D @ (np.eye(N_OUT) - PI_D), 0, atol=1e-12)
    # resolvent, a polynomial, the exponential and the square root
    z = 0.4 + 0.9j
    resolvent = V @ np.linalg.inv(X - z * np.eye(N_IN)) @ dag(V)
    assert_allclose((X_D - z * PI_D) @ resolvent, PI_D, atol=1e-12)
    assert_allclose(X_D @ X_D @ X_D - 2 * X_D, V @ (X @ X @ X - 2 * X) @ dag(V), atol=1e-12)
    assert_allclose(expm(1j * X_D) @ PI_D, V @ expm(1j * X) @ dag(V), atol=1e-12)
    assert_allclose(hermitian_power(dag(V) @ X_D @ V, 0.5), hermitian_power(X, 0.5), atol=1e-12)

    # the spectral projections are transported, and the three distributions agree
    blocks = spectral_blocks(X)
    assert [round(float(value), 12) for value, _ in blocks] == [0.3, 1.2]
    total, weighted = np.zeros((N_OUT, N_OUT), complex), np.zeros((N_OUT, N_OUT), complex)
    graph = graph_map(V, 0.0)
    for value, basis in blocks:
        E = basis @ dag(basis)
        E_D = V @ E @ dag(V)
        assert_allclose(E_D @ E_D, E_D, atol=1e-12)
        assert_allclose(X_D @ E_D, value * E_D, atol=1e-12)
        total += E_D
        weighted += value * E_D
        p = np.vdot(psi, E @ psi).real
        assert_allclose(np.vdot(V @ psi, E_D @ V @ psi).real, p, atol=1e-12)
        assert_allclose(np.vdot(graph @ psi, paired(E) @ graph @ psi).real, p, atol=1e-12)
    assert_allclose(total, PI_D, atol=1e-12)
    assert_allclose(weighted, X_D, atol=1e-12)

    # extended by zero to the full output space, the outcome zero collects the defect and a
    # zero eigenvalue of the input observable
    X0 = frame @ np.diag([0.0, 0.7, 0.7]) @ dag(frame)
    extended = V @ X0 @ dag(V)
    values, vectors = np.linalg.eigh(extended)
    kernel = vectors[:, np.abs(values) < 1e-9]
    assert kernel.shape[1] == 1 + (N_OUT - N_IN)
    E0 = projection(frame[:, 0])
    assert_allclose(kernel @ dag(kernel), V @ E0 @ dag(V) + (np.eye(N_OUT) - PI_D), atol=1e-12)
    assert_allclose(np.sort(values), [0, 0, 0, 0.7, 0.7], atol=1e-12)

    # an operator that is not normal keeps its eigenvalues, and an effect stays an effect
    T = random_complex(N_IN, N_IN)
    assert_allclose(np.sort_complex(np.linalg.eigvals(dag(V) @ (V @ T @ dag(V)) @ V)),
                    np.sort_complex(np.linalg.eigvals(T)), atol=1e-10)
    effect = frame @ np.diag([0.0, 0.45, 1.0]) @ dag(frame)
    transported = np.linalg.eigvalsh(V @ effect @ dag(V))
    assert transported.min() > -1e-12 and transported.max() < 1 + 1e-12
    assert np.linalg.eigvalsh(PI_D - V @ effect @ dag(V)).min() > -1e-12


def test_graph_map_adjoint_and_the_two_reflection_eigenspaces():
    x, y = random_complex(N_IN), random_complex(N_OUT)
    for phase in PHASES:
        G = graph_map(V, phase)
        assert_allclose(dag(G) @ G, np.eye(N_IN), atol=1e-12)
        # adjoint on the active carrier and on the full direct sum
        for output in (PI_D @ y, y):
            assert_allclose(dag(G) @ np.concatenate((x, output)),
                            (x + np.exp(-1j * phase) * dag(V) @ output) / np.sqrt(2), atol=1e-12)
        blocks = np.block([[np.eye(N_IN), np.exp(-1j * phase) * dag(V)],
                           [np.exp(1j * phase) * V, PI_D]]) / 2
        assert_allclose(G @ dag(G), blocks, atol=1e-12)
        # two graphs at different phases
        for other in PHASES:
            assert_allclose(dag(G) @ graph_map(V, other),
                            (1 + np.exp(1j * (other - phase))) / 2 * np.eye(N_IN), atol=1e-12)

    aligned, opposite = graph_map(V, 0.0), graph_map(V, np.pi)
    P_plus, P_minus = aligned @ dag(aligned), opposite @ dag(opposite)
    assert_allclose(P_plus, (ACTIVE + REFLECTION) / 2, atol=1e-12)
    assert_allclose(P_minus, (ACTIVE - REFLECTION) / 2, atol=1e-12)
    assert_allclose(P_plus + P_minus, ACTIVE, atol=1e-12)
    assert_allclose(P_plus @ P_minus, 0, atol=1e-12)
    assert_allclose(REFLECTION @ aligned, aligned, atol=1e-12)
    assert_allclose(REFLECTION @ opposite, -opposite, atol=1e-12)
    assert rank(P_plus) == rank(P_minus) == N_IN and rank(ACTIVE) == 2 * N_IN


def sector_coordinates():
    """The unitary from S + D to C^2 (x) S, written on the full direct sum."""
    return np.block([[np.eye(N_IN), np.zeros((N_IN, N_OUT))],
                     [np.zeros((N_IN, N_IN)), dag(V)]])


def test_paired_representation_and_graph_compression():
    X, Y = random_complex(N_IN, N_IN), random_complex(N_IN, N_IN)
    assert_allclose(paired(X @ Y), paired(X) @ paired(Y), atol=1e-12)
    assert_allclose(paired(dag(X)), dag(paired(X)), atol=1e-12)
    assert_allclose(paired(X + 2j * Y), paired(X) + 2j * paired(Y), atol=1e-12)
    assert_allclose(paired(np.eye(N_IN)), ACTIVE, atol=1e-12)
    # X is read from either diagonal block
    assert_allclose(paired(X)[:N_IN, :N_IN], X, atol=1e-12)
    assert_allclose(dag(V) @ paired(X)[N_IN:, N_IN:] @ V, X, atol=1e-12)
    assert_allclose(paired(X) @ GRADING, GRADING @ paired(X), atol=1e-12)

    U = sector_coordinates()
    assert_allclose(U @ dag(U), np.eye(2 * N_IN), atol=1e-12)
    assert_allclose(dag(U) @ U, ACTIVE, atol=1e-12)
    assert_allclose(U @ paired(X) @ dag(U), np.kron(np.eye(2), X), atol=1e-12)
    for phase in PHASES:
        G = graph_map(V, phase)
        assert_allclose(paired(X) @ G, G @ X, atol=1e-12)
        assert_allclose(dag(G) @ paired(X) @ G, X, atol=1e-12)
        # the range of the graph map and its projection in the sector coordinates
        sector = np.array([1, np.exp(1j * phase)]) / np.sqrt(2)
        assert_allclose(U @ G, np.kron(sector.reshape(2, 1), np.eye(N_IN)), atol=1e-12)
        assert_allclose(U @ G @ dag(G) @ dag(U), np.kron(projection(sector), np.eye(N_IN)), atol=1e-12)

    # A block-diagonal operator that is not paired: the compression is the average of the
    # two blocks, and the graph is not invariant.
    G = graph_map(V, 0.72)
    B, L = random_complex(N_IN, N_IN), random_complex(N_IN, N_IN)
    T = block_diag(B, V @ L @ dag(V))
    assert_allclose(dag(G) @ T @ G, (B + L) / 2, atol=1e-12)
    leak = (np.eye(N_IN + N_OUT) - G @ dag(G)) @ T @ G
    assert_allclose(leak, graph_map(V, 0.72 + np.pi) @ (B - L) / 2, atol=1e-12)
    assert np.linalg.norm(leak) > 0.1


def test_graph_density_blocks_compressions_and_partial_trace():
    raw = random_complex(N_IN, N_IN)
    rho = raw @ dag(raw)
    rho /= np.trace(rho).real
    rho_D = V @ rho @ dag(V)
    X = random_complex(N_IN, N_IN)
    U = sector_coordinates()
    P_in = block_diag(np.eye(N_IN), np.zeros((N_OUT, N_OUT)))
    P_out = block_diag(np.zeros((N_IN, N_IN)), np.eye(N_OUT))
    for phase in PHASES:
        G = graph_map(V, phase)
        rho_graph = G @ rho @ dag(G)
        expected = np.block([[rho, np.exp(-1j * phase) * rho @ dag(V)],
                             [np.exp(1j * phase) * V @ rho, rho_D]]) / 2
        assert_allclose(rho_graph, expected, atol=1e-12)
        assert_allclose(np.trace(rho_graph).real, 1, atol=1e-12)
        assert np.linalg.eigvalsh(rho_graph).min() > -1e-12
        # compressions to the two sectors have trace one half
        assert_allclose(P_in @ rho_graph @ P_in, block_diag(rho / 2, np.zeros((N_OUT, N_OUT))), atol=1e-12)
        assert_allclose(P_out @ rho_graph @ P_out, block_diag(np.zeros((N_IN, N_IN)), rho_D / 2), atol=1e-12)
        # in the sector coordinates the state is a product; the partial trace returns rho
        product = U @ rho_graph @ dag(U)
        sector = np.array([1, np.exp(1j * phase)]) / np.sqrt(2)
        assert_allclose(product, np.kron(projection(sector), rho), atol=1e-12)
        reduced = product[:N_IN, :N_IN] + product[N_IN:, N_IN:]
        assert_allclose(reduced, rho, atol=1e-12)
        for operator, state in ((X, rho), (V @ X @ dag(V), rho_D), (paired(X), rho_graph)):
            assert_allclose(np.trace(state @ operator), np.trace(rho @ X), atol=1e-12)
        # removing the cross-sector blocks gives the equal mixture of two opposite graphs
        opposite = graph_map(V, phase + np.pi)
        dephased = P_in @ rho_graph @ P_in + P_out @ rho_graph @ P_out
        assert_allclose(dephased, (rho_graph + opposite @ rho @ dag(opposite)) / 2, atol=1e-12)
    # two phases differ only in the off-diagonal blocks
    difference = graph_map(V, 0.72) @ rho @ dag(graph_map(V, 0.72)) - graph_map(V, 4.1) @ rho @ dag(graph_map(V, 4.1))
    assert_allclose(P_in @ difference @ P_in, 0, atol=1e-12)
    assert_allclose(P_out @ difference @ P_out, 0, atol=1e-12)
    assert np.linalg.norm(difference) > 0.1


def weak_bridge_example():
    """A map with a kernel and reference vectors with A x = a y, x not an eigenvector of A*A."""
    left = np.linalg.qr(random_complex(N_OUT, N_OUT))[0]
    right = np.linalg.qr(random_complex(N_IN, N_IN))[0]
    A = left[:, :2] @ np.diag([0.5, 1.7]) @ dag(right[:, :2])
    x = unit(0.6 * right[:, 0] + 0.7 * right[:, 1] + 0.3 * right[:, 2])
    a = np.linalg.norm(A @ x)
    return A, x, A @ x / a, a


def test_restricted_reference_map():
    A, x, y, a = weak_bridge_example()
    P_x, P_y = projection(x), projection(y)
    u0 = np.outer(y, x.conj())
    restricted = A @ P_x
    assert_allclose(restricted, a * u0, atol=1e-12)
    assert rank(restricted) == 1
    assert_allclose(np.linalg.norm(restricted, 2), a, atol=1e-12)
    assert_allclose(dag(restricted) @ restricted, a ** 2 * P_x, atol=1e-12)
    assert_allclose(restricted @ dag(restricted), a ** 2 * P_y, atol=1e-12)
    # the positive square root of a positive operator with a kernel
    values, vectors = np.linalg.eigh(dag(restricted) @ restricted)
    assert values.min() > -1e-12
    modulus = vectors @ np.diag(np.sqrt(values.clip(min=0))) @ dag(vectors)
    assert_allclose(modulus, a * P_x, atol=1e-7)
    assert_allclose(u0 @ modulus, restricted, atol=1e-12)
    assert_allclose(u0 @ x, y, atol=1e-12)
    assert_allclose(dag(u0) @ y, x, atol=1e-12)
    assert_allclose(dag(u0) @ u0, P_x, atol=1e-12)
    assert_allclose(u0 @ dag(u0), P_y, atol=1e-12)
    # the bridge without representatives, and its converse
    assert_allclose(P_y @ restricted, restricted, atol=1e-12)
    assert_allclose((np.eye(N_OUT) - P_y) @ restricted, 0, atol=1e-12)
    assert_allclose(restricted @ dag(restricted) / np.linalg.norm(A @ x) ** 2, P_y, atol=1e-12)
    other = unit(random_complex(N_OUT))               # a line that does not contain A x
    assert np.linalg.norm((np.eye(N_OUT) - projection(other)) @ restricted) > 0.1
    # the adjoint of the restriction extends P_x A*
    assert_allclose(dag(restricted), P_x @ dag(A), atol=1e-12)

    # partial isometries between the two lines differ by a phase; the positive coefficient selects one
    for beta in (0.0, 0.4, np.pi / 2, np.pi, 5.0):
        w = np.exp(1j * beta) * u0
        assert_allclose(dag(w) @ w, P_x, atol=1e-12)
        assert_allclose(w @ dag(w), P_y, atol=1e-12)
        assert_allclose(np.linalg.norm(restricted - a * w, 2), a * abs(1 - np.exp(1j * beta)), atol=1e-12)
    # a common phase of the two representatives changes nothing; a phase on one of them breaks the bridge
    alpha = 1.3
    assert_allclose(np.outer(np.exp(1j * alpha) * y, (np.exp(1j * alpha) * x).conj()), u0, atol=1e-12)
    assert_allclose(A @ (np.exp(1j * alpha) * x), a * np.exp(1j * alpha) * y, atol=1e-12)
    assert np.linalg.norm(A @ (np.exp(1j * alpha) * x) - a * y) > 0.1

    # the full polar map does not restrict to u0 here
    u, modulus_A = polar_parts(A)
    assert np.linalg.norm(u @ P_x - u0) > 0.1
    assert np.linalg.norm(modulus_A @ x - a * x) > 0.1

    # A = diag(1, 2, 3, ...) and x with components proportional to n^-2: closed forms
    a_squared = zeta(2) / zeta(4)
    assert_allclose(a_squared, 15 / np.pi ** 2, rtol=1e-12)
    fidelity = zeta(3) ** 2 / (zeta(2) * zeta(4))
    assert_allclose(fidelity, 0.8116, atol=5e-5)
    n = np.arange(1, 200001, dtype=float)
    c = n ** -2.0 / np.sqrt(zeta(4))
    assert_allclose(np.sum(c ** 2), 1, atol=1e-12)
    assert_allclose(np.sum((n * c) ** 2), a_squared, atol=1e-5)                      # x is in Dom A
    assert_allclose(np.sum(n * c ** 2) ** 2 / np.sum((n * c) ** 2), fidelity, atol=1e-5)
    y_series = n * c / np.sqrt(a_squared)
    partial = np.cumsum((n * y_series) ** 2)                                         # ||A y||^2, truncated
    assert_allclose(partial[[999, 99999]], np.array([1000, 100000]) / (zeta(4) * a_squared), rtol=1e-12)


def plane_operators(x, y):
    u0 = np.outer(y, x.conj())
    R0 = off_diagonal(u0)
    return u0, R0, GRADING @ R0, block_diag(projection(x), projection(y))


def plane_vector(x, y, phase):
    return np.concatenate((x, np.exp(1j * phase) * y)) / np.sqrt(2)


def test_plane_operators_and_the_circle_of_neutral_rays():
    x, y = unit(random_complex(N_IN)), unit(random_complex(N_OUT))
    u0, R0, I0, P_xy = plane_operators(x, y)
    Q0 = -1j * I0
    J = GRADING
    assert_allclose(dag(R0), R0, atol=1e-12)
    assert_allclose(R0 @ R0, P_xy, atol=1e-12)
    assert_allclose(J @ R0, -R0 @ J, atol=1e-12)
    assert_allclose(dag(I0), -I0, atol=1e-12)
    assert_allclose(I0 @ I0, -P_xy, atol=1e-12)
    assert_allclose(dag(Q0), Q0, atol=1e-12)
    assert_allclose(Q0 @ Q0, P_xy, atol=1e-12)
    assert_allclose(J @ R0, 1j * Q0, atol=1e-12)
    assert_allclose(R0 @ Q0, 1j * J @ P_xy, atol=1e-12)
    assert_allclose(Q0 @ J, 1j * R0, atol=1e-12)
    for operator in (R0, I0, Q0):
        assert_allclose(operator @ P_xy, operator, atol=1e-12)
        assert_allclose(P_xy @ operator, operator, atol=1e-12)
    # Pauli matrices in the basis x + 0, 0 + y
    basis = np.column_stack((np.concatenate((x, np.zeros(N_OUT))), np.concatenate((np.zeros(N_IN), y))))
    pauli = {"z": np.diag([1, -1]), "x": np.array([[0, 1], [1, 0]]), "y": np.array([[0, -1j], [1j, 0]])}
    for operator, name in ((J, "z"), (R0, "x"), (Q0, "y")):
        assert_allclose(dag(basis) @ operator @ basis, pauli[name], atol=1e-12)

    z_plus, z_minus = plane_vector(x, y, 0.0), plane_vector(x, y, np.pi)
    assert_allclose(R0 @ z_plus, z_plus, atol=1e-12)
    assert_allclose(R0 @ z_minus, -z_minus, atol=1e-12)
    assert_allclose(J @ z_plus, z_minus, atol=1e-12)
    assert_allclose(J @ z_minus, z_plus, atol=1e-12)
    assert_allclose(I0 @ z_plus, z_minus, atol=1e-12)
    assert_allclose(I0 @ z_minus, -z_plus, atol=1e-12)
    # a Witt pair: orthonormal for the Hilbert product, hyperbolic for the Krein form
    krein = lambda v, w: np.vdot(v, J @ w)
    assert_allclose([np.vdot(z_plus, z_plus), np.vdot(z_minus, z_minus), np.vdot(z_plus, z_minus)], [1, 1, 0], atol=1e-12)
    assert_allclose([krein(z_plus, z_plus), krein(z_minus, z_minus), krein(z_plus, z_minus)], [0, 0, 1], atol=1e-12)

    P_in = block_diag(np.eye(N_IN), np.zeros((N_OUT, N_OUT)))
    P_out = np.eye(N_IN + N_OUT) - P_in
    for phase in PHASES:
        z = plane_vector(x, y, phase)
        assert_allclose([np.vdot(z, z), krein(z, z)], [1, 0], atol=1e-12)
        assert_allclose([np.vdot(z, R0 @ z), np.vdot(z, Q0 @ z), np.vdot(z, J @ z)],
                        [np.cos(phase), np.sin(phase), 0], atol=1e-12)
        # half of the norm squared in each sector, and the same two sector rays at every phase
        assert_allclose([np.linalg.norm(P_in @ z) ** 2, np.linalg.norm(P_out @ z) ** 2], [0.5, 0.5], atol=1e-12)
        assert_allclose(projection(unit((P_in @ z)[:N_IN])), projection(x), atol=1e-12)
        assert_allclose(projection(unit((P_out @ z)[N_IN:])), projection(y), atol=1e-12)
        # within the plane the Krein-orthogonal complement of the line is the line itself
        orthogonal = plane_vector(x, y, phase + np.pi)
        assert_allclose(np.vdot(z, orthogonal), 0, atol=1e-12)
        assert_allclose(krein(z, z), 0, atol=1e-12)
        assert abs(krein(z, orthogonal)) > 0.9
        # the projection onto the line
        assert_allclose(projection(z), (P_xy + np.cos(phase) * R0 + np.sin(phase) * Q0) / 2, atol=1e-12)

    # sector Hamiltonians with x and y as eigenvectors: equal mean energies on the circle,
    # and the evolution moves the phase at the rate of the energy difference
    E_x, E_y, hbar = -0.4, 1.1, 1.0
    H_in = E_x * projection(x) + 2.0 * (np.eye(N_IN) - projection(x))
    H_out = E_y * projection(y) + 3.0 * (np.eye(N_OUT) - projection(y))
    H = block_diag(H_in, H_out)
    for phase in PHASES:
        z = plane_vector(x, y, phase)
        assert_allclose(np.vdot(z, H @ z).real, (E_x + E_y) / 2, atol=1e-12)
        for t in (0.3, 2.0):
            evolved = expm(-1j * t * H / hbar) @ z
            target = np.exp(-1j * E_x * t / hbar) * plane_vector(x, y, phase + (E_x - E_y) * t / hbar)
            assert_allclose(evolved, target, atol=1e-12)
    # an operator commuting with J has the same mean at every phase
    block = block_diag(random_complex(N_IN, N_IN), random_complex(N_OUT, N_OUT))
    means = [np.vdot(plane_vector(x, y, phase), block @ plane_vector(x, y, phase)) for phase in PHASES]
    assert_allclose(means, means[0], atol=1e-12)

    # every unit vector of the plane with zero Krein square lies on the circle
    for _ in range(20):
        common, relative = RNG.uniform(0, 2 * np.pi, size=2)
        vector = np.exp(1j * common) * plane_vector(x, y, relative)
        alpha, beta = np.vdot(basis[:, 0], vector), np.vdot(basis[:, 1], vector)
        assert_allclose([abs(alpha), abs(beta)], [np.sqrt(0.5)] * 2, atol=1e-12)
        assert_allclose(projection(vector), projection(plane_vector(x, y, np.angle(beta / alpha))), atol=1e-12)
    for _ in range(20):
        alpha, beta = unit(random_complex(2))
        vector = basis @ np.array([alpha, beta])
        assert_allclose(krein(vector, vector), abs(alpha) ** 2 - abs(beta) ** 2, atol=1e-12)
    # the plane itself is not neutral, so a neutral subspace of it has dimension at most one
    assert_allclose(np.linalg.eigvalsh(dag(basis) @ J @ basis), [-1, 1], atol=1e-12)


def full_polar_operators(u):
    reflection = off_diagonal(u)
    return reflection, GRADING @ reflection, block_diag(dag(u) @ u, u @ dag(u))


def test_plane_and_the_full_polar_map():
    # strong compatibility: the plane reduces both polar operators
    left = np.linalg.qr(random_complex(N_OUT, N_OUT))[0]
    right = np.linalg.qr(random_complex(N_IN, N_IN))[0]
    A = left[:, :2] @ np.diag([0.5, 1.7]) @ dag(right[:, :2])
    x, y, a = right[:, 1], left[:, 1], 1.7
    assert_allclose(A @ x, a * y, atol=1e-12)
    u, _ = polar_parts(A)
    reflection, structure, active = full_polar_operators(u)
    u0, R0, I0, P_xy = plane_operators(x, y)
    assert_allclose(u @ x, y, atol=1e-12)
    for full, plane in ((reflection, R0), (structure, I0)):
        assert_allclose(full @ P_xy, plane, atol=1e-12)
        assert_allclose(P_xy @ full, plane, atol=1e-12)
    assert_allclose(P_xy @ active, P_xy, atol=1e-12)
    V_active = u @ right[:, :2]                       # the unitary from S to D in a basis of S
    for phase in PHASES:
        on_graph = np.concatenate((right[:, :2], np.exp(1j * phase) * V_active)) / np.sqrt(2)
        assert_allclose(on_graph @ dag(right[:, :2]) @ x, plane_vector(x, y, phase), atol=1e-12)

    # weak bridge only: the compressions are the plane operators times the square root of the fidelity
    for case in range(3):
        if case < 2:
            A, x, y, a = weak_bridge_example()
        else:                                          # the unequal-cell path of Section 4
            A = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]]) / 2
            x, y, a = np.ones(2) / np.sqrt(2), np.ones(5) / np.sqrt(5), np.sqrt(5 / 8)
        n_in = A.shape[1]
        grading = block_diag(np.eye(n_in), -np.eye(N_OUT))
        u, modulus = polar_parts(A)
        assert_allclose(A @ x, a * y, atol=1e-12)
        r = dag(u) @ y
        assert_allclose(r, modulus @ x / a, atol=1e-12)
        fidelity = abs(np.vdot(x, r)) ** 2
        assert np.vdot(x, r).real > 0 and abs(np.vdot(x, r).imag) < 1e-12
        assert fidelity < 0.95
        reflection = off_diagonal(u)
        structure = grading @ reflection
        active = block_diag(dag(u) @ u, u @ dag(u))
        u0 = np.outer(y, x.conj())
        R0 = off_diagonal(u0)
        I0 = grading @ R0
        P_xy = block_diag(projection(x), projection(y))
        assert_allclose(P_xy @ reflection @ P_xy, np.sqrt(fidelity) * R0, atol=1e-12)
        assert_allclose(P_xy @ structure @ P_xy, np.sqrt(fidelity) * I0, atol=1e-12)
        assert np.linalg.norm(P_xy @ reflection - reflection @ P_xy) > 0.1
        # the balanced vector and the two graphs of the full polar map
        z_plus = np.concatenate((x, y)) / np.sqrt(2)
        x_S = dag(u) @ u @ x
        in_support = np.linalg.norm(x_S) ** 2
        aligned = (active + reflection) / 2
        opposite = (active - reflection) / 2
        assert_allclose(np.vdot(z_plus, aligned @ z_plus).real, (in_support + 1) / 4 + np.sqrt(fidelity) / 2, atol=1e-12)
        assert_allclose(np.vdot(z_plus, opposite @ z_plus).real, (in_support + 1) / 4 - np.sqrt(fidelity) / 2, atol=1e-12)
        assert_allclose(np.vdot(z_plus, (np.eye(n_in + N_OUT) - active) @ z_plus).real, (1 - in_support) / 2, atol=1e-12)
        # the two components are the graph vectors of (x_S + r)/2 and (x_S - r)/2
        assert_allclose(aligned @ z_plus, np.concatenate(((x_S + r) / 2, u @ (x_S + r) / 2)) / np.sqrt(2), atol=1e-12)
        assert_allclose(opposite @ z_plus, np.concatenate(((x_S - r) / 2, -u @ (x_S - r) / 2)) / np.sqrt(2), atol=1e-12)
        # three parts of the output space: the ground line, the rest of D, and the defect
        P_y, P_D = projection(y), u @ dag(u)
        parts = (P_y, P_D - P_y, np.eye(N_OUT) - P_D)
        assert_allclose(P_D @ y, y, atol=1e-12)
        for index, part in enumerate(parts):
            assert_allclose(part @ part, part, atol=1e-12)
            for later in parts[index + 1:]:
                assert_allclose(part @ later, 0, atol=1e-12)
        assert [rank(part) for part in parts] == [1, rank(P_D) - 1, N_OUT - rank(P_D)]
        assert_allclose(dag(A) @ parts[2], 0, atol=1e-12)
        # V carries the complement of r in S onto the complement of y in D
        P_S, P_r = dag(u) @ u, projection(r)
        assert_allclose(u @ (P_S - P_r) @ dag(u), P_D - P_y, atol=1e-12)
        if case == 2:
            assert_allclose(fidelity, 0.9, atol=1e-12)
            assert_allclose(in_support, 1, atol=1e-12)
            assert_allclose(np.vdot(z_plus, aligned @ z_plus).real, (1 + np.sqrt(0.9)) / 2, atol=1e-12)
            assert_allclose((1 + np.sqrt(0.9)) / 2, 0.9743, atol=5e-5)
        else:
            assert in_support < 0.95
