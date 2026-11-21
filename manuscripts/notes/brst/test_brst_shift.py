"""Checks for the note "The shift example as a BRST model".

Every statement of the note that can be tested in finite dimensions is checked here.
The shift is modelled exactly: input C^N with basis f_0..f_{N-1}, output C^{N+1} with basis
g_0..g_N, A f_n = g_{n+1}. Its energies involve only the levels 0, 1 and eps, so the
truncation introduces no error.

Run:  python test_brst_shift.py      (or with pytest)
"""
import numpy as np
from scipy.linalg import expm, null_space, orth

RNG = np.random.default_rng(20261121)
TOL = 1e-10


def rc(*shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def herm(k):
    x = rc(k, k)
    return (x + x.conj().T) / 2


def block(*blocks):
    size = sum(b.shape[0] for b in blocks)
    out = np.zeros((size, size), complex)
    i = 0
    for b in blocks:
        k = b.shape[0]
        out[i:i + k, i:i + k] = b
        i += k
    return out


def proj(v):
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def polar_data(A):
    """Polar partial isometry V, projections onto S = (ker A)^perp and D = closure Ran A."""
    U, s, Vh = np.linalg.svd(A)
    r = int((s > TOL).sum())
    V = U[:, :r] @ Vh[:r]
    PS = Vh[:r].conj().T @ Vh[:r]
    PD = U[:, :r] @ U[:, :r].conj().T
    return V, PS, PD, r


def grading(m, n):
    return np.diag([1.0] * m + [-1.0] * n)


def charge(A):
    """Q_A = P_G J with G the aligned graph of the polar map."""
    n, m = A.shape
    V, PS, PD, _ = polar_data(A)
    PG = 0.5 * np.block([[PS, V.conj().T], [V, PD]])
    return PG @ grading(m, n), PG


def cohomology(Q, J):
    """Representatives of ker Q / Ran Q (orthogonal complement of Ran Q in ker Q) and the
    signature of the induced Krein form."""
    K = null_space(Q, rcond=TOL)
    R = orth(Q, rcond=TOL)
    H = K @ null_space(R.conj().T @ K, rcond=TOL) if R.size else K
    w = np.linalg.eigvalsh(H.conj().T @ J @ H) if H.size else np.array([])
    return H, int((w > TOL).sum()), int((w < -TOL).sum())


def in_span(basis, v):
    if basis.size == 0:
        return np.allclose(v, 0, atol=1e-9)
    P = basis @ basis.conj().T
    return np.allclose(P @ v, v, atol=1e-9)


# ---- the shift model -------------------------------------------------------------------

N = 6


def shift(N=N):
    A = np.zeros((N + 1, N))
    for k in range(N):
        A[k + 1, k] = 1
    return A


def hamiltonians(eps, N=N):
    H_psi = np.eye(N) - proj(np.eye(N)[0])                                 # I - P_{f0}
    g = np.eye(N + 1)
    H_phi = np.eye(N + 1) - proj(g[1]) + (eps - 1) * proj(g[0])          # I - P_{g1} + (eps-1) P_{g0}
    return H_psi, H_phi


def rotated_shift(alpha, N=N):
    U = np.eye(N + 1)
    c, s = np.cos(alpha), np.sin(alpha)
    U[0, 0], U[2, 0], U[0, 2], U[2, 2] = c, s, -s, c                     # U g0 = cos g0 + sin g2
    return U @ shift(N)


# ---- Proposition 1: the charge, its kernel and its cohomology ------------------------------

def test_charge_is_nilpotent_krein_selfadjoint_with_graph_range():
    for A in (rc(5, 3), rc(5, 2) @ rc(2, 3), shift(), rc(4, 4)):
        n, m = A.shape
        J = grading(m, n)
        Q, PG = charge(A)
        V, PS, PD, r = polar_data(A)
        assert np.allclose(Q @ Q, 0, atol=TOL)
        assert np.allclose(J @ Q.conj().T @ J, Q, atol=TOL)                 # Krein self-adjoint
        assert np.allclose(PG @ PG, PG, atol=TOL) and np.allclose(PG, PG.conj().T, atol=TOL)
        assert np.linalg.matrix_rank(Q, tol=TOL) == r                       # Ran Q = G, dim G = dim S
        # ker Q = G (+) (ker A (+) D^perp)
        kerA = np.eye(m) - PS
        defect = np.eye(n) - PD
        N0 = np.block([[kerA, np.zeros((m, n))], [np.zeros((n, m)), defect]])
        assert np.allclose(Q @ N0, 0, atol=TOL)
        assert np.linalg.matrix_rank(null_space(Q, rcond=TOL), tol=TOL) == r + (m - r) + (n - r)
        # Ran Q is Krein-orthogonal to ker Q, so the form descends to the cohomology
        K = null_space(Q, rcond=TOL)
        assert np.allclose((Q @ rc(m + n, 3)).conj().T @ J @ K, 0, atol=1e-9)
        H, plus, minus = cohomology(Q, J)
        assert (plus, minus) == (m - r, n - r)


def test_witt_pairs():
    A = rc(5, 3)
    n, m = A.shape
    J = grading(m, n)
    Q, _ = charge(A)
    V, PS, _, _ = polar_data(A)
    for _ in range(3):
        psi = PS @ rc(m)
        psi /= np.linalg.norm(psi)
        aligned = np.r_[psi, V @ psi] / np.sqrt(2)
        companion = np.r_[psi, -V @ psi] / np.sqrt(2)
        assert np.allclose(companion, J @ aligned)
        assert np.allclose(Q @ aligned, 0, atol=TOL)
        assert np.allclose(Q @ companion, aligned, atol=TOL)
        assert abs(np.vdot(aligned, J @ aligned)) < TOL
        assert abs(np.vdot(companion, J @ companion)) < TOL
        assert np.isclose(np.vdot(aligned, J @ companion), 1)


# ---- Lemma 2 and Corollary 3: even invariant operators ------------------------------------

def even_invariant(A, X_psi_S, X_ker, X_def):
    """Build X = X_psi (+) X_phi from data on S, on ker A and on D^perp (Lemma 2)."""
    n, m = A.shape
    V, PS, PD, _ = polar_data(A)
    kerA, defect = np.eye(m) - PS, np.eye(n) - PD
    X_psi = PS @ X_psi_S @ PS + kerA @ X_ker @ kerA
    X_phi = V @ X_psi @ V.conj().T + defect @ X_def @ defect
    return X_psi, X_phi


def test_even_operator_invariant_iff_lemma_conditions():
    for A in (rc(5, 3), rc(5, 2) @ rc(2, 3)):
        n, m = A.shape
        Q, _ = charge(A)
        V, PS, PD, _ = polar_data(A)
        X_psi, X_phi = even_invariant(A, rc(m, m), rc(m, m), rc(n, n))     # need not be self-adjoint
        X = block(X_psi, X_phi)
        assert np.allclose(X @ Q, Q @ X, atol=1e-9)
        assert np.allclose(X_phi @ V, V @ X_psi, atol=1e-9)
        assert np.allclose(X_psi @ V.conj().T, V.conj().T @ X_phi, atol=1e-9)
        # breaking any one condition breaks invariance
        broken = [block(X_psi, X_phi + 0.1 * (np.eye(n) - PD) @ rc(n, n) @ PD),       # mixes D with defect
                  block(X_psi + 0.1 * PS @ rc(m, m) @ PS, X_phi)]                    # breaks intertwining
        if np.linalg.matrix_rank(np.eye(m) - PS, tol=TOL):                           # only if ker A != 0
            broken.append(block(X_psi + 0.1 * (np.eye(m) - PS) @ rc(m, m) @ PS, X_phi))
        for bad in broken:
            assert not np.allclose(bad @ Q, Q @ bad, atol=1e-6)
        # action on the cohomology: X_psi on ker A and X_phi on the defect
        J = grading(m, n)
        H, _, _ = cohomology(Q, J)
        XH = X @ H
        assert np.allclose(H @ np.linalg.lstsq(H, XH, rcond=None)[0], XH, atol=1e-9)


def test_paired_effects_are_invariant_and_vanish_on_the_physical_space():
    for A in (rc(5, 3), shift()):
        n, m = A.shape
        Q, _ = charge(A)
        V, PS, _, _ = polar_data(A)
        J = grading(m, n)
        H, _, _ = cohomology(Q, J)
        E0 = rc(m, m)
        E = PS @ (E0 @ E0.conj().T) @ PS
        E /= np.linalg.eigvalsh(E).max()
        X = block(E, V @ E @ V.conj().T)
        assert np.allclose(X @ Q, Q @ X, atol=1e-9)
        assert np.allclose(X @ H, 0, atol=1e-9)


def test_invariant_hamiltonian_and_physical_evolution():
    A = rc(5, 3)
    n, m = A.shape
    Q, _ = charge(A)
    V, PS, PD, _ = polar_data(A)
    defect = np.eye(n) - PD
    H_psi = herm(m)
    H_ex = herm(n)
    H_phi = V @ H_psi @ V.conj().T + defect @ H_ex @ defect
    H = block(H_psi, H_phi)
    assert np.allclose(H @ Q, Q @ H, atol=1e-9)
    J = grading(m, n)
    Hc, _, _ = cohomology(Q, J)
    for t in (0.3, 1.7):
        Ut = expm(-1j * t * H)
        assert np.allclose(Ut @ Q, Q @ Ut, atol=1e-9)
        # restricted to the physical space the evolution is that of the excess Hamiltonian
        lower = Hc[m:]
        assert np.allclose(Ut @ Hc, np.r_[np.zeros((m, Hc.shape[1])),
                                          expm(-1j * t * defect @ H_ex @ defect) @ lower], atol=1e-9)


# ---- Section 4: the shift ------------------------------------------------------------------

def test_shift_physical_space_is_the_missing_line():
    A = shift()
    n, m = A.shape
    J = grading(m, n)
    Q, _ = charge(A)
    H, plus, minus = cohomology(Q, J)
    assert (plus, minus) == (0, 1)
    g0 = np.r_[np.zeros(m), np.eye(n)[0]]
    assert in_span(H, g0)
    assert np.isclose(np.vdot(g0, J @ g0), -1)


def test_shift_dynamics_is_invariant_for_every_eps_and_the_line_evolves_alone():
    A = shift()
    n, m = A.shape
    Q, _ = charge(A)
    g0 = np.r_[np.zeros(m), np.eye(n)[0]]
    for eps in (0.4, 1.0, 2.0, 2.5):
        H_psi, H_phi = hamiltonians(eps)
        H = block(H_psi, H_phi)
        assert np.allclose(H_phi @ A, A @ H_psi)
        assert np.allclose(H @ Q, Q @ H, atol=TOL)
        assert np.allclose(expm(-1j * 0.8 * H) @ g0, np.exp(-1j * 0.8 * eps) * g0)


def test_rotated_shift_is_gauge_consistent_only_at_eps_one():
    for alpha in (0.3, 1.2):
        B = rotated_shift(alpha)
        n, m = B.shape
        QB, _ = charge(B)
        for eps, consistent in ((1.0, True), (2.5, False), (0.4, False)):
            H_psi, H_phi = hamiltonians(eps)
            H = block(H_psi, H_phi)
            assert np.isclose(np.linalg.norm(H_phi @ B - B @ H_psi, 2), abs(eps - 1) * abs(np.sin(alpha)))
            assert np.allclose(H @ QB, QB @ H, atol=TOL) == consistent


def test_spectral_obstruction_for_families():
    """Proposition 4: for eps not in {0, 1}, any family whose joint charge commutes with the
    dynamics keeps g0 in its joint defect. Families containing a non-intertwiner fail the test;
    families of intertwiners pass and keep g0."""
    A = shift()
    n, m = A.shape
    g0 = np.eye(n)[0]
    # an intertwiner: maps f0 to the output ground and the excited inputs into span{g2..gN}
    C = np.zeros((n, m), complex)
    C[1, 0] = np.exp(0.4j)
    C[2:, 1:] = np.linalg.qr(rc(n - 2, m - 1))[0][:, :m - 1]
    channel_map = np.zeros((n, m))
    channel_map[0, 1] = 1                                                   # f1 -> g0 (covariant channel)
    for eps in (2.5, 0.4, 1.0):
        H_psi, H_phi = hamiltonians(eps)
        for members, label in (((A, C), "intertwiners"),
                               ((A, rotated_shift(0.7)), "with rotated shift"),
                               ((A, channel_map), "with channel map")):
            B = np.hstack(members)
            QB, _ = charge(B)
            Hj = block(*([H_psi] * len(members)), H_phi)
            invariant = np.allclose(Hj @ QB, QB @ Hj, atol=TOL)
            _, _, PD = polar_data(B)[:3]
            g0_in_defect = np.allclose(PD @ g0, 0, atol=1e-9)
            if label == "intertwiners":
                assert invariant and g0_in_defect
            elif eps != 1.0:
                assert not invariant                                       # obstruction
            if invariant and eps not in (0.0, 1.0):
                assert g0_in_defect                                        # Proposition 4
    # at eps = 1 the rotated shift is admissible and the joint family has no defect at all
    H_psi, H_phi = hamiltonians(1.0)
    B = np.hstack((A, rotated_shift(0.7)))
    QB, _ = charge(B)
    Hj = block(H_psi, H_psi, H_phi)
    assert np.allclose(Hj @ QB, QB @ Hj, atol=TOL)
    assert np.linalg.matrix_rank(B, tol=TOL) == n


def test_redundant_presentations_make_the_space_indefinite_until_quotiented():
    A = shift()
    n, m = A.shape
    C = np.zeros((n, m), complex)
    C[1, 0] = 1j
    C[2:, 1:] = np.linalg.qr(rc(n - 2, m - 1))[0][:, :m - 1]
    B = np.hstack((A, C))
    QB, _ = charge(B)
    JB = grading(2 * m, n)
    _, plus, minus = cohomology(QB, JB)
    kernel_dim = 2 * m - np.linalg.matrix_rank(B, tol=TOL)
    assert kernel_dim > 0 and (plus, minus) == (kernel_dim, 1)
    # treat redundant presentations as gauge: restrict B to (ker B)^perp
    _, s, Vh = np.linalg.svd(B)
    r = int((s > TOL).sum())
    Br = B @ Vh[:r].conj().T
    Qr, _ = charge(Br)
    Hr, plus, minus = cohomology(Qr, grading(r, n))
    assert (plus, minus) == (0, 1)
    assert in_span(Hr, np.r_[np.zeros(r), np.eye(n)[0]])


def test_no_invariant_operator_moves_gauge_into_the_physical_space():
    A = shift()
    n, m = A.shape
    J = grading(m, n)
    Q, _ = charge(A)
    R = orth(Q, rcond=TOL)
    Hc, _, _ = cohomology(Q, J)
    # invariant operators keep Ran Q inside Ran Q
    X_psi, X_phi = even_invariant(A, rc(m, m), rc(m, m), rc(n, n))
    X = block(X_psi, X_phi)
    assert np.allclose(X @ Q, Q @ X, atol=1e-9)
    for v in (R @ rc(R.shape[1]) for _ in range(3)):
        assert in_span(R, X @ v)
    # the energy-transferring operator |g0><g1| is covariant but not invariant
    g = np.eye(n)
    K = block(np.zeros((m, m)), np.outer(g[0], g[1]))
    eps = 2.5
    H_psi, H_phi = hamiltonians(eps)
    assert np.allclose(H_phi @ K[m:, m:] - K[m:, m:] @ H_phi, eps * K[m:, m:])   # [H_phi, K] = eps K
    assert not np.allclose(K @ Q, Q @ K, atol=1e-6)
    assert np.allclose(block(H_psi, H_phi) @ Q, Q @ block(H_psi, H_phi), atol=TOL)   # the evolution premise holds
    z_plus = np.r_[np.eye(m)[0], g[1]] / np.sqrt(2)                         # ground event, exact
    assert in_span(R, z_plus)
    image = K @ z_plus                                                     # = (0 + g0)/sqrt2
    assert in_span(Hc, image) and not np.allclose(image, 0)                # a nonzero physical class
    # the second premise is independent of the Kraus representation: an isometric recombination
    # of operators that preserve Ran Q and ker Q preserves them too
    ops = [block(*even_invariant(A, rc(m, m), rc(m, m), rc(n, n))) for _ in range(2)]
    iso = np.linalg.qr(rc(3, 2))[0]
    kerQ = null_space(Q, rcond=TOL)
    for j in range(3):
        Kj = iso[j, 0] * ops[0] + iso[j, 1] * ops[1]
        assert all(in_span(R, Kj @ R[:, i]) for i in range(R.shape[1]))
        assert all(in_span(kerQ, Kj @ kerQ[:, i]) for i in range(kerQ.shape[1]))


def test_the_ground_is_pure_gauge():
    A = shift()
    n, m = A.shape
    J = grading(m, n)
    Q, PG = charge(A)
    R = orth(Q, rcond=TOL)
    f0 = np.r_[np.eye(m)[0], np.zeros(n)]
    g1 = np.r_[np.zeros(m), np.eye(n)[1]]
    z_plus = (f0 + g1) / np.sqrt(2)
    assert in_span(R, z_plus)                                              # exact
    assert np.allclose(Q @ g1, -z_plus / np.sqrt(2))                       # luminous ground alone: not closed
    assert np.allclose(Q @ f0, z_plus / np.sqrt(2))                        # vacuum alone: not closed
    Hc, _, _ = cohomology(Q, J)
    assert np.allclose(Hc.conj().T @ f0, 0) and np.allclose(Hc.conj().T @ g1, 0)
    # keeping the ground pair out of the gauge makes the physical space indefinite
    G = orth(PG, rcond=TOL)
    Gp = orth(G - np.outer(z_plus, z_plus.conj()) @ G, rcond=TOL)
    Qp = Gp @ Gp.conj().T @ J
    assert np.allclose(Qp @ Qp, 0, atol=TOL) and np.allclose(J @ Qp.conj().T @ J, Qp, atol=TOL)
    _, plus, minus = cohomology(Qp, J)
    assert (plus, minus) == (1, 2)


def test_evolution_acts_on_the_quotient_iff_it_commutes_with_the_charge():
    A = rc(5, 3)
    n, m = A.shape
    Q, PG = charge(A)
    V, PS, PD, _ = polar_data(A)
    defect = np.eye(n) - PD
    H_psi = herm(m)
    good = block(H_psi, V @ H_psi @ V.conj().T + defect @ herm(n) @ defect)
    bad = block(H_psi, V @ (H_psi + 0.2 * herm(m)) @ V.conj().T + defect @ herm(n) @ defect)
    R = orth(Q, rcond=TOL)
    for H, invariant in ((good, True), (bad, False)):
        preserved = all(np.allclose((np.eye(m + n) - R @ R.conj().T) @ expm(-1j * t * H) @ R, 0, atol=1e-9)
                        for t in (0.4, -1.3, 2.2))
        commutes = all(np.allclose(expm(-1j * t * H) @ Q, Q @ expm(-1j * t * H), atol=1e-9)
                       for t in (0.4, -1.3, 2.2))
        assert preserved == commutes == invariant


def test_it_is_the_polar_map_that_must_intertwine():
    H = np.diag([0.0, 1.0])
    A = np.array([[2.0, 0.5], [0.5, 1.0]])                                 # positive invertible, mixes energies
    assert not np.allclose(H @ A, A @ H)
    V, _, _, _ = polar_data(A)
    assert np.allclose(V, np.eye(2))
    Q, _ = charge(A)
    assert np.allclose(block(H, H) @ Q, Q @ block(H, H))


def test_spectral_obstruction_with_a_degenerate_eigenvalue():
    m, n, lam = 3, 6, 2.5
    H_in = np.diag([0.0, 1.0, 1.0])
    W = np.linalg.qr(rc(n, m))[0]                                          # isometry with defect of dimension 3
    defect = np.eye(n) - W @ W.conj().T
    Dperp = orth(defect, rcond=TOL)
    mix = np.linalg.qr(rc(3, 3))[0]
    H_out = W @ H_in @ W.conj().T + Dperp @ mix @ np.diag([lam, lam, 4.0]) @ mix.conj().T @ Dperp.conj().T
    Q, _ = charge(W)
    assert np.allclose(block(H_in, H_out) @ Q, Q @ block(H_in, H_out), atol=1e-9)
    vals, vecs = np.linalg.eigh(H_out)
    eig = vecs[:, np.isclose(vals, lam)]
    assert eig.shape[1] == 2                                               # lam is degenerate
    assert np.allclose(W.conj().T @ eig, 0, atol=1e-9)                      # every eigenvector in the defect


def test_resonant_family_keeps_redundancies_until_quotiented():
    A = shift()
    n, m = A.shape
    B = np.hstack((A, rotated_shift(0.7)))
    assert np.allclose(B @ np.r_[np.eye(m)[0], -np.eye(m)[0]], 0)          # (f0, -f0) is redundant
    QB, _ = charge(B)
    assert cohomology(QB, grading(2 * m, n))[1:] == (2 * m - n, 0)         # nonzero without any excess
    _, s, Vh = np.linalg.svd(B)
    r = int((s > TOL).sum())
    Qr, _ = charge(B @ Vh[:r].conj().T)
    assert cohomology(Qr, grading(r, n))[1:] == (0, 0)
    # the quarter-turn prepares the missing line exactly
    assert np.allclose(rotated_shift(np.pi / 2) @ np.eye(m)[1], -np.eye(n)[0])
    # direct return of g0: zero through the shift, one through the quarter-turn
    g0 = np.eye(n)[0]
    for M, ret in ((A, 0.0), (rotated_shift(np.pi / 2), 1.0)):
        V, _, _, _ = polar_data(M)
        assert np.isclose(np.linalg.norm(V.conj().T @ g0) ** 2, ret)


def test_relational_ground_needs_strong_compatibility():
    """On the unequal-cell path (weak but not strong compatibility) the manuscript's relational
    ground (x + y)/sqrt2 is not closed, and the ground plane is not invariant."""
    A = np.repeat(np.eye(2), [1, 4], axis=0) / 2
    x, y = np.ones(2) / np.sqrt(2), np.ones(5) / np.sqrt(5)
    Q, _ = charge(A)
    V, _, _, _ = polar_data(A)
    assert np.isclose((x @ V.conj().T @ y) ** 2, 0.9)
    z0 = np.r_[x, y] / np.sqrt(2)
    assert np.isclose(np.linalg.norm(Q @ z0) ** 2, (1 - np.sqrt(0.9)) / 2)
    assert abs(np.vdot(z0, grading(2, 5) @ z0)) < TOL                      # Krein-null
    # exact for the rank-one charge of the restricted ground map u0 = |y><x|
    Q0, _ = charge(np.outer(y, x))
    assert in_span(orth(Q0, rcond=TOL), z0)
    plane = np.column_stack([np.r_[x, np.zeros(5)], np.r_[np.zeros(2), y]])
    P = plane @ plane.T
    assert not np.allclose(P @ Q @ P, Q @ P)
    # the full-graph event over u*y is exact, and the physical space is orthogonal to y
    event = np.r_[V.conj().T @ y, y] / np.sqrt(2)
    assert in_span(orth(Q, rcond=TOL), event)
    Hc, _, _ = cohomology(Q, grading(2, 5))
    assert np.allclose(Hc[2:].conj().T @ y, 0)


def test_member_and_joint_conservation_are_independent():
    """Remark 'Members and families'. (a) Two injective members with the same ground, range and
    polar map I each have a conserved charge, while the joint charge is not conserved; with an
    extra output level lam, Proposition 6 applied to each member still puts its eigenvector in
    the common defect. (b) A joint charge can be conserved while neither member's charge is."""
    H = np.diag([0.0, 1.0, 2.0])
    A1 = np.eye(3)
    A2 = np.eye(3)
    A2[1:, 1:] = [[2, 1], [1, 2]]
    def conserved(M, H_in, H_out, copies):
        Q, _ = charge(M)
        T = block(*([H_in] * copies), H_out)
        return np.allclose(T @ Q, Q @ T, atol=TOL)
    assert conserved(A1, H, H, 1) and conserved(A2, H, H, 1)
    assert np.allclose(A2 @ np.eye(3)[0], A1 @ np.eye(3)[0])                # same ground
    assert not conserved(np.hstack([A1, A2]), H, H, 2)
    lam = 3.5
    H_out = np.diag([0.0, 1.0, 2.0, lam])
    W = np.eye(4)[:, :3]
    members = (W @ A1, W @ A2)
    assert all(conserved(M, H, H_out, 1) for M in members)
    B = np.hstack(members)
    assert not conserved(B, H, H_out, 2)
    assert np.allclose(B.conj().T @ np.eye(4)[3], 0)                        # in the common defect
    # (b) with a shared lower ground, so that R3 holds for both members
    H3 = np.diag([0.0, 1.0, 2.0])
    ug = np.zeros(6); ug[[0, 3]] = 1 / np.sqrt(2)
    u0 = np.zeros(6); u0[1], u0[4] = 0.9, np.sqrt(0.19)
    u1 = np.zeros(6); u1[2], u1[5] = 0.3, np.sqrt(0.91)
    U = np.column_stack([ug, u0, u1])
    absB = np.eye(3)
    absB[1:, 1:] = [[1.0, 0.5], [0.5, 1.0]]
    B = absB @ U.conj().T                                                   # W|B| with W u_g = h0, W u_j = h_{j+1}
    assert np.allclose(polar_data(B)[0], U.conj().T)
    assert conserved(B, H3, H3, 2)
    for c in range(2):
        Ac = B[:, 3 * c:3 * c + 3]
        assert np.allclose(Ac[:, 0], np.eye(3)[0] / np.sqrt(2))             # R3: e0 -> h0 / sqrt2
        assert not conserved(Ac, H3, H3, 1)
    # both examples use three input energies, unlike the shift


def test_shift_member_and_joint_conservation_coincide_under_R3():
    """Proposition 9: in the shift, for presentations with A f0 in C g1, conservation of every
    member's charge, conservation of the joint charge and intertwining by every member are
    equivalent, at resonance and away from it. Without that condition a conserved charge does not
    force intertwining."""
    A = shift()
    n, m = A.shape

    def conserved(M, H_in, H_out, copies):
        Q, _ = charge(M)
        T = block(*([H_in] * copies), H_out)
        return np.allclose(T @ Q, Q @ T, atol=1e-9)

    def intertwines(M, H_in, H_out):
        return np.allclose(H_out @ M, M @ H_in, atol=1e-9)

    def intertwiner(eps):
        W = np.zeros((n, m), complex)
        W[1, 0] = np.exp(1j * RNG.normal())
        targets = [0] + list(range(2, n)) if eps == 1.0 else list(range(2, n))
        W[targets, 1:] = np.linalg.qr(rc(len(targets), m - 1))[0][:, :m - 1]
        return W

    def excited_positive():
        Y = rc(m - 1, m - 1)
        P = np.eye(m, dtype=complex)
        P[1:, 1:] = Y @ Y.conj().T + 0.1 * np.eye(m - 1)
        return P

    for eps in (1.0, 2.5, 0.4):
        H_psi, H_phi = hamiltonians(eps)
        for _ in range(12):
            members = []
            for kind in RNG.integers(0, 3, size=3):
                if kind == 0:
                    members.append(intertwiner(eps) @ np.diag(RNG.uniform(0.5, 2, m)))
                elif kind == 1:
                    members.append(rotated_shift(RNG.uniform(0.2, 1.4)))
                else:
                    members.append(intertwiner(eps) @ excited_positive())
            for M in members:
                assert np.allclose(M[:, 0], M[1, 0] * np.eye(n)[1])        # A f0 in C g1
            B = np.hstack(members)
            a = all(conserved(M, H_psi, H_phi, 1) for M in members)
            b = conserved(B, H_psi, H_phi, len(members))
            c = all(intertwines(M, H_psi, H_phi) for M in members)
            assert a == b == c
    # the ground condition cannot be dropped
    H_psi, H_phi = hamiltonians(2.5)
    P = np.eye(m)
    P[:2, :2] = [[2.0, 0.5], [0.5, 1.0]]
    assert conserved(A @ P, H_psi, H_phi, 1) and not intertwines(A @ P, H_psi, H_phi)
    assert not np.allclose((A @ P)[:, 0], (A @ P)[1, 0] * np.eye(n)[1])


def test_scope_of_the_intervention_premise():
    """Section 4, Interventions. The premise is sufficient, not necessary; for sector-preserving
    Kraus operators it is commutation with the charge and protects the whole matched part; for
    sector-mixing ones it protects exact events only; a channel from inputs to appearances,
    extended by zero, fails it under R4."""
    A = shift()
    n, m = A.shape
    J = grading(m, n)
    Q, _ = charge(A)
    R = orth(Q, rcond=TOL)
    Kq = null_space(Q, rcond=TOL)

    def preserves(X, basis):
        return np.allclose((np.eye(X.shape[0]) - basis @ basis.conj().T) @ X @ basis, 0, atol=1e-9)

    f0 = np.r_[np.eye(m)[0], np.zeros(n)]
    g1 = np.r_[np.zeros(m), np.eye(n)[1]]
    g0 = np.r_[np.zeros(m), np.eye(n)[0]]
    zp, zm = (f0 + g1) / np.sqrt(2), (f0 - g1) / np.sqrt(2)
    # J: does not act on the quotient, yet preserves the missing-line weight of every state
    assert not preserves(J, R)
    v = rc(m + n)
    assert np.isclose(abs(np.vdot(g0, J @ v)), abs(np.vdot(g0, v)))
    # sector-preserving operators: acting on the quotient <=> commuting; invariant ones keep S (+) D
    matched = orth(np.column_stack([np.r_[np.eye(m)[:, j], np.zeros(n)] for j in range(m)] +
                                   [np.r_[np.zeros(m), A[:, j]] for j in range(m)]), rcond=TOL)
    for _ in range(3):
        X_psi, X_phi = even_invariant(A, rc(m, m), rc(m, m), rc(n, n))
        X = block(X_psi, X_phi)
        assert preserves(X, R) and preserves(X, Kq) and np.allclose(X @ Q, Q @ X, atol=1e-9)
        assert preserves(X, matched)
        Y = block(X_psi + 0.3 * rc(m, m), X_phi)
        assert not preserves(Y, R) and not np.allclose(Y @ Q, Q @ Y, atol=1e-6)
    # sector-mixing Kraus operators: exact events protected, the bare input is not
    K1 = np.outer(g0, zm.conj())
    K2 = np.eye(m + n) - np.outer(zm, zm.conj())
    assert np.allclose(K1.conj().T @ K1 + K2.conj().T @ K2, np.eye(m + n))
    assert all(preserves(X, R) and preserves(X, Kq) for X in (K1, K2))
    assert not np.allclose(K2 @ Q, Q @ K2)                                  # acts on the quotient without commuting
    weight = lambda u: sum(abs(np.vdot(g0, X @ u)) ** 2 for X in (K1, K2))
    assert np.isclose(weight(zp), 0) and np.allclose(K2 @ zp, zp)
    assert np.isclose(weight(f0), 0.5)
    # a channel from inputs to appearances, extended by zero, fails the premise under R4
    L = np.zeros((n, m))
    L[0, 1] = 1                                                             # f1 -> g0
    X = np.zeros((m + n, m + n))
    X[m:, :m] = L
    assert not preserves(X, R)


def test_carrier_implementations_and_the_defect_bound():
    """Section 4, Interventions. The reflection R + P_N gives the same output as the zero
    extension of V on bare inputs and fixes exact events. A trace-preserving operation that
    preserves exact events sends at most half of a bare input into the defect, and none when the
    whole output lies in the appearance sector; the sector-mixing example attains one half. The
    reset channel |g1><f_j| has a joint row that is not closable."""
    A = shift()
    n, m = A.shape
    d = m + n
    J = grading(m, n)
    Q, PG = charge(A)
    V, PS, PD, _ = polar_data(A)
    Rpol = np.block([[np.zeros((m, m)), V.conj().T], [V, np.zeros((n, n))]])
    PN = np.eye(d) - block(PS, PD)
    G = orth(PG, rcond=TOL)
    Gm = orth(J @ PG @ J, rcond=TOL)
    Nn = orth(PN, rcond=TOL)
    kerQ = null_space(Q, rcond=TOL)

    def preserves(X, basis):
        return np.allclose((np.eye(X.shape[0]) - basis @ basis.conj().T) @ X @ basis, 0, atol=1e-9)

    # two carrier implementations of the same action on bare inputs
    swap = Rpol + PN
    zero_ext = np.zeros((d, d))
    zero_ext[m:, :m] = V
    assert np.allclose(swap.conj().T @ swap, np.eye(d))
    for _ in range(3):
        s = PS @ rc(m)
        b = np.r_[s, np.zeros(n)]
        assert np.allclose(swap @ b, zero_ext @ b)
    assert np.allclose(swap @ G, G)                                          # fixes every exact event
    assert preserves(swap, G) and preserves(swap, kerQ) and not np.allclose(swap @ Q, Q @ swap)
    assert not preserves(zero_ext, G)
    # the bound p <= 1/2 for random exact-preserving operations (Stinespring isometries)
    e = 3
    Eid = np.eye(e)
    defect = block(np.zeros((m, m)), np.eye(n) - PD)
    for _ in range(100):
        cols = []
        for space, k in ((np.kron(G, Eid), G.shape[1]),
                         (np.kron(np.column_stack([G, Nn]), Eid), Nn.shape[1]),
                         (np.eye(d * e), Gm.shape[1])):
            for _ in range(k):
                v = space @ rc(space.shape[1])
                for c in cols:
                    v = v - np.vdot(c, v) * c
                cols.append(v / np.linalg.norm(v))
        Ust = np.column_stack(cols) @ np.column_stack([G, Nn, Gm]).conj().T
        assert np.allclose(Ust.conj().T @ Ust, np.eye(d), atol=1e-9)
        kraus = [Ust.reshape(d, e, d)[:, j, :] for j in range(e)]
        assert all(preserves(K, G) and preserves(K, kerQ) for K in kraus)
        s = PS @ rc(m)
        s /= np.linalg.norm(s)
        b = np.r_[s, np.zeros(n)]
        assert sum(np.linalg.norm(defect @ K @ b) ** 2 for K in kraus) <= 0.5 + 1e-9
    # an output entirely in the appearance sector: ||(P_G (x) I) U b||^2 = (1 - p)/2, so p = 0
    for _ in range(20):
        h = rc(n, e)
        h /= np.linalg.norm(h)
        Ub = np.vstack([np.zeros((m, e)), h])
        p = np.linalg.norm((np.eye(n) - PD) @ h) ** 2
        assert np.isclose(np.linalg.norm(PG @ Ub) ** 2, (1 - p) / 2)
    # the sector-mixing example attains p = 1/2 and leaves a quarter in the physical sector
    f0 = np.r_[np.eye(m)[0], np.zeros(n)]
    g1 = np.r_[np.zeros(m), np.eye(n)[1]]
    g0 = np.r_[np.zeros(m), np.eye(n)[0]]
    zm = (f0 - g1) / np.sqrt(2)
    K1, K2 = np.outer(g0, zm.conj()), np.eye(d) - np.outer(zm, zm.conj())
    assert np.isclose(sum(np.linalg.norm(defect @ K @ f0) ** 2 for K in (K1, K2)), 0.5)
    physical = block(np.eye(m), np.zeros((n, n)))
    assert np.isclose(sum(np.linalg.norm(physical @ K @ f0) ** 2 for K in (K1, K2)), 0.25)
    # reset channel: Kraus operators |g1><f_j| sum to the identity; the joint row is not closable
    M = 40
    kraus = [np.outer(np.eye(M + 1)[1], np.eye(M)[j]) for j in range(M)]
    assert np.allclose(sum(K.T @ K for K in kraus), np.eye(M))
    row = np.hstack(kraus)
    for k in (4, 16, 36):
        psi = np.zeros(M * M)
        for j in range(k):
            psi[j * M + j] = 1 / k                                           # f_j / k in slot j
        assert np.isclose(np.linalg.norm(psi), 1 / np.sqrt(k))               # tends to zero
        assert np.allclose(row @ psi, np.eye(M + 1)[1])                      # image stays g1


def test_ghost_number_is_half_the_polar_reflection():
    for A in (rc(5, 3), rc(5, 2) @ rc(2, 3), shift()):
        n, m = A.shape
        J = grading(m, n)
        Q, PG = charge(A)
        V, PS, PD, _ = polar_data(A)
        PGm = 0.5 * np.block([[PS, -V.conj().T], [-V, PD]])
        R = np.block([[np.zeros((m, m)), V.conj().T], [V, np.zeros((n, n))]])
        assert np.allclose(R, PG - PGm, atol=TOL)
        Ngh = R / 2
        assert np.allclose(Ngh @ Q - Q @ Ngh, Q, atol=TOL)
        assert np.allclose(J @ Ngh.conj().T @ J, -Ngh, atol=TOL)
        assert not np.allclose(J @ Q @ J, -Q) and not np.allclose(J @ Q @ J, Q)   # J does not grade Q
        # the relative phase is read by R: mean cos(phi) on graph states
        s = PS @ rc(m)
        s /= np.linalg.norm(s)
        for phi in (0.0, 1.1, np.pi):
            g = np.r_[s, np.exp(1j * phi) * V @ s] / np.sqrt(2)
            assert np.isclose(np.vdot(g, R @ g).real, np.cos(phi))


def test_identification_keeps_one_copy_and_cohomology_keeps_none():
    A = rc(5, 3)
    n, m = A.shape
    V, PS, PD, r = polar_data(A)
    PGm = 0.5 * np.block([[PS, -V.conj().T], [-V, PD]])
    active = np.block([[PS, np.zeros((m, n))], [np.zeros((n, m)), PD]])
    assert np.linalg.matrix_rank(active, tol=TOL) - np.linalg.matrix_rank(PGm, tol=TOL) == r   # one copy
    Q, _ = charge(A)
    Hc, _, _ = cohomology(Q, grading(m, n))
    assert np.allclose(active @ Hc, 0, atol=1e-9)                                             # none survive
    s = PS @ rc(m)
    s /= np.linalg.norm(s)
    E0 = rc(m, m)
    E = PS @ (E0 @ E0.conj().T) @ PS
    E /= np.linalg.eigvalsh(E).max()
    g = np.r_[s, V @ s] / np.sqrt(2)
    X = block(E, V @ E @ V.conj().T)
    assert np.isclose(np.vdot(g, X @ g), np.vdot(s, E @ s))                                   # matched probability


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print("passed", t.__name__)
    print(len(tests), "passed")
