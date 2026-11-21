"""Checks for the note "Several presentations: redundancy and the reduced space".

Run:  python test_brst_families.py      (or with pytest)
"""
import numpy as np
from scipy.linalg import null_space, orth

RNG = np.random.default_rng(20261122)
TOL = 1e-9


def rc(*shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def polar(d):
    U, s, Vh = np.linalg.svd(d)
    r = int((s > TOL).sum())
    return U[:, :r] @ Vh[:r], Vh[:r].conj().T @ Vh[:r], U[:, :r] @ U[:, :r].conj().T


def chain_charge(dims, maps, degrees):
    """Q = sum_k P_{G_k} J for the chain C_0 -> C_1 -> ... with J = (-1)^degree, and the ghost
    number N = (1/2) sum_k R_k built from the polar reflections of the maps."""
    off = np.cumsum([0] + list(dims))
    n = off[-1]
    J = np.diag(np.concatenate([[(-1.0) ** degrees[i]] * dims[i] for i in range(len(dims))]))
    PG = np.zeros((n, n), complex)
    R = np.zeros((n, n), complex)
    for i, d in enumerate(maps):
        V, PS, PD = polar(d)
        a, b = slice(off[i], off[i + 1]), slice(off[i + 1], off[i + 2])
        PG[a, a] += PS / 2
        PG[a, b] += V.conj().T / 2
        PG[b, a] += V / 2
        PG[b, b] += PD / 2
        R[a, b] += V.conj().T
        R[b, a] += V
    return PG @ J, J, R / 2


def cohomology(Q, J):
    K = null_space(Q, rcond=TOL)
    Rg = orth(Q, rcond=TOL)
    H = K @ null_space(Rg.conj().T @ K, rcond=TOL) if Rg.size else K
    w = np.linalg.eigvalsh(H.conj().T @ J @ H) if H.size else np.array([])
    return H, int((w > TOL).sum()), int((w < -TOL).sum())


def chain_cohomology_dims(dims, maps):
    out = []
    for i in range(len(dims)):
        ker = dims[i] - (np.linalg.matrix_rank(maps[i], tol=TOL) if i < len(maps) else 0)
        im = np.linalg.matrix_rank(maps[i - 1], tol=TOL) if i > 0 else 0
        out.append(ker - im)
    return out


def random_chain(degrees):
    """C_k = A_k (+) H_k (+) B_k with d_k : B_k -> A_{k+1} invertible, in a random basis."""
    L = len(degrees)
    h = [int(RNG.integers(0, 3)) for _ in range(L)]
    b = [int(RNG.integers(1, 3)) for _ in range(L - 1)] + [0]
    a = [0] + b[:-1]
    dims = [a[k] + h[k] + b[k] for k in range(L)]
    bases = [np.linalg.qr(rc(n, n))[0] for n in dims]
    maps = []
    for k in range(L - 1):
        src = bases[k][:, a[k] + h[k]:]
        dst = bases[k + 1][:, :a[k + 1]]
        maps.append(dst @ rc(a[k + 1], b[k]) @ src.conj().T)
    return dims, maps, h


# ---- Theorem 1 -------------------------------------------------------------------------

def test_chain_charge_theorem():
    for degrees in ([0, 1], [-1, 0, 1], [-2, -1, 0, 1], [-1, 0, 1, 2]):
        for _ in range(4):
            dims, maps, h = random_chain(degrees)
            for k in range(len(maps) - 1):
                assert np.allclose(maps[k + 1] @ maps[k], 0)
            Q, J, N = chain_charge(dims, maps, degrees)
            assert np.allclose(Q @ Q, 0, atol=TOL)
            assert np.allclose(J @ Q.conj().T @ J, Q, atol=TOL)
            assert chain_cohomology_dims(dims, maps) == h
            _, plus, minus = cohomology(Q, J)
            assert plus == sum(h[k] for k in range(len(h)) if degrees[k] % 2 == 0)
            assert minus == sum(h[k] for k in range(len(h)) if degrees[k] % 2 != 0)
            assert np.allclose(N @ Q - Q @ N, Q, atol=TOL)
            assert np.allclose(J @ N.conj().T @ J, -N, atol=TOL)


def test_two_term_chain_is_the_single_presentation_charge():
    A = rc(5, 2) @ rc(2, 3)
    Q, J, _ = chain_charge([3, 5], [A], [0, 1])
    _, plus, minus = cohomology(Q, J)
    assert (plus, minus) == (3 - 2, 5 - 2)


# ---- Corollary 3: resolving redundancy ---------------------------------------------------

def test_resolution_leaves_exactly_the_common_excess():
    for k_members, m, n, rank in ((2, 3, 5, 2), (3, 2, 5, 2), (2, 4, 6, 3)):
        members = [rc(n, rank) @ rc(rank, m) for _ in range(k_members)]
        B = np.hstack(members)
        K = null_space(B, rcond=TOL)
        defect = n - np.linalg.matrix_rank(B, tol=TOL)
        # single charge: kernel survives with positive sign, excess with negative sign
        Q1, J1, _ = chain_charge([B.shape[1], n], [B], [0, 1])
        assert cohomology(Q1, J1)[1:] == (K.shape[1], defect)
        # resolved: the redundancies pair with a space of degree -1
        Q, J, _ = chain_charge([K.shape[1], B.shape[1], n], [K, B], [-1, 0, 1])
        H, plus, minus = cohomology(Q, J)
        assert (plus, minus) == (0, defect)
        # the surviving directions are exactly the common defect of the family
        common = null_space(np.vstack([M.conj().T for M in members]), rcond=TOL)
        off = K.shape[1] + B.shape[1]
        assert np.allclose(H[:off], 0, atol=1e-8)
        if defect:
            P = common @ common.conj().T
            assert np.allclose(P @ H[off:], H[off:], atol=1e-8)


# ---- Proposition 4: redundancy explained by shared appearances (three presentations) --------

def intersection(*bases):
    """Orthonormal basis of the intersection of the column spaces, with an absolute tolerance.
    The complementary projectors are stacked, so a vector at angle t from the intersection has
    singular value of order t rather than t^2 (summing them would give t^2 / 2, and two lines
    1e-5 apart would be merged). A relative tolerance fails when the stack is pure roundoff, for
    example for two copies of the full space."""
    n = bases[0].shape[0]
    M = np.vstack([np.eye(n) - B @ B.conj().T for B in bases])
    _, s, Vh = np.linalg.svd(M)
    s = np.r_[s, np.zeros(n - len(s))]
    return Vh.conj().T[:, s < 1e-8]


def cech_three(Ds, n):
    """Chain  (D1 n D2 n D3) -> (+)_{pairs} (Di n Dj) -> (+)_c D_c -> C^n  in degrees -2..1."""
    pairs = [(0, 1), (0, 2), (1, 2)]
    Ip = [intersection(Ds[i], Ds[j]) for i, j in pairs]
    It = intersection(*Ds)
    dD = [D.shape[1] for D in Ds]
    dP = [x.shape[1] for x in Ip]
    oD = np.cumsum([0] + dD)
    oP = np.cumsum([0] + dP)
    sigma = np.hstack(Ds)
    delta = np.zeros((oD[-1], oP[-1]), complex)
    for q, (i, j) in enumerate(pairs):
        w = Ip[q]
        delta[oD[i]:oD[i + 1], oP[q]:oP[q + 1]] = Ds[i].conj().T @ w
        delta[oD[j]:oD[j + 1], oP[q]:oP[q + 1]] = -Ds[j].conj().T @ w
    eps = np.zeros((oP[-1], It.shape[1]), complex)
    for (a, b), sgn in (((1, 2), 1), ((0, 2), -1), ((0, 1), 1)):
        p = pairs.index((a, b))
        eps[oP[p]:oP[p + 1], :] = sgn * Ip[p].conj().T @ It
    assert np.allclose(sigma @ delta, 0) and np.allclose(delta @ eps, 0)
    dims = [It.shape[1], oP[-1], oD[-1], n]
    maps = [eps, delta, sigma]
    return dims, maps, [-2, -1, 0, 1], dD, dP, It.shape[1]


def collective_redundancy(Ds, n):
    dims, maps, degrees, dD, dP, dT = cech_three(Ds, n)
    span = np.linalg.matrix_rank(np.hstack(Ds), tol=TOL)
    return sum(dD) - sum(dP) + dT - span, n - span


def line(*v):
    v = np.asarray(v, complex)
    return (v / np.linalg.norm(v)).reshape(-1, 1)


def test_cech_cohomology_for_three_presentations():
    angles = (0, np.pi / 3, 2 * np.pi / 3)
    cases = {
        "three lines in C^2": ([line(np.cos(t), np.sin(t)) for t in angles], 2, (1, 0)),
        "three lines in a plane of C^3": ([line(np.cos(t), np.sin(t), 0) for t in angles], 3, (1, 1)),
        "three random planes in C^3": ([orth(rc(3, 2)) for _ in range(3)], 3, (0, 0)),
        "three random lines in C^3": ([orth(rc(3, 1)) for _ in range(3)], 3, (0, 0)),
        "a line inside a plane, and a third line": ([line(1, 0, 0), orth(np.array([[1, 0], [0, 1], [0, 0]], complex)),
                                                     line(1, 1, 1)], 3, (0, 0)),
    }
    for label, (Ds, n, expected) in cases.items():
        dims, maps, degrees, dD, dP, dT = cech_three(Ds, n)
        hd = chain_cohomology_dims(dims, maps)
        assert hd[0] == 0 and hd[1] == 0, label                       # degrees -2 and -1 always vanish
        assert (hd[2], hd[3]) == expected, label                      # collective redundancy, excess
        assert collective_redundancy(Ds, n) == expected, label        # inclusion-exclusion formula
        Q, J, _ = chain_charge(dims, maps, degrees)
        _, plus, minus = cohomology(Q, J)
        assert (plus, minus) == expected, label                        # signs: + collective, - excess


def test_two_presentations_have_no_collective_redundancy():
    for _ in range(5):
        n = 5
        D1, D2 = orth(rc(n, int(RNG.integers(1, 4)))), orth(rc(n, int(RNG.integers(1, 4))))
        I = intersection(D1, D2)
        sigma = np.hstack([D1, D2])
        delta = np.vstack([D1.conj().T @ I, -D2.conj().T @ I])
        assert np.linalg.matrix_rank(null_space(sigma, rcond=TOL), tol=TOL) == I.shape[1]
        assert np.allclose(sigma @ delta, 0)


# ---- the shift family --------------------------------------------------------------------

def shift(N=6):
    A = np.zeros((N + 1, N))
    for k in range(N):
        A[k + 1, k] = 1
    return A


def hamiltonians(eps, N=6):
    H_psi = np.diag([0.0] + [1.0] * (N - 1))
    H_phi = np.diag([eps, 0.0] + [1.0] * (N - 1))
    return H_psi, H_phi


def test_shift_family_resolved_keeps_the_missing_line_and_invariant_dynamics():
    N = 6
    A = shift(N)
    n, m = A.shape
    C = np.zeros((n, m), complex)                                     # a second intertwiner
    C[1, 0] = np.exp(0.3j)
    C[2:, 1:] = np.linalg.qr(rc(n - 2, m - 1))[0][:, :m - 1]
    B = np.hstack([A, C])
    K = null_space(B, rcond=TOL)
    assert K.shape[1] > 0
    Q, J, _ = chain_charge([K.shape[1], 2 * m, n], [K, B], [-1, 0, 1])
    H, plus, minus = cohomology(Q, J)
    assert (plus, minus) == (0, 1)
    g0 = np.zeros(K.shape[1] + 2 * m + n)
    g0[K.shape[1] + 2 * m] = 1
    assert np.allclose(H @ H.conj().T @ g0, g0, atol=1e-8)
    for eps in (2.5, 0.4):
        H_psi, H_phi = hamiltonians(eps, N)
        Hm = np.kron(np.eye(2), H_psi)
        assert np.allclose(H_phi @ B, B @ Hm)                           # the joint map intertwines
        Hr = K.conj().T @ Hm @ K                                        # ker B is invariant
        assert np.allclose(Hm @ K, K @ Hr, atol=1e-9)
        total = np.zeros_like(Q)
        r = K.shape[1]
        total[:r, :r] = Hr
        total[r:r + 2 * m, r:r + 2 * m] = Hm
        total[r + 2 * m:, r + 2 * m:] = H_phi
        assert np.allclose(total @ Q, Q @ total, atol=1e-9)


def test_intertwining_chains_have_invariant_charges():
    """Remark after Theorem 1: if every map of a chain intertwines the given dynamics, the chain
    charge commutes with it. Checked on the shared-appearance chain of three invariant ranges."""
    n = 5
    H_phi = np.diag([0.0, 0.0, 1.0, 1.0, 1.0])
    e = np.eye(n)
    Ds = [orth(np.column_stack([e[0], e[2]])),
          orth(np.column_stack([e[0] + e[1], e[3]])),
          orth(np.column_stack([e[1], e[2] + e[3] + e[4]]))]
    for D in Ds:                                                       # each range reduces H_phi
        P = D @ D.conj().T
        assert np.allclose(P @ H_phi, H_phi @ P)
    dims, maps, degrees, dD, dP, dT = cech_three(Ds, n)
    Q, J, _ = chain_charge(dims, maps, degrees)
    # restricted dynamics on every space of the chain
    pairs = [(0, 1), (0, 2), (1, 2)]
    Ip = [intersection(Ds[i], Ds[j]) for i, j in pairs]
    It = intersection(*Ds)
    blocks = [It.conj().T @ H_phi @ It] + [W.conj().T @ H_phi @ W for W in Ip] + \
             [D.conj().T @ H_phi @ D for D in Ds] + [H_phi]
    total = np.zeros_like(Q)
    i = 0
    for b in blocks:
        k = b.shape[0]
        total[i:i + k, i:i + k] = b
        i += k
    assert i == Q.shape[0]
    off = np.cumsum([0] + dims)
    for k, d in enumerate(maps):                                       # each map intertwines
        Hs = total[off[k]:off[k + 1], off[k]:off[k + 1]]
        Ht = total[off[k + 1]:off[k + 2], off[k + 1]:off[k + 2]]
        assert np.allclose(Ht @ d, d @ Hs, atol=1e-9)
    assert np.allclose(total @ Q, Q @ total, atol=1e-9)


def test_intersection_helper_handles_identical_and_full_spaces():
    full = np.linalg.qr(rc(3, 3))[0]
    assert intersection(full, full).shape[1] == 3
    sub = orth(rc(5, 2))
    assert intersection(sub, sub).shape[1] == 2
    assert intersection(sub, sub, sub).shape[1] == 2
    other = orth(rc(5, 2))
    assert intersection(sub, other).shape[1] == 0
    line_in_plane = orth(np.column_stack([sub[:, 0]]))
    assert intersection(line_in_plane, sub).shape[1] == 1


def test_intersection_helper_separates_nearby_lines():
    """Two distinct lines 1e-5 radians apart meet only in zero, and the three-presentation chain
    built with a third line passes its own checks."""
    t = 1e-5
    L1, L2 = line(1, 0, 0), line(np.cos(t), np.sin(t), 0)
    assert intersection(L1, L2).shape[1] == 0
    assert intersection(L1, L1).shape[1] == 1
    for L3, expected in ((line(0, 0, 1), (0, 0)), (line(0, 1, 0), (1, 1))):
        Ds = [L1, L2, L3]
        dims, maps, degrees, dD, dP, dT = cech_three(Ds, 3)
        assert dP == [0, 0, 0] and dT == 0
        assert chain_cohomology_dims(dims, maps)[2:] == list(expected)
        assert collective_redundancy(Ds, 3) == expected
        Q, J, _ = chain_charge(dims, maps, degrees)
        assert cohomology(Q, J)[1:] == expected


def test_bare_redundancy_has_no_class_and_its_aligned_pair_is_exact():
    """Section 5: a redundant input k is a boundary of the resolved chain. In the graph charge the
    bare vector (0, k, 0) is not closed, Q(0, k, 0) = (k, k, 0)/2, while the aligned pair
    (k, k, 0)/sqrt2 is exact."""
    N = 6
    A = shift(N)
    n, m = A.shape
    C = np.zeros((n, m), complex)
    C[1, 0] = 1
    C[2:, 1:] = np.linalg.qr(rc(n - 2, m - 1))[0][:, :m - 1]
    B = np.hstack([A, C])
    K = null_space(B, rcond=TOL)
    r = K.shape[1]
    Q, J, _ = chain_charge([r, 2 * m, n], [K, B], [-1, 0, 1])
    exact = orth(Q, rcond=TOL)
    for _ in range(3):
        k = K @ rc(r)
        lift = K.conj().T @ k                                         # k as an element of degree -1
        bare = np.r_[np.zeros(r), k, np.zeros(n)]
        pair = np.r_[lift, k, np.zeros(n)] / np.sqrt(2)
        assert np.allclose(Q @ bare, np.r_[lift, k, np.zeros(n)] / 2, atol=1e-9)
        assert not np.allclose(Q @ bare, 0, atol=1e-6)                 # not closed: no class
        assert np.allclose(exact @ exact.conj().T @ pair, pair, atol=1e-9)   # exact: zero class


def test_collective_redundancy_with_a_shared_ground():
    """Every pair shares g, yet three distinct lines in a plane orthogonal to g give one
    cancellation that the shared g cannot explain."""
    g = np.array([1, 0, 0], complex)
    lines = [np.array([0, np.cos(t), np.sin(t)], complex) for t in (0.0, np.pi / 3, 2 * np.pi / 3)]
    Ds = [orth(np.column_stack([g, L])) for L in lines]
    assert all(intersection(Ds[i], Ds[j]).shape[1] == 1 for i, j in ((0, 1), (0, 2), (1, 2)))
    assert intersection(*Ds).shape[1] == 1
    assert collective_redundancy(Ds, 3) == (1, 0)
    dims, maps, degrees, dD, dP, dT = cech_three(Ds, 3)
    assert chain_cohomology_dims(dims, maps) == [0, 0, 1, 0]
    Q, J, _ = chain_charge(dims, maps, degrees)
    assert cohomology(Q, J)[1:] == (1, 0)


def rotated_shift(alpha, N=6):
    U = np.eye(N + 1)
    c, s = np.cos(alpha), np.sin(alpha)
    U[0, 0], U[2, 0], U[0, 2], U[2, 2] = c, s, -s, c
    return U @ shift(N)


def test_shift_families_at_resonance():
    """At eps = 1 a second intertwiner can prepare g0. Families then have no common defect, and
    the resolved reduced space is zero; the single joint charge keeps the redundancies."""
    N = 6
    A = shift(N)
    n, m = A.shape
    H_psi, H_phi = hamiltonians(1.0, N)
    half_turn = rotated_shift(np.pi / 2, N)
    assert np.allclose(half_turn @ np.eye(m)[1], -np.eye(n)[0])           # f1 -> -g0
    for alpha in (0.7, np.pi / 2):
        B = np.hstack([A, rotated_shift(alpha, N)])
        assert np.allclose(H_phi @ B, B @ np.kron(np.eye(2), H_psi))       # both intertwine at eps = 1
        assert np.linalg.matrix_rank(B, tol=TOL) == n                      # no common defect
        Q1, J1, _ = chain_charge([2 * m, n], [B], [0, 1])
        kernel_dim = 2 * m - n
        assert cohomology(Q1, J1)[1:] == (kernel_dim, 0)                  # redundancies remain
        K = null_space(B, rcond=TOL)
        Q, J, _ = chain_charge([K.shape[1], 2 * m, n], [K, B], [-1, 0, 1])
        assert cohomology(Q, J)[1:] == (0, 0)                              # resolved: zero


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print("passed", t.__name__)
    print(len(tests), "passed")
