"""Four spin-1/2 fermions in a two-dimensional harmonic trap with Gaussian attraction
(LG, sec:orbital-ground-source), in a shell-truncated oscillator basis.

Units hbar = M = omega = 1, so the oscillator length is 1 and xi = a^2:

    H = H_ho - lam * W - Omega_B * L_tot,
    W = (1 / (2 pi xi)) * sum_{i<j} exp(-|r_i - r_j|^2 / (2 xi))   (a positive operator).

The basis keeps every Cartesian oscillator orbital with n_x + n_y <= Nmax.  Each shell
is rotation invariant, so the truncation preserves rotation, spin and time-reversal
symmetry exactly.  Operators act on the S_z = 0 sector (two spin-up, two spin-down);
every spin multiplet has a member there, so its spectrum contains every energy level.
"""
from functools import lru_cache
from itertools import combinations
from math import factorial, pi, sqrt

import numpy as np
import scipy.sparse as sp
from numpy.polynomial.hermite import hermgauss, hermval
from scipy.sparse.linalg import eigsh


def gaussian_integrals_1d(Nmax, xi):
    """I[a, b, c, d] = int phi_a(x1) phi_b(x2) exp(-(x1 - x2)^2 / (2 xi)) phi_c(x1) phi_d(x2).

    Exact Gauss-Hermite quadrature in u = (x1 + x2)/sqrt2 and w = (x1 - x2)/sqrt2, where
    the integrand is a polynomial times exp(-u^2) exp(-w^2 (1 + 1/xi))."""
    K = 2 * Nmax + 8
    t, wt = hermgauss(K)
    scale = 1 / np.sqrt(1 + 1 / xi)
    U, Wp = np.meshgrid(t, t * scale, indexing="ij")
    weight = np.outer(wt, wt) * scale
    x1, x2 = (U + Wp) / sqrt(2), (U - Wp) / sqrt(2)

    def poly(k, x):
        c = np.zeros(k + 1)
        c[k] = 1
        return hermval(x, c) / sqrt(2.0**k * factorial(k) * sqrt(pi))

    P1 = [poly(k, x1) for k in range(Nmax + 1)]
    P2 = [poly(k, x2) for k in range(Nmax + 1)]
    out = np.zeros((Nmax + 1,) * 4)
    for a in range(Nmax + 1):
        for b in range(Nmax + 1):
            for c in range(Nmax + 1):
                for d in range(Nmax + 1):
                    if (a + b + c + d) % 2 == 0:
                        out[a, b, c, d] = np.sum(weight * P1[a] * P1[c] * P2[b] * P2[d])
    return out


@lru_cache(maxsize=None)
def build(Nmax, xi):
    """Sparse operators on the S_z = 0 sector: Hho, W, L, Qx, S2 and the dimension."""
    orbitals = [(nx, N - nx) for N in range(Nmax + 1) for nx in range(N + 1)]
    n = len(orbitals)
    index = {o: k for k, o in enumerate(orbitals)}
    I1 = gaussian_integrals_1d(Nmax, xi)
    w2 = np.zeros((n,) * 4)                                   # <ab|W|cd>, a,c particle 1
    for a, (ax, ay) in enumerate(orbitals):
        for b, (bx, by) in enumerate(orbitals):
            for c, (cx, cy) in enumerate(orbitals):
                for d, (dx, dy) in enumerate(orbitals):
                    w2[a, b, c, d] = I1[ax, bx, cx, dx] * I1[ay, by, cy, dy] / (2 * pi * xi)

    def x2_element(k, kp):                                    # <kp| x^2 |k> in one dimension
        if kp == k:
            return k + 0.5
        if kp == k + 2:
            return 0.5 * sqrt((k + 1) * (k + 2))
        if kp == k - 2:
            return 0.5 * sqrt(k * (k - 1))
        return 0.0

    energy = np.array([nx + ny + 1.0 for nx, ny in orbitals])
    lz = np.zeros((n, n), complex)                            # L = x p_y - y p_x
    x2m, y2m = np.zeros((n, n)), np.zeros((n, n))
    for k, (nx, ny) in enumerate(orbitals):
        if nx >= 1 and (nx - 1, ny + 1) in index:
            lz[index[(nx - 1, ny + 1)], k] += 1j * sqrt(nx * (ny + 1))
        if ny >= 1 and (nx + 1, ny - 1) in index:
            lz[index[(nx + 1, ny - 1)], k] -= 1j * sqrt((nx + 1) * ny)
        for kp, (mx, my) in enumerate(orbitals):
            if my == ny:
                x2m[kp, k] = x2_element(nx, mx)
            if mx == nx:
                y2m[kp, k] = x2_element(ny, my)

    # two same-spin fermions: determinants |r s> = a_r^* a_s^* |0>, r < s
    dets = list(combinations(range(n), 2))
    det_index = {dt: k for k, dt in enumerate(dets)}
    D = len(dets)

    def hop(a, c):                                            # a_a^* a_c on those determinants
        rows, cols, vals = [], [], []
        for k, (r, s) in enumerate(dets):
            if c == r:
                other, sign = s, 1.0
            elif c == s:
                other, sign = r, -1.0
            else:
                continue
            if a == other:
                continue
            new = (a, other) if a < other else (other, a)
            rows.append(det_index[new])
            cols.append(k)
            vals.append(sign if a < other else -sign)
        return sp.csr_matrix((vals, (rows, cols)), shape=(D, D))

    hops = [[hop(a, c) for c in range(n)] for a in range(n)]

    def one_body(o):
        out = sp.csr_matrix((D, D), dtype=complex if np.iscomplexobj(o) else float)
        for a in range(n):
            for c in range(n):
                if o[a, c] != 0:
                    out = out + o[a, c] * hops[a][c]
        return out

    eye = sp.identity(D, format="csr")

    def both_spins(op):
        return sp.kron(op, eye, format="csr") + sp.kron(eye, op, format="csr")

    same = np.zeros((D, D))                                   # direct minus exchange
    for k, (p, q) in enumerate(dets):
        for l, (r, s) in enumerate(dets):
            same[k, l] = w2[p, q, r, s] - w2[p, q, s, r]
    W = both_spins(sp.csr_matrix(same))
    for a in range(n):                                        # opposite-spin part
        for c in range(n):
            partner = one_body(w2[a, :, c, :])
            if partner.nnz:
                W = W + sp.kron(hops[a][c], partner, format="csr")

    # S^2 on the S_z = 0 sector as S_+^* S_+, with S_+ mapping into (3 up, 1 down)
    triples = list(combinations(range(n), 3))
    triple_index = {dt: k for k, dt in enumerate(triples)}
    rows, cols, vals = [], [], []
    for ku, (u1, u2) in enumerate(dets):
        for kd, (d1, d2) in enumerate(dets):
            for o, rest, sign in ((d1, d2, 1.0), (d2, d1, -1.0)):
                if o in (u1, u2):
                    continue
                up = tuple(sorted((u1, u2, o)))
                rows.append(triple_index[up] * n + rest)
                cols.append(ku * D + kd)
                vals.append(sign * (-1.0) ** up.index(o))
    Splus = sp.csr_matrix((vals, (rows, cols)), shape=(len(triples) * n, D * D))

    return {
        "Hho": both_spins(one_body(np.diag(energy))),
        "W": W.tocsr(),
        "L": both_spins(one_body(lz)),
        "Qx": both_spins(one_body((x2m - y2m) / 2)),
        "S2": (Splus.T @ Splus).tocsr(),
        "dim": D * D,
    }


def lowest(H, k):
    """Lowest k eigenpairs (dense for small problems, Lanczos otherwise)."""
    if H.shape[0] <= 2500:
        w, v = np.linalg.eigh(H.toarray())
        return w[:k], v[:, :k]
    w, v = eigsh(H, k=k, which="SA", tol=1e-12)
    order = np.argsort(w)
    return w[order], v[:, order]


def first_order(lam, xi):
    """First-order projected energies of LG eq:orbital-projected-energies (ell = 1)."""
    v = -lam / (8 * pi)
    core = 4 * v * (5 * xi + 2) / (1 + xi) ** 2
    return {
        "g": core + 4 * v * (xi**2 + xi + 1) / (1 + xi) ** 3,
        "e": core + 2 * v * (2 * xi**2 + 2 * xi + 1) / (1 + xi) ** 3,
        "W": core + 4 * v * xi / (1 + xi) ** 2,
        "d_a": 1 / (4 * pi * (1 + xi) ** 3),
    }
