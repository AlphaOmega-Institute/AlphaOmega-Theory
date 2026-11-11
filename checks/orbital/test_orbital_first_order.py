"""First-order (degenerate perturbation theory) checks for the orbital parent:
the projected interaction matrix on the s^2 p^2 space, the angular-momentum
identities and the quadrupole compressions (LG, sec:orbital-ground-source)."""
from itertools import product
from math import pi

import numpy as np

from lib.fock import modes, dag, vacuum, number
from orbital.model import first_order

LABELS = {
    "LG": [
        "sec:orbital-ground-source", "eq:orbital-three-orbitals", "eq:orbital-derived-current",
        "eq:orbital-projected-energies", "eq:orbital-quadrupoles",
    ],
}

# Oscillator orbitals s, p_x, p_y with ell = 1: polynomial (x power, y power, coefficient)
# multiplying exp(-r^2/2)/sqrt(pi).
POLY = {"s": (0, 0, 1.0), "x": (1, 0, np.sqrt(2)), "y": (0, 1, np.sqrt(2))}
ORBITALS = ["s", "x", "y"]
SPIN_ORBITALS = [(o, s) for o in ORBITALS for s in (0, 1)]
OPS = dict(zip(SPIN_ORBITALS, modes(6)))


def _moment(px1, px2, xi):
    """int x1^px1 x2^px2 exp(-x1^2 - x2^2 - (x1 - x2)^2/(2 xi)) by Wick contraction."""
    A = 2 * np.array([[1 + 1 / (2 * xi), -1 / (2 * xi)], [-1 / (2 * xi), 1 + 1 / (2 * xi)]])
    S = np.linalg.inv(A)
    idx = [0] * px1 + [1] * px2
    if len(idx) % 2:
        return 0.0

    def wick(l):
        if not l:
            return 1.0
        first, rest = l[0], l[1:]
        return sum(S[first, rest[k]] * wick(rest[:k] + rest[k + 1:]) for k in range(len(rest)))

    return 2 * pi / np.sqrt(np.linalg.det(A)) * wick(idx)


def _two_body(a, b, c, d, lam, xi):
    pa, pb, pc, pd = POLY[a], POLY[b], POLY[c], POLY[d]
    coefficient = pa[2] * pb[2] * pc[2] * pd[2] / pi**2
    jx = _moment(pa[0] + pc[0], pb[0] + pd[0], xi)
    jy = _moment(pa[1] + pc[1], pb[1] + pd[1], xi)
    return -(lam / (2 * pi * xi)) * coefficient * jx * jy


def _interaction(lam, xi):
    V = np.zeros((64, 64))
    for A_, B_, C_, D_ in product(SPIN_ORBITALS, repeat=4):
        if A_[1] != C_[1] or B_[1] != D_[1]:
            continue
        value = _two_body(A_[0], B_[0], C_[0], D_[0], lam, xi)
        if value:
            V += 0.5 * value * dag(OPS[A_]) @ dag(OPS[B_]) @ OPS[D_] @ OPS[C_]
    return V


def _valence():
    """g, e, m and the triplet W on the s^2 p^2 space; c_i = p_x with spin i, d_i = p_y."""
    core = dag(OPS[("s", 0)]) @ dag(OPS[("s", 1)]) @ vacuum(6)
    c = [OPS[("x", 0)], OPS[("x", 1)]]
    d = [OPS[("y", 0)], OPS[("y", 1)]]
    C = dag(c[0]) @ dag(c[1]) @ core
    D = dag(d[0]) @ dag(d[1]) @ core
    g, e = (C + D) / np.sqrt(2), (C - D) / np.sqrt(2)
    m = (dag(c[0]) @ dag(d[1]) + dag(d[0]) @ dag(c[1])) @ core / np.sqrt(2)
    W = [dag(c[0]) @ dag(d[0]) @ core, dag(c[1]) @ dag(d[1]) @ core,
         (dag(c[0]) @ dag(d[1]) - dag(d[0]) @ dag(c[1])) @ core / np.sqrt(2)]
    return c, d, g, e, m, W


def test_projected_interaction_matrix():
    """The first-order matrix on s^2 p^2 is diagonal in (g, e, m, W) with the displayed
    entries; the gaps above g are exactly lam d_a (singlets) and 2 lam d_a (triplet)."""
    _, _, g, e, m, W = _valence()
    basis = [g, e, m] + W
    lam = 1.0
    for xi in (0.3, 1.0, 2.7):
        V = _interaction(lam, xi)
        matrix = np.array([[np.vdot(u, V @ w) for w in basis] for u in basis])
        formulas = first_order(lam, xi)
        expected = [formulas["g"], formulas["e"], formulas["e"]] + [formulas["W"]] * 3
        assert np.allclose(matrix, np.diag(expected), atol=1e-12)
        assert np.isclose(formulas["e"] - formulas["g"], lam * formulas["d_a"])
        assert np.isclose(formulas["W"] - formulas["g"], 2 * lam * formulas["d_a"])


def test_angular_momentum_on_valence_shell():
    """L_tot = 2 Y on s^2 p^2 (hbar = 1): L g = 0, L e = 2i m, L m = -2i e."""
    c, d, g, e, m, W = _valence()
    # one-body L = x p_y - y p_x gives L|p_x> = i|p_y>, L|p_y> = -i|p_x>
    L = sum(1j * dag(dk) @ ck - 1j * dag(ck) @ dk for ck, dk in zip(c, d))
    T = sum(dag(ck) @ dk for ck, dk in zip(c, d))
    Y = (T - dag(T)) / 2j
    assert np.allclose(L, 2 * Y)
    assert np.allclose(L @ g, 0)
    assert np.allclose(L @ e, 2j * m)
    assert np.allclose(L @ m, -2j * e)
    for w in W:
        assert np.allclose(L @ w, 0)


def test_quadrupole_compressions():
    """On the valence shell the quadrupoles Q_x, Q_y compress to Z and X; in particular
    <e| Q_x |g> = 1 and <m| Q_x |g> = 0."""
    def element(a, b, px, py):                       # <a| x^px y^py |b>
        pa, pb = POLY[a], POLY[b]

        def gauss(k):                                # int x^k exp(-x^2) dx
            return 0.0 if k % 2 else np.sqrt(pi) * np.prod(np.arange(k - 1, 0, -2), initial=1.0) / 2 ** (k / 2)

        return pa[2] * pb[2] / pi * gauss(pa[0] + pb[0] + px) * gauss(pa[1] + pb[1] + py)

    qx = np.array([[0.5 * (element(p, q, 2, 0) - element(p, q, 0, 2)) for q in ORBITALS] for p in ORBITALS])
    qy = np.array([[element(p, q, 1, 1) for q in ORBITALS] for p in ORBITALS])
    assert np.allclose(qx, np.diag([0.0, 0.5, -0.5]))
    assert np.allclose(qy, [[0, 0, 0], [0, 0, 0.5], [0, 0.5, 0]])

    def second_quantized(o):
        return sum(o[i, k] * dag(OPS[(ORBITALS[i], s)]) @ OPS[(ORBITALS[k], s)]
                   for i in range(3) for k in range(3) for s in (0, 1))

    c, d, g, e, m, W = _valence()
    Z = (number(c[0]) + number(c[1]) - number(d[0]) - number(d[1])) / 2
    T = sum(dag(ck) @ dk for ck, dk in zip(c, d))
    X = (T + dag(T)) / 2
    basis = np.column_stack([g, e, m] + W)
    for full, compressed in ((second_quantized(qx), Z), (second_quantized(qy), X)):
        assert np.allclose(basis.conj().T @ full @ basis, basis.conj().T @ compressed @ basis)
    assert np.isclose(e @ second_quantized(qx) @ g, 1.0)
    assert np.isclose(m @ second_quantized(qx) @ g, 0.0)
