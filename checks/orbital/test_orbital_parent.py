"""Checks of the orbital parent beyond first order (LG, prop:orbital-full-ground and
the selected-carrier construction), by diagonalizing the four-fermion Hamiltonian in a
shell-truncated oscillator basis.  Units hbar = M = omega = 1.

The standard run uses shells n_x + n_y <= 3 (dimension 2025).  `--full` repeats the
spectral checks with shells <= 4 (dimension 11025) to show the numbers are converged.
"""
import os
from math import pi, sqrt

import numpy as np

from lib import skip
from orbital.model import build, lowest, first_order

LABELS = {
    "LG": [
        "eq:orbital-parent", "eq:orbital-bias-parameters", "eq:orbital-minimal-coupling",
        "prop:orbital-full-ground", "eq:orbital-low-carrier", "eq:orbital-typed-preparation",
        "eq:orbital-filtered-readout", "eq:orbital-low-response", "eq:orbital-unfiltered-response",
    ],
}

XI = 1.0
LAM = 0.2


def _zero_bias(Nmax, lam):
    ops = build(Nmax, XI)
    H0 = (ops["Hho"] - lam * ops["W"]).astype(complex)
    w, v = lowest(H0, 6)
    return ops, H0, w, v


def _spectral_structure(Nmax):
    d_a = first_order(1.0, XI)["d_a"]
    deviations = []
    for lam in (0.05, 0.2):
        ops, H0, w, v = _zero_bias(Nmax, lam)
        L2 = [np.linalg.norm(ops["L"] @ v[:, k]) ** 2 for k in range(4)]
        S2 = [np.vdot(v[:, k], ops["S2"] @ v[:, k]).real for k in range(4)]
        # ground: L = 0 singlet; then the L = +-2 singlets (degenerate); then the triplet
        assert np.allclose(L2, [0, 4, 4, 0], atol=1e-6)
        assert np.allclose(S2, [0, 0, 0, 2], atol=1e-6)
        assert np.isclose(w[1], w[2], atol=1e-9)
        gaps = (w[1:4] - w[0]) / (lam * d_a)
        deviations.append(abs(gaps[0] - 1))
        assert abs(gaps[0] - 1) < 0.1 * lam and abs(gaps[2] - 2) < 0.2 * lam
        # ground energy: first order plus an O(lam^2) remainder
        assert abs(w[0] - 6 - first_order(lam, XI)["g"]) < 0.02 * lam**2
        # the next level is one trap quantum up, as the min-max bound requires
        norm_W = 6 / (2 * pi * XI)
        assert w[4] - 6 >= 1 - lam * norm_W - 1e-9
        assert w[4] - w[0] > 0.9
    assert deviations[0] < deviations[1]                      # the correction shrinks with lam


def test_ground_and_excitation_structure():
    """Unique L = 0 singlet ground; L = +-2 singlets at gap -> lam d_a and the triplet at
    gap -> 2 lam d_a as lam -> 0; the seventh level stays a trap quantum above."""
    _spectral_structure(3)


def test_bias_shifts_and_window():
    """With H = H_0 - Omega_B L and Omega_B = chi lam d_a / 2: the ground energy and vector
    are unchanged and the singlet gaps are (E_2 - E_0) -+ 2 Omega_B, for |chi| < 1.
    Outside the window (chi = 1.2) the L = 2 state becomes the ground."""
    ops, H0, w0, _ = _zero_bias(3, LAM)
    d_a = first_order(1.0, XI)["d_a"]
    gap = w0[1] - w0[0]
    for chi in (0.5, 0.95):
        Om = chi * LAM * d_a / 2
        w, v = lowest(H0 - Om * ops["L"], 4)
        assert abs(w[0] - w0[0]) < 1e-9
        assert np.linalg.norm(ops["L"] @ v[:, 0]) < 1e-7
        assert np.allclose(w[1:3] - w[0], [gap - 2 * Om, gap + 2 * Om], atol=1e-9)
        assert w[1] - w[0] > 0
    Om = 1.2 * LAM * d_a / 2
    w, v = lowest(H0 - Om * ops["L"], 2)
    assert w[0] < w0[0] - 1e-6
    assert np.isclose(np.linalg.norm(ops["L"] @ v[:, 0]), 2.0, atol=1e-6)


def _carrier(chi=0.5):
    ops, H0, w0, _ = _zero_bias(3, LAM)
    Om = chi * LAM * first_order(1.0, XI)["d_a"] / 2
    HB = H0 - Om * ops["L"]
    w, v = lowest(HB, 4)
    Lz = [np.vdot(v[:, k], ops["L"] @ v[:, k]).real for k in range(4)]
    G = v[:, 0]
    Fp = v[:, int(np.argmin([abs(x - 2) for x in Lz]))]
    Fm = v[:, int(np.argmin([abs(x + 2) for x in Lz]))]
    QG = ops["Qx"] @ G
    Fp = Fp * np.exp(1j * np.angle(np.vdot(Fp, QG)))          # make <F|Qx|G> real positive
    Fm = Fm * np.exp(1j * np.angle(np.vdot(Fm, QG)))
    e = (Fp + Fm) / sqrt(2)
    m = (Fp - Fm) / (1j * sqrt(2))
    return ops, H0, HB, Om, w, G, e, m, QG


def test_filtered_quadrupole_readout():
    """On span{G, e, m}: P Q_x P / ||P Q_x G|| = |G><e| + |e><G|, so Z^2 = AA^* with
    A|0> = G, A|1> = e; H_Psi = A^* H A is diagonal; Ker A^* is the line through m."""
    ops, H0, HB, Om, w, G, e, m, QG = _carrier()
    B = np.column_stack([G, e, m])
    assert np.allclose(B.conj().T @ B, np.eye(3), atol=1e-9)
    Qc = B.conj().T @ (ops["Qx"] @ B)
    q = np.linalg.norm(B.conj().T @ QG)
    Z = Qc / q
    target = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    assert np.allclose(Z, target, atol=1e-8)
    A = np.array([[1, 0], [0, 1], [0, 0]], dtype=complex)     # in the basis (G, e, m)
    assert np.allclose(Z @ Z, A @ A.conj().T, atol=1e-8)
    Hc = B.conj().T @ (HB @ B)
    HPsi = A.conj().T @ Hc @ A
    assert np.allclose(HPsi, np.diag(np.diag(HPsi)), atol=1e-9) and HPsi[0, 0].real < HPsi[1, 1].real
    assert abs(Hc[2, 1]) > 1e-6                               # the prepared plane is not invariant
    # the low transition element stays close to its valence value 1
    assert abs(np.vdot(e, QG).real - 1) < 0.02


def test_low_carrier_dynamics():
    """Starting from e: defect weight sin^2(2 Omega_B t), and the ground correlation
    exp(-i (E_2 - E_0) t) cos(2 Omega_B t)."""
    ops, H0, HB, Om, w, G, e, m, QG = _carrier()
    B = np.column_stack([G, e, m])
    Hc = B.conj().T @ (HB @ B)
    evals, V = np.linalg.eigh(Hc)
    E0 = w[0]
    E2 = (Hc[1, 1]).real
    for t in (5.0, 80.0, pi / (4 * Om)):
        psi = V @ (np.exp(-1j * evals * t) * (V.conj().T @ np.array([0, 1, 0], dtype=complex)))
        assert np.isclose(abs(psi[2]) ** 2, np.sin(2 * Om * t) ** 2, atol=1e-8)
        assert np.isclose(np.exp(1j * E0 * t) * psi[1], np.exp(-1j * (E2 - E0) * t) * np.cos(2 * Om * t), atol=1e-8)
    full_departure = V @ (np.exp(-1j * evals * pi / (4 * Om)) * (V.conj().T @ np.array([0, 1, 0], dtype=complex)))
    assert np.isclose(abs(full_departure[2]) ** 2, 1.0, atol=1e-8)


def test_unfiltered_quadrupole_response():
    """For e_Q = Q_x G / ||Q_x G||: C_B(t) = C_0(t) cos(2 Omega_B t).  Only part of e_Q lies
    in the low carrier, and |C_0| can fall below one at sampled times."""
    ops, H0, HB, Om, w, G, e, m, QG = _carrier()
    eQ = QG / np.linalg.norm(QG)
    w0, V0 = np.linalg.eigh(H0.toarray())
    wB, VB = np.linalg.eigh(HB.toarray())
    p0, pB = np.abs(V0.conj().T @ eQ) ** 2, np.abs(VB.conj().T @ eQ) ** 2
    E0 = w[0]
    returns = []
    for t in (3.0, 40.0, 400.0, 2500.0):
        c0 = np.sum(p0 * np.exp(-1j * t * (w0 - E0)))
        cB = np.sum(pB * np.exp(-1j * t * (wB - E0)))
        assert np.isclose(cB, c0 * np.cos(2 * Om * t), atol=1e-9)
        returns.append(abs(c0))
    B = np.column_stack([G, e, m])
    fraction = np.linalg.norm(B.conj().T @ eQ) ** 2
    assert 0.3 < fraction < 0.5
    assert min(returns) < 0.5


def test_basis_convergence_full():
    """Slow: repeat the spectral checks with one more oscillator shell."""
    if not os.environ.get("CHECKS_FULL"):
        skip("run with --full")
    _spectral_structure(4)
    d_a = first_order(1.0, XI)["d_a"]
    gaps = []
    for Nmax in (3, 4):
        _, _, w, _ = _zero_bias(Nmax, LAM)
        gaps.append((w[1] - w[0]) / (LAM * d_a))
    assert abs(gaps[0] - gaps[1]) < 2e-3
