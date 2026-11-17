"""Four fermions in the two-dimensional trap, by an independent determinant calculation.

Units hbar = M = omega = 1, so the oscillator length is 1 and xi = a^2.  This file shares
no code with orbital/model.py.  Single-particle states are Cartesian oscillator functions
(n_x, n_y) with spin; the Gaussian pair attraction factorises into an x part and a y part,
each evaluated exactly by Gauss-Hermite quadrature in the centre and relative coordinates.
The interaction is applied to the six determinants of the s^2 p^2 shell without forming a
Hamiltonian matrix.

It checks, for UM Section 8:
  * the interaction has no matrix element from the shell to a configuration with an odd
    number of oscillator quanta (evaluated with no parity screening);
  * the first-order shell energies, both the absolute ground value and the two gaps;
  * the second-order coefficients in the manuscript's table;
  * which single-particle cutoff completes each level of the second-order sum.
The standard run uses a cutoff of six quanta per orbital; `--full` repeats it at ten.
"""
import itertools
import os

import numpy as np

from lib import skip

LABELS = {
    "UM": ["eq:orbital-shell-ritz", "eq:orbital-second-order", "eq:orbital-shell-accuracy"],
    "LG": ["eq:orbital-projected-energies", "eq:orbital-three-orbitals"],
}


def _pair_tensor(nmax, xi, nodes=64):
    """T[i, j, k, l] = integral of h_i(x1) h_j(x2) h_k(x1) h_l(x2) exp(-(x1 - x2)^2 / (2 xi))."""
    t, w = np.polynomial.hermite.hermgauss(nodes)
    beta = 1 + 1 / xi
    centre, relative = np.meshgrid(t, t / np.sqrt(beta), indexing="ij")
    weight = (w[:, None] * w[None, :] / np.sqrt(beta)).ravel()
    x1, x2 = ((centre + relative) / np.sqrt(2)).ravel(), ((centre - relative) / np.sqrt(2)).ravel()

    def polynomials(x):
        p = np.zeros((nmax + 1, x.size))
        p[0] = np.pi**-0.25
        if nmax > 0:
            p[1] = np.sqrt(2.0) * x * p[0]
        for n in range(1, nmax):
            p[n + 1] = np.sqrt(2.0 / (n + 1)) * x * p[n] - np.sqrt(n / (n + 1)) * p[n - 1]
        return p

    p1, p2 = polynomials(x1), polynomials(x2)
    return np.einsum("im,jm,km,lm,m->ijkl", p1, p2, p1, p2, weight, optimize=True)


class Shell:
    """V applied to the six shell determinants, for lambda = 1."""

    def __init__(self, cutoff, xi):
        self.xi = xi
        self.tensor = _pair_tensor(cutoff, xi)
        self.orbitals = [(nx, n - nx) for n in range(cutoff + 1) for nx in range(n + 1)]
        self.quanta = [nx + ny for nx, ny in self.orbitals]
        self.count = 2 * len(self.orbitals)                    # spin orbital p = 2 * orbital + spin
        ix, iy = self.orbitals.index((1, 0)), self.orbitals.index((0, 1))
        valence = [2 * ix, 2 * ix + 1, 2 * iy, 2 * iy + 1]
        self.shell = [tuple(sorted((0, 1) + pair)) for pair in itertools.combinations(valence, 2)]
        self.images = [self._apply(det) for det in self.shell]

    def element(self, p, q, r, s):
        """<p q| v |r s> with particle 1 in p, r and particle 2 in q, s."""
        if (p & 1) != (r & 1) or (q & 1) != (s & 1):
            return 0.0
        a, b, c, d = (self.orbitals[i >> 1] for i in (p, q, r, s))
        return -self.tensor[a[0], b[0], c[0], d[0]] * self.tensor[a[1], b[1], c[1], d[1]] / (2 * np.pi * self.xi)

    def _apply(self, det, screen_parity=True, limit=None):
        out = {}
        last = self.count if limit is None else min(self.count, limit)
        for i in range(4):
            for j in range(i + 1, 4):
                r, s = det[i], det[j]
                rest = [o for k, o in enumerate(det) if k not in (i, j)]
                parity = (self.quanta[r >> 1] + self.quanta[s >> 1]) & 1
                for p in range(last):
                    if p in rest:
                        continue
                    for q in range(p + 1, last):
                        if q in rest:
                            continue
                        if screen_parity and ((self.quanta[p >> 1] + self.quanta[q >> 1]) & 1) != parity:
                            continue
                        amplitude = self.element(p, q, r, s) - self.element(p, q, s, r)
                        if amplitude == 0.0:
                            continue
                        new = sorted(rest + [p, q])
                        sign = (-1) ** (i + j - 1) * (-1) ** (new.index(p) + new.index(q) - 1)
                        key = tuple(new)
                        out[key] = out.get(key, 0.0) + sign * amplitude
        return out

    def total_quanta(self, det):
        return sum(self.quanta[p >> 1] for p in det)

    def first_order(self):
        return np.array([[self.images[j].get(self.shell[i], 0.0) for j in range(6)] for i in range(6)])

    def outside(self):
        """Rows: determinants outside the shell reached by V; returns their levels, highest orbitals and amplitudes."""
        keys = sorted(set().union(*[set(image) for image in self.images]) - set(self.shell))
        amplitudes = np.array([[image.get(key, 0.0) for image in self.images] for key in keys])
        levels = np.array([self.total_quanta(key) - 2 for key in keys])
        highest = np.array([max(self.quanta[p >> 1] for p in key) for key in keys])
        return levels, highest, amplitudes


def _second_order(shell, cutoff=None):
    levels, highest, amplitudes = shell.outside()
    keep = levels > 0 if cutoff is None else (levels > 0) & (highest <= cutoff)
    correction = -(amplitudes[keep].T / levels[keep]) @ amplitudes[keep]
    _, basis = np.linalg.eigh(shell.first_order())
    rotated = basis.T @ correction @ basis
    assert np.abs(rotated - np.diag(np.diag(rotated))).max() < 1e-12
    return np.diag(rotated)


def test_no_coupling_to_odd_configurations():
    for xi in (2.0, 0.5):
        shell = Shell(3, xi)                                    # orbitals up to three quanta hold every level-one target
        strongest_even = 0.0
        for det in shell.shell:
            image = shell._apply(det, screen_parity=False)
            for key, amplitude in image.items():
                if (shell.total_quanta(key) - 2) % 2:
                    assert abs(amplitude) < 1e-14
                elif key not in shell.shell:
                    strongest_even = max(strongest_even, abs(amplitude))
        assert strongest_even > 1e-3


def test_first_order_shell_energies():
    for xi in (0.5, 1.0, 2.0):
        values = np.linalg.eigvalsh(Shell(2, xi).first_order())
        v = -1 / (8 * np.pi)
        ground = 4 * v * (5 * xi + 2) / (1 + xi) ** 2 + 4 * v * (xi**2 + xi + 1) / (1 + xi) ** 3
        d_a = 1 / (4 * np.pi * (1 + xi) ** 3)
        np.testing.assert_allclose(values, ground + d_a * np.array([0, 1, 1, 2, 2, 2]), atol=1e-12)


def test_second_order_coefficients_and_level_completeness():
    shell = Shell(6, 2.0)
    np.testing.assert_allclose(_second_order(shell), [-0.00162813, -0.00170517, -0.00170517,
                                                      -0.00177006, -0.00177006, -0.00177006], atol=6e-9)
    levels, highest, amplitudes = shell.outside()
    assert not np.any(levels % 2)
    # A level k quanta above the shell needs orbitals up to k + 1 and no higher.
    for step, strength in ((2, 0.01972668), (4, 0.00179189)):
        level = levels == step
        assert highest[level].max() == step + 1
        assert np.isclose(np.linalg.norm(amplitudes[level]) ** 2, strength, atol=5e-9)
    assert highest[levels == 6].max() == 6                      # truncated: the full level needs cutoff seven
    assert np.isclose(np.linalg.norm(amplitudes[levels == 6]) ** 2, 0.00020630, atol=5e-9)
    # singlet gap: the two singlet levels move down by nearly the same amount
    coefficients = _second_order(shell)
    assert abs((coefficients[1] - coefficients[0]) + 0.00007704) < 2e-8


def test_second_order_at_cutoff_ten_full():
    if not os.environ.get("CHECKS_FULL"):
        skip("run with --full")
    shell = Shell(10, 2.0)
    np.testing.assert_allclose(_second_order(shell), [-0.00162846, -0.00170561, -0.00170561,
                                                      -0.00177062, -0.00177062, -0.00177062], atol=6e-9)
    np.testing.assert_allclose(_second_order(shell, cutoff=6)[[0, 1, 3]], [-0.00162813, -0.00170517, -0.00177006], atol=6e-9)
    levels, highest, amplitudes = shell.outside()
    level = levels == 6
    assert highest[level].max() == 7
    assert np.isclose(np.linalg.norm(amplitudes[level]) ** 2, 0.00021752, atol=5e-9)
    coefficients = _second_order(shell)
    gap_change = coefficients[1] - coefficients[0]
    assert abs(gap_change + 0.00007715) < 2e-8
    first_gap = 1 / (4 * np.pi * 27)
    assert abs(abs(gap_change) * (8 * np.pi / 243) / first_gap - 0.0027) < 1e-4
