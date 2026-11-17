"""Effective Hamiltonians, shell error bounds, and rotated orbital detectors (hbar = 1)."""
import numpy as np
from scipy.linalg import expm

from orbital.model import build, first_order, lowest
from pair.model import PairModel

LABELS = {
    "UM": ["eq:prepared-schur", "prop:shell-accuracy", "eq:rotated-orbital-readout"],
}


def test_pair_schur_resolvent_and_silent_levels():
    """Moving the silent energies changes the ambient spectrum while preserving the port resolvent."""
    silent_block = np.array([[3, 1j, 1], [-1j, -2, 0.5], [1, 0.5, 4]])
    for theta in (0, np.pi / 18, np.pi / 8):
        pair = PairModel(theta)
        shell = pair.F2
        port = shell.conj().T @ pair.A
        other = shell.conj().T @ np.column_stack((pair.dvec, pair.W))
        h = shell.conj().T @ pair.H @ shell
        silent = shell.conj().T @ pair.W
        moved = h + silent @ silent_block @ silent.conj().T
        assert not np.allclose(np.linalg.eigvalsh(h), np.linalg.eigvalsh(moved))
        coupling = other.conj().T @ h @ port
        direct = port.conj().T @ h @ port
        complement = other.conj().T @ h @ other
        for z in (0.19j, 0.7 + 0.15j):
            sigma = coupling.conj().T @ np.linalg.solve(complement - z * np.eye(4), coupling)
            reduced = np.linalg.inv(direct - z * np.eye(2) - sigma)
            for ambient in (h, moved):
                measured = port.conj().T @ np.linalg.solve(ambient - z * np.eye(6), port)
                np.testing.assert_allclose(measured, reduced, atol=2e-12)
            excitation = np.vdot(pair.ePsi, reduced @ pair.ePsi)
            expected = pair.u**2 / (pair.j - z) + pair.v**2 / (-z)
            np.testing.assert_allclose(excitation, expected, atol=2e-12)


def test_sharp_ordered_shell_energies():
    """Place an independent shell level below the extremizer or between its two variational energies."""
    gap = 1.7
    for width in (0.12, 0.8, 1.4):
        sine = width / (2 * gap)
        cosine = np.sqrt(1 - sine**2)
        active = -width / 2 * np.array([[1 - sine, cosine], [cosine, 1 + sine]])
        bound = width**2 / (4 * gap)
        for spectator in (-width, -width / 2 + bound / 2):
            h = np.diag([0.0, spectator, gap])
            h[np.ix_([0, 2], [0, 2])] += active
            ritz = np.linalg.eigvalsh(h[:2, :2])
            exact = np.linalg.eigvalsh(h)[:2]
            expected = [0, bound] if spectator == -width else [bound / 2, bound / 2]
            np.testing.assert_allclose(ritz - exact, expected, atol=2e-14)


def test_sharp_shell_leakage_and_accumulating_amplitude_error():
    """A two-level reflection saturates leakage; at complete return a phase error can persist."""
    gap = 1.7
    for ratio in (0.1, 0.5, 0.9):
        width = ratio * gap
        off_diagonal = np.sqrt(1 - ratio**2)
        interaction = (width / 2) * np.array([[ratio, off_diagonal], [off_diagonal, -ratio]])
        h = np.diag([0.0, gap]) + interaction
        peak_time = np.pi / (gap * off_diagonal)
        for multiple in (0.3, 1, 2):
            t = multiple * peak_time
            amplitudes = expm(-1j * h * t)[:, 0]
            leakage = abs(amplitudes[1])**2
            drift = abs(amplitudes[0] - np.exp(-1j * h[0, 0] * t))
            np.testing.assert_allclose(leakage, ratio**2 * np.sin(multiple * np.pi / 2)**2,
                                       atol=2e-14)
            assert drift <= width**2 * t / (2 * gap) + 2e-14
            if multiple == 2:
                assert drift > 0.01


def test_rotated_preparation_with_fixed_and_rotating_orbital_detectors():
    """Rotate the exact filtered quadrupole excitation of a truncated four-fermion parent."""
    xi, lam = 1.0, 0.2
    omega_b = lam * first_order(1.0, xi)["d_a"] / 4
    ops = build(2, xi)
    h = ops["Hho"] - lam * ops["W"] - omega_b * ops["L"]
    energies, low = lowest(h, 3)
    angular = low.conj().T @ (ops["L"] @ low)
    np.testing.assert_allclose(np.diag(angular), [0, 2, -2], atol=1e-10)
    excitation = low.conj().T @ (ops["Qx"] @ low[:, 0])
    excitation /= np.linalg.norm(excitation)
    ground = np.array([1., 0., 0.])
    for angle in np.deg2rad([0, 10, 22.5, 40]):
        turned = expm(-1j * angle * angular) @ excitation
        np.testing.assert_allclose(abs(turned)**2, [0, 0.5, 0.5], atol=2e-10)
        for phase in (0, 0.43, np.pi / 2):
            t = phase / (2 * omega_b)
            evolved = np.exp(-1j * (energies - energies[0]) * t) * turned
            for detector, offset in ((turned, 0), (excitation, 2 * angle)):
                plane = np.column_stack((ground, detector))
                missing = 1 - np.linalg.norm(plane.conj().T @ evolved)**2
                np.testing.assert_allclose(missing, np.sin(phase - offset)**2, atol=2e-10)
