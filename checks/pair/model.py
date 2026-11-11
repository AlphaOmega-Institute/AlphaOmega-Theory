"""The four-mode pair-transfer model of LG, section "A finite pair-transfer realization".

Modes c1, c2, d1, d2; pair annihilators P_c = c2 c1, P_d = d2 d1; Hamiltonian
H = -j (P_c^* P_d + P_d^* P_c).  Everything is an explicit matrix on the
16-dimensional Fock space, with hbar = 1.
"""
import numpy as np

from lib.fock import modes, dag, vacuum, number
from lib.linalg import eigenspace


class PairModel:
    def __init__(self, theta=0.3, j=1.0):
        self.theta, self.j = theta, j
        c1, c2, d1, d2 = modes(4)
        self.c, self.d = [c1, c2], [d1, d2]
        self.vac = vacuum(4)
        self.N = sum(number(a) for a in (c1, c2, d1, d2))
        self.Pc, self.Pd = c2 @ c1, d2 @ d1
        self.H = self.exchange_hamiltonian(0.0)
        # two-particle sector F_2 (six-dimensional): columns of F2 are an orthonormal basis
        self.F2 = eigenspace(self.N, 2.0)
        self.P2 = self.F2 @ self.F2.conj().T
        # pair states and the vectors g, e, m
        self.C = dag(self.Pc) @ self.vac
        self.D = dag(self.Pd) @ self.vac
        self.g = (self.C + self.D) / np.sqrt(2)
        self.e = (self.C - self.D) / np.sqrt(2)
        self.m = (dag(c1) @ dag(d2) + dag(d1) @ dag(c2)) @ self.vac / np.sqrt(2)
        self.w0 = (dag(c1) @ dag(d2) - dag(d1) @ dag(c2)) @ self.vac / np.sqrt(2)
        self.W = np.column_stack([dag(c1) @ dag(d1) @ self.vac, dag(c2) @ dag(d2) @ self.vac, self.w0])
        self.Phi = np.column_stack([self.g, self.e, self.m])       # output carrier
        # one-body orbital currents acting identically on both internal labels
        self.T = dag(c1) @ d1 + dag(c2) @ d2
        self.X = (self.T + dag(self.T)) / 2
        self.Y = (self.T - dag(self.T)) / 2j
        self.Z = (number(c1) + number(c2) - number(d1) - number(d2)) / 2
        # rotated orbital modes and the preparation channel at angle theta
        self.u, self.v = np.cos(2 * theta), np.sin(2 * theta)
        self.fdag, self.hdag = self.rotated_creators(theta)
        self.A = self.preparation(theta)                           # 16 x 2, basis (C_Psi, D_Psi)
        self.Ztheta = (sum(x @ dag(x) for x in self.fdag) - sum(x @ dag(x) for x in self.hdag)) / 2
        self.gPsi = np.array([1.0, 1.0]) / np.sqrt(2)
        self.ePsi = np.array([1.0, -1.0]) / np.sqrt(2)
        self.b = self.u * self.e + self.v * self.m
        self.dvec = -self.v * self.e + self.u * self.m

    def exchange_hamiltonian(self, s_c, s_d=None):
        """H = -s_c P_c^*P_c - s_d P_d^*P_d - j (P_c^*P_d + P_d^*P_c)."""
        s_d = s_c if s_d is None else s_d
        Pc, Pd = self.Pc, self.Pd
        return -(s_c * dag(Pc) @ Pc + s_d * dag(Pd) @ Pd + self.j * (dag(Pc) @ Pd + dag(Pd) @ Pc))

    def rotated_creators(self, theta):
        co, si = np.cos(theta), np.sin(theta)
        f = [co * dag(c) + si * dag(d) for c, d in zip(self.c, self.d)]
        h = [-si * dag(c) + co * dag(d) for c, d in zip(self.c, self.d)]
        return f, h

    def preparation(self, theta):
        """A|C_Psi> = f1^* f2^* |vac>,  A|D_Psi> = h1^* h2^* |vac>."""
        f, h = self.rotated_creators(theta)
        return np.column_stack([f[0] @ f[1] @ self.vac, h[0] @ h[1] @ self.vac])

    def on_carrier(self, op):
        """Compression of an operator to the output carrier span{g, e, m}."""
        return self.Phi.conj().T @ op @ self.Phi
