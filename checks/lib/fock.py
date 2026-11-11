"""Fermionic mode operators as explicit matrices (Jordan-Wigner construction)."""
import numpy as np

_I2 = np.eye(2)
_Z = np.diag([1.0, -1.0])
_LOWER = np.array([[0.0, 1.0], [0.0, 0.0]])  # maps occupied -> empty


def modes(n):
    """Annihilation operators a_0, ..., a_{n-1} on the 2**n-dimensional Fock space.

    The basis is the tensor product of (empty, occupied) for each mode, so the
    vector with a single 1 in position 0 is the empty state.
    """
    out = []
    for k in range(n):
        factors = [_Z] * k + [_LOWER] + [_I2] * (n - 1 - k)
        M = factors[0]
        for f in factors[1:]:
            M = np.kron(M, f)
        out.append(M)
    return out


def dag(M):
    return M.conj().T


def vacuum(n):
    v = np.zeros(2**n)
    v[0] = 1.0
    return v


def number(a):
    return dag(a) @ a


def satisfies_car(ops, atol=1e-12):
    """True when the operators obey the canonical anticommutation relations."""
    dim = ops[0].shape[0]
    for i, a in enumerate(ops):
        for k, b in enumerate(ops):
            if not np.allclose(a @ b + b @ a, 0, atol=atol):
                return False
            target = np.eye(dim) if i == k else np.zeros((dim, dim))
            if not np.allclose(a @ dag(b) + dag(b) @ a, target, atol=atol):
                return False
    return True
