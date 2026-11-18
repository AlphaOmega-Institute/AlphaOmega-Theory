"""Relaxation toward a ground vector y: the complete semigroup statement, the distance to
the ground, the absence of a unitary substitute, successive moments, and the state
returned through a pairing.

    D_t(X) = L X L + (1 - s^2) Tr(Q X) P,   L = P + s Q,   s = exp(-gamma t / 2),
    P the projection onto y,  Q = I - P.
"""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm, null_space, polar

from lib.linalg import hermitian_power

LABELS = {
    "LG": ["thm:exact-rigpa-density", "def:ground-nonaltering-channel", "thm:no-unitary-deconditioning",
           "def:ground-attractive-family", "thm:canonical-deconditioning-semigroup",
           "eq:canonical-deconditioning-channel", "eq:deconditioning-exact-contraction",
           "eq:deconditioning-trace-distance-bound", "eq:deconditioning-generator",
           "eq:ground-deviation-lyapunov", "rem:reduced-irreversibility",
           "def:momentary-ground-covariance", "prop:exact-realization-absorbing",
           "thm:momentary-ground-stabilization", "def:realization-event-stabilization",
           "def:conditional-projective-actualization", "prop:return-operation-trace",
           "thm:deconditioning-contracts-excess", "eq:restricted-return-after-deconditioning",
           "eq:global-return-after-deconditioning"],
    "UM": ["eq:relaxation-deviation", "prop:relaxation-pure-distance", "eq:relaxation-pure-distance",
           "prop:relaxation-semigroup", "eq:relaxation-generator", "eq:relaxation-decay", "eq:relaxation-distance",
           "eq:relaxation-blocks", "eq:relaxation-jumps", "eq:relaxation-sharp-distance",
           "prop:relaxation-carried", "eq:relaxation-carried", "eq:relaxation-unitary", "eq:relaxation-dilation",
           "eq:relaxation-dilation-trace", "eq:relaxation-ground-carried", "eq:relaxation-contraction",
           "prop:relaxation-moments", "eq:relaxation-product", "eq:relaxation-return",
           "eq:relaxation-paired-weight", "eq:relaxation-return-bound", "eq:relaxation-restricted-return",
           "eq:moment-deconditioning", "eq:moment-excess-decay", "prop:two-moment-deconditioning",
           "eq:restricted-ground-orientation", "eq:actual-ground-partner"],
}


def _unit(rng, n):
    v = rng.normal(size=n) + 1j * rng.normal(size=n)
    return v / np.linalg.norm(v)


def _state(rng, n, rank=None):
    raw = rng.normal(size=(n, rank or n)) + 1j * rng.normal(size=(n, rank or n))
    rho = raw @ raw.conj().T
    return rho / np.trace(rho).real


def _projection(y):
    return np.outer(y, y.conj())


def _trace_norm(X):
    return np.linalg.svd(X, compute_uv=False).sum()


def _relax(X, y, s):
    P = _projection(y)
    Q = np.eye(len(y)) - P
    L = P + s * Q
    return L @ X @ L + (1 - s**2) * np.trace(Q @ X) * P


def _generator(X, y, gamma):
    P = _projection(y)
    Q = np.eye(len(y)) - P
    return gamma * (np.trace(Q @ X) * P - (Q @ X + X @ Q) / 2)


def _deviation(sigma, y):
    return 1 - (y.conj() @ sigma @ y).real


def _as_matrix(linear_map, n):
    """Matrix of a linear map on n-by-n matrices, in the basis of matrix units."""
    units = np.eye(n * n).reshape(n * n, n, n).astype(complex)
    return np.array([linear_map(unit).reshape(-1) for unit in units]).T


def test_distance_from_a_pure_state():
    rng = np.random.default_rng(31)
    for _ in range(400):
        n = int(rng.integers(2, 7))
        y = _unit(rng, n)
        sigma = _state(rng, n, rank=int(rng.integers(1, n + 1)))
        deviation = _deviation(sigma, y)
        distance = _trace_norm(sigma - _projection(y))
        assert 2 * deviation - 1e-12 <= distance <= 2 * np.sqrt(deviation) + 1e-12
        # equality above for pure states
        phi = _unit(rng, n)
        assert_allclose(_trace_norm(_projection(phi) - _projection(y)),
                        2 * np.sqrt(1 - abs(np.vdot(y, phi))**2), atol=1e-12)
        # equality below for states that commute with the projection
        P = _projection(y)
        Q = np.eye(n) - P
        pinched = P @ sigma @ P + Q @ sigma @ Q
        assert_allclose(_trace_norm(pinched - P), 2 * _deviation(pinched, y), atol=1e-12)
    # the upper bound is what makes zero weight outside the ground force the ground state
    y = np.array([1, 0, 0], dtype=complex)
    for small in (1e-2, 1e-4, 1e-8):
        mixed = np.diag([1 - small, small / 2, small / 2]).astype(complex)
        assert _trace_norm(mixed - _projection(y)) <= 2 * np.sqrt(small)


def test_semigroup_continuity_and_generator():
    rng = np.random.default_rng(32)
    for _ in range(60):
        n = int(rng.integers(2, 6))
        y = _unit(rng, n)
        gamma, t, r = rng.uniform(0.2, 2), rng.uniform(0.05, 3), rng.uniform(0, 3)
        s, s_r = np.exp(-gamma * t / 2), np.exp(-gamma * r / 2)
        X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))       # any matrix
        # composition law, and the same formula as the four blocks of LG's statement
        assert_allclose(_relax(_relax(X, y, s_r), y, s), _relax(X, y, s * s_r), atol=1e-12)
        P = _projection(y)
        Q = np.eye(n) - P
        four_blocks = P @ X @ P + s * (P @ X @ Q + Q @ X @ P) + s**2 * Q @ X @ Q + (1 - s**2) * np.trace(Q @ X) * P
        assert_allclose(_relax(X, y, s), four_blocks, atol=1e-12)
        # trace preserved; continuity in norm on the trace class
        assert_allclose(np.trace(_relax(X, y, s)), np.trace(X), atol=1e-12)
        assert _trace_norm(_relax(X, y, s) - X) <= 4 * (1 - s**2) * _trace_norm(X) + 1e-12
        # bounded generator, and the relaxation as its exponential
        assert _trace_norm(_generator(X, y, gamma)) <= 2 * gamma * _trace_norm(X) + 1e-12
        G = _as_matrix(lambda Z: _generator(Z, y, gamma), n)
        D = _as_matrix(lambda Z: _relax(Z, y, s), n)
        assert_allclose(expm(t * G), D, atol=1e-10)
        # the generator in the form with one jump operator for each direction orthogonal to y
        basis = null_space(y.conj()[None, :])
        jumps = [np.sqrt(gamma) * np.outer(y, e.conj()) for e in basis.T]
        lindblad = sum(J @ X @ J.conj().T - (J.conj().T @ J @ X + X @ J.conj().T @ J) / 2 for J in jumps)
        assert_allclose(lindblad, _generator(X, y, gamma), atol=1e-12)
        # at a positive time the fixed operators are the multiples of the projection
        values, vectors = np.linalg.eig(D)
        fixed = np.flatnonzero(np.abs(values - 1) < 1e-9)
        assert len(fixed) == 1
        fixed_operator = vectors[:, fixed[0]].reshape(n, n)
        assert_allclose(fixed_operator / np.trace(fixed_operator), P, atol=1e-9)
        assert_allclose(_relax(P, y, s), P, atol=1e-13)
        # a density operator is sent to a density operator
        out = _relax(_state(rng, n), y, s)
        assert np.linalg.eigvalsh((out + out.conj().T) / 2).min() > -1e-12


def test_weight_outside_the_ground_decays_exactly():
    rng = np.random.default_rng(33)
    for _ in range(200):
        n = int(rng.integers(2, 7))
        y = _unit(rng, n)
        sigma = _state(rng, n)
        gamma, t = rng.uniform(0.2, 2), rng.uniform(0, 4)
        s = np.exp(-gamma * t / 2)
        relaxed = _relax(sigma, y, s)
        assert_allclose(_deviation(relaxed, y), np.exp(-gamma * t) * _deviation(sigma, y), atol=1e-12)
        # its rate of change is -gamma times itself
        rate = np.trace((np.eye(n) - _projection(y)) @ _generator(relaxed, y, gamma)).real
        assert_allclose(rate, -gamma * _deviation(relaxed, y), atol=1e-12)
        # a state off the ground is off the ground at every finite time
        assert _deviation(relaxed, y) > 0
        # the reset to the ground is the value s = 0, which no finite time attains
        assert_allclose(_relax(sigma, y, 0.0), _projection(y), atol=1e-12)


def test_distance_to_the_ground_under_relaxation():
    rng = np.random.default_rng(34)
    closest = 0.0
    for _ in range(1500):
        n = int(rng.integers(2, 7))
        y = _unit(rng, n)
        sigma = _state(rng, n, rank=int(rng.integers(1, n + 1)))
        gamma, t = rng.uniform(0.2, 2), rng.uniform(0, 6)
        s = np.exp(-gamma * t / 2)
        deviation = _deviation(sigma, y)
        distance = _trace_norm(_relax(sigma, y, s) - _projection(y))
        assert 2 * s**2 * deviation - 1e-12 <= distance <= 2 * s * np.sqrt(deviation) + 1e-12
        closest = max(closest, distance / (2 * s * np.sqrt(deviation)))
        # the sharper bound, an equality for pure inputs; it makes the bound above strict at positive times
        sharper = 2 * s * np.sqrt(deviation) * np.sqrt(1 - (1 - s**2) * deviation)
        assert distance <= sharper + 1e-12
        if t > 0.05:
            assert sharper < 2 * s * np.sqrt(deviation)
        phi = _unit(rng, n)
        weight = _deviation(_projection(phi), y)
        assert_allclose(_trace_norm(_relax(_projection(phi), y, s) - _projection(y)),
                        2 * s * np.sqrt(weight) * np.sqrt(1 - (1 - s**2) * weight), atol=1e-12)
        # without coherence with the ground the distance decays at the full rate
        P = _projection(y)
        Q = np.eye(n) - P
        pinched = P @ sigma @ P + Q @ sigma @ Q
        assert_allclose(_trace_norm(_relax(pinched, y, s) - P), 2 * s**2 * deviation, atol=1e-12)
    assert closest > 0.9
    # a pure input at angle theta from the ground
    y = np.array([1, 0, 0], dtype=complex)
    e = np.array([0, 0.6, 0.8j])
    for s in (1.0, 0.7, 0.2):
        ratios = []
        for theta in (1.2, 0.4, 0.05, 0.002):
            phi = np.cos(theta) * y + np.sin(theta) * e
            distance = _trace_norm(_relax(_projection(phi), y, s) - _projection(y))
            assert_allclose(distance, 2 * s * np.sin(theta) * np.sqrt(1 - (1 - s**2) * np.sin(theta)**2), atol=1e-13)
            ratios.append(distance / (2 * s * np.sin(theta)))
        # the bound 2 s sqrt(deviation) is approached as the input approaches the ground
        assert np.all(np.diff(ratios) > -1e-12) and ratios[-1] > 1 - 1e-5
        if s == 1.0:
            assert_allclose(ratios, 1.0)                       # equality at time zero
    # one fixed pure input already shows the rate gamma/2: distance / s tends to 2 sqrt(d (1 - d))
    for theta in (0.3, 0.9, 1.4):
        phi = np.cos(theta) * y + np.sin(theta) * e
        weight = np.sin(theta)**2
        late = _trace_norm(_relax(_projection(phi), y, 1e-5) - _projection(y)) / 1e-5
        assert_allclose(late, 2 * np.sqrt(weight * (1 - weight)), rtol=1e-6)


def _carrying_channel(rng, y, y_next, env):
    """A random channel with T(P_y) = P_y': an isometry W into (next space) x (environment)
    with W y = y' (x) k. Every channel that carries the ground has this form."""
    n, m = len(y), len(y_next)
    first = np.kron(y_next, _unit(rng, env))
    rest = rng.normal(size=(m * env, n - 1)) + 1j * rng.normal(size=(m * env, n - 1))
    rest = rest - np.outer(first, first.conj() @ rest)
    rest, _ = np.linalg.qr(rest)
    W = np.column_stack([first, rest]) @ np.column_stack([y, null_space(y.conj()[None, :])]).conj().T
    assert_allclose(W.conj().T @ W, np.eye(n), atol=1e-10)
    branches = [W.reshape(m, env, n)[:, a, :] for a in range(env)]             # one Kraus operator for each record
    return (lambda sigma: np.einsum("iaja->ij", (W @ sigma @ W.conj().T).reshape(m, env, m, env))), branches


def test_a_channel_carries_the_ground_exactly_when_no_weight_outside_it_grows():
    rng = np.random.default_rng(38)
    largest_gain, some_strict_loss = -1.0, False
    for _ in range(600):
        n, m = int(rng.integers(2, 5)), int(rng.integers(2, 5))
        env = int(rng.integers(-(-n // m), 5))
        y, y_next = _unit(rng, n), _unit(rng, m)
        channel, branches = _carrying_channel(rng, y, y_next, env)
        assert_allclose(channel(_projection(y)), _projection(y_next), atol=1e-10)
        # every outcome of an instrument with this channel leaves the ground on the next ground
        assert_allclose(sum(K.conj().T @ K for K in branches), np.eye(n), atol=1e-10)
        for K in branches:
            image = K @ y
            assert_allclose(image, np.vdot(y_next, image) * y_next, atol=1e-10)
        for _ in range(4):
            sigma = _state(rng, n, rank=int(rng.integers(1, n + 1)))
            out = channel(sigma)
            assert_allclose(np.trace(out).real, 1.0, atol=1e-10)
            gain = _deviation(out, y_next) - _deviation(sigma, y)
            largest_gain = max(largest_gain, gain)
            some_strict_loss = some_strict_loss or gain < -1e-3
    assert largest_gain < 1e-10 and some_strict_loss
    # a channel that does not carry the ground raises the weight outside it for the ground itself
    y, other = _unit(rng, 4), _unit(rng, 4)
    assert _deviation(_relax(_projection(y), other, 0.5), y) > 0
    # a conditional step can fix the ground and raise the weight outside it
    y = np.array([1, 0, 0], dtype=complex)
    e = np.array([0, 0.6, 0.8j])
    K = _projection(y) / 2 + (np.eye(3) - _projection(y))
    assert np.linalg.eigvalsh(np.eye(3) - K.conj().T @ K).min() > -1e-12        # one outcome of an instrument
    fixed = K @ _projection(y) @ K
    assert_allclose(fixed / np.trace(fixed), _projection(y), atol=1e-12)
    for theta in (0.1, 0.7, 1.3):
        phi = np.cos(theta) * y + np.sin(theta) * e
        out = K @ _projection(phi) @ K
        out = out / np.trace(out)
        weight = np.sin(theta)**2
        assert_allclose(_deviation(out, y), 4 * weight / (1 + 3 * weight), atol=1e-12)
        assert _deviation(out, y) > weight
    for _ in range(50):                                        # the same formula for mixed states
        sigma = _state(rng, 3, rank=int(rng.integers(1, 4)))
        weight = _deviation(sigma, y)
        out = K @ sigma @ K
        assert_allclose(np.trace(out).real, (1 + 3 * weight) / 4, atol=1e-12)
        assert_allclose(_deviation(out / np.trace(out), y), 4 * weight / (1 + 3 * weight), atol=1e-12)
    # an outcome that takes the ground to another ray belongs to an instrument whose channel moves the ground
    X = 0.6 * (_projection(y) + np.outer(e, y.conj()))                         # X y is not a multiple of y
    rest = hermitian_power(np.eye(3) - X.conj().T @ X, 0.5)
    total = X @ _projection(y) @ X.conj().T + rest @ _projection(y) @ rest.conj().T
    assert_allclose(np.trace(total).real, 1.0, atol=1e-9)
    assert _deviation(total, y) > 0.1


def test_no_unitary_substitute_and_a_one_step_dilation():
    rng = np.random.default_rng(35)
    for _ in range(40):
        n = int(rng.integers(2, 6))
        y = _unit(rng, n)
        basis = null_space(y.conj()[None, :])                 # orthonormal basis of the complement of y
        sigma = _state(rng, n)
        # a unitary that keeps the ground ray leaves the weight outside it unchanged
        inner, _ = np.linalg.qr(rng.normal(size=(n - 1, n - 1)) + 1j * rng.normal(size=(n - 1, n - 1)))
        U = np.exp(1j * rng.uniform(0, 6)) * _projection(y) + basis @ inner @ basis.conj().T
        assert_allclose(U @ U.conj().T, np.eye(n), atol=1e-12)
        assert_allclose(U @ _projection(y) @ U.conj().T, _projection(y), atol=1e-12)
        assert_allclose(_deviation(U @ sigma @ U.conj().T, y), _deviation(sigma, y), atol=1e-12)
        # environment C + (complement of y), reference vector (1, 0); the isometry of the text
        s = rng.uniform(0.05, 0.95)
        P = _projection(y)
        L = P + s * (np.eye(n) - P)
        reference = np.eye(n)[:, 0]
        iota = np.kron(L, reference[:, None])
        for j in range(n - 1):
            record = np.eye(n)[:, j + 1]
            iota = iota + np.sqrt(1 - s**2) * np.kron(np.outer(y, basis[:, j].conj()), record[:, None])
        assert_allclose(iota.conj().T @ iota, np.eye(n), atol=1e-12)
        joint = (iota @ sigma @ iota.conj().T).reshape(n, n, n, n)
        assert_allclose(np.einsum("iaja->ij", joint), _relax(sigma, y, s), atol=1e-12)
        # extend it to a unitary of system and environment
        start = np.kron(np.eye(n), reference[:, None])          # h -> h (x) reference
        total = iota @ start.conj().T + null_space(iota.conj().T) @ null_space(start.conj().T).conj().T
        assert_allclose(total.conj().T @ total, np.eye(n * n), atol=1e-10)
        prepared = np.kron(sigma, np.outer(reference, reference))
        evolved = (total @ prepared @ total.conj().T).reshape(n, n, n, n)
        assert_allclose(np.einsum("iaja->ij", evolved), _relax(sigma, y, s), atol=1e-10)
        # the record: the environment is found outside its reference with the probability lost by the complement
        assert_allclose(1 - np.einsum("iaib->ab", evolved)[0, 0].real, (1 - s**2) * _deviation(sigma, y), atol=1e-10)


def test_relaxation_in_successive_moments():
    rng = np.random.default_rng(36)
    dims = [2, 3, 3, 4, 5, 5, 6]
    grounds = [_unit(rng, d) for d in dims]
    steps = []
    for n in range(len(dims) - 1):
        # an isometry into the next space that carries this ground to the next ground
        here = np.column_stack([grounds[n], null_space(grounds[n].conj()[None, :])])
        there = np.column_stack([grounds[n + 1], null_space(grounds[n + 1].conj()[None, :])])
        steps.append(there[:, :dims[n]] @ here.conj().T)
        assert_allclose(steps[-1] @ grounds[n], grounds[n + 1], atol=1e-12)
    factors = rng.uniform(0.2, 1.0, size=len(steps))
    factors[2] = 1.0                                          # a step without progress
    sigma = _state(rng, dims[0])
    ground_state = _projection(grounds[0])
    deviations = [_deviation(sigma, grounds[0])]
    for n, (carry, q) in enumerate(zip(steps, factors)):
        sigma = carry @ _relax(sigma, grounds[n], np.sqrt(q)) @ carry.conj().T
        ground_state = carry @ _relax(ground_state, grounds[n], np.sqrt(q)) @ carry.conj().T
        assert_allclose(np.trace(sigma).real, 1.0, atol=1e-12)
        deviations.append(_deviation(sigma, grounds[n + 1]))
        assert_allclose(deviations[-1], q * deviations[-2], atol=1e-12)
        assert_allclose(ground_state, _projection(grounds[n + 1]), atol=1e-12)        # once at the ground, always
        assert _trace_norm(sigma - _projection(grounds[n + 1])) <= 2 * np.sqrt(deviations[-1]) + 1e-12
    assert_allclose(deviations[-1], np.prod(factors) * deviations[0], atol=1e-12)
    # steps that carry the ground never let that weight grow, whatever the state
    assert np.all(np.diff(deviations) <= 1e-12)
    # two steps that carry the ground, each leaving one excited state alone, together reach the ground
    g, a, b = np.eye(3)
    first = [_projection(g) + _projection(b), np.outer(g, a)]
    second = [_projection(g) + _projection(a), np.outer(g, b)]
    apply = lambda kraus, state: sum(K @ state @ K.conj().T for K in kraus)
    for kraus, untouched in ((first, b), (second, a)):
        assert_allclose(sum(K.conj().T @ K for K in kraus), np.eye(3), atol=1e-14)
        assert_allclose(apply(kraus, _projection(g)), _projection(g), atol=1e-14)
        assert_allclose(apply(kraus, _projection(untouched)), _projection(untouched), atol=1e-14)
    for _ in range(20):
        assert_allclose(apply(second, apply(first, _state(rng, 3))), _projection(g), atol=1e-12)
    # exponents 1/n have a divergent sum and exponents 1/n^2 a convergent one
    n = np.arange(1, 200001)
    assert np.exp(-np.sum(1 / n)) < 1e-5
    assert_allclose(np.exp(-np.sum(1 / n**2)), np.exp(-np.pi**2 / 6), rtol=1e-5)
    assert np.exp(-np.pi**2 / 6) > 0.19


def test_return_after_relaxation():
    rng = np.random.default_rng(37)
    for _ in range(200):
        m = int(rng.integers(1, 5))
        n = m + int(rng.integers(1, 4))
        V, _ = np.linalg.qr(rng.normal(size=(n, m)) + 1j * rng.normal(size=(n, m)))     # isometry onto D
        P_D = V @ V.conj().T
        p = _unit(rng, m)
        y = V @ p                                             # target in D with partner p
        sigma = _state(rng, n, rank=int(rng.integers(1, n + 1)))
        # the return operation and its trace; for a pure state it is the ray of V* phi
        assert_allclose(np.trace(V.conj().T @ sigma @ V).real, np.trace(P_D @ sigma).real, atol=1e-12)
        phi = _unit(rng, n)
        back = V.conj().T @ phi
        returned = V.conj().T @ _projection(phi) @ V / np.trace(P_D @ _projection(phi)).real
        assert_allclose(returned, _projection(back / np.linalg.norm(back)), atol=1e-12)
        # an instrument that prepares a fixed state after the same outcome: same probability, another state
        prepared = _state(rng, m)
        assert_allclose(np.trace(np.trace(P_D @ sigma) * prepared), np.trace(V.conj().T @ sigma @ V), atol=1e-12)
        if m > 1:
            returned_state = V.conj().T @ sigma @ V / np.trace(P_D @ sigma).real
            assert not np.allclose(prepared, returned_state, atol=1e-6)
            values, vectors = np.linalg.eigh(prepared - returned_state)
            telling = _projection(vectors[:, -1])                  # an effect with different expectations
            assert abs(np.trace(telling @ (prepared - returned_state))) > 1e-6
            assert_allclose(np.trace(np.eye(m) @ (prepared - returned_state)), 0, atol=1e-12)    # the trivial effect tells nothing
        gamma, t = rng.uniform(0.2, 2), rng.uniform(0.01, 5)
        s = np.exp(-gamma * t / 2)
        relaxed = _relax(sigma, y, s)
        outside = 1 - np.trace(P_D @ sigma).real
        weight = np.trace(P_D @ relaxed).real
        assert_allclose(weight, 1 - s**2 * outside, atol=1e-12)
        assert outside <= _deviation(sigma, y) + 1e-12
        # global return: normalized compression, close to the partner of the target
        conditional = V.conj().T @ relaxed @ V / weight
        bound = 2 * s * np.sqrt(_deviation(sigma, y)) + s**2 * outside
        assert _trace_norm(conditional - _projection(p)) <= bound + 1e-12
        # restricted return through u0 = |y><x|: exactly the reference state, with this probability
        x = _unit(rng, m)
        u0 = np.outer(y, x.conj())
        restricted = u0.conj().T @ relaxed @ u0
        probability = 1 - s**2 * _deviation(sigma, y)
        assert probability > 0
        assert_allclose(restricted, probability * _projection(x), atol=1e-12)
    # the unequal-cell path: fine site 2 selected, then relaxed halfway toward the fine ground
    copying = np.repeat(np.eye(2), [1, 4], axis=0) / 2
    V, _ = polar(copying)
    y = np.ones(5) / np.sqrt(5)
    x = np.ones(2) / np.sqrt(2)
    p = V.T @ y
    assert_allclose(p, np.array([1, 2]) / np.sqrt(5), atol=1e-12)
    start = _projection(np.eye(5)[1])
    assert_allclose(1 - np.trace(V @ V.T @ start), 3 / 4)
    assert_allclose(_deviation(start, y), 4 / 5)
    half = _relax(start, y, np.sqrt(0.5))
    assert_allclose(np.trace(V @ V.T @ half).real, 5 / 8)
    assert_allclose((y @ half @ y).real, 3 / 5)
    # the returned state passes the coarse ground test with a probability that tends to 9/10
    for s in (np.sqrt(0.5), 0.3, 0.05, 0.001):
        relaxed = _relax(start, y, s)
        conditional = V.T @ relaxed @ V / np.trace(V @ V.T @ relaxed).real
        exact = ((6 - s)**2 / 200 + 18 * (1 - s**2) / 25) / (1 - 3 * s**2 / 4)
        assert_allclose((x @ conditional @ x).real, exact, atol=1e-12)
    assert_allclose(exact, 9 / 10, atol=1e-4)
    # the restricted return gives the coarse ground itself, here with probability 3/5
    u0 = np.outer(y, x)
    assert_allclose(u0.T @ half @ u0, 3 / 5 * _projection(x), atol=1e-12)
