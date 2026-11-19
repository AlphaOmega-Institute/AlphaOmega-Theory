"""Checks for the sharp frame on differentiated content: the frame relations and its
grading, vanishing on the ground, the part of a state the frame can miss, unsharp
effects, the spectral construction and its converse, and normalized frame components."""
import numpy as np
from numpy.testing import assert_allclose

from lib.fock import dag
from lib.linalg import rank

LABELS = {
    "LG": ["def:derived-duality-frame", "thm:duality-frame-vanishes-ground",
           "choice:duality-generator", "q:projection-or-effect-duality-frame",
           "cor:no-independent-subject-ket"],
    "UM": ["sec:determination", "subsec:generator-families", "eq:ground-line-kept",
           "sec:neutral-graphs-phase"],
}

RNG = np.random.default_rng(5335)
N = 6


def random_complex(*shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def unit(vector):
    return vector / np.linalg.norm(vector)


def projection(vector):
    return np.outer(vector, vector.conj())


def ground_and_complement():
    y = unit(random_complex(N))
    others = np.linalg.qr(np.column_stack([y, random_complex(N, N - 1)]))[0][:, 1:]
    return y, others                                    # orthonormal basis of y-perp


def frame(others, n_plus, n_minus):
    """Projections onto the first n_plus and the next n_minus vectors orthogonal to y."""
    plus = others[:, :n_plus]
    minus = others[:, n_plus:n_plus + n_minus]
    return plus @ dag(plus), minus @ dag(minus)


def test_frame_relations_grading_and_vanishing_on_the_ground():
    y, others = ground_and_complement()
    Q = np.eye(N) - projection(y)
    for n_plus, n_minus in ((2, 2), (1, 3), (2, 0), (0, 0)):
        Pp, Pm = frame(others, n_plus, n_minus)
        Pd = Pp + Pm
        assert_allclose(Pp @ Pm, 0, atol=1e-12)
        assert_allclose(Pd @ Pd, Pd, atol=1e-12)
        assert np.linalg.eigvalsh(Q - Pd).min() > -1e-12               # P_d <= I - P_y
        G = Pp - Pm
        assert_allclose(dag(G), G, atol=1e-12)
        assert_allclose(G @ G, Pd, atol=1e-12)
        for P in (Pp, Pm, Pd):
            assert_allclose(P @ y, 0, atol=1e-12)
        # components of a state: only its differentiated part enters, and the frame can miss some
        phi = unit(random_complex(N))
        eta = Q @ phi
        assert_allclose(Pp @ phi, Pp @ eta, atol=1e-12)
        assert_allclose(Pm @ phi, Pm @ eta, atol=1e-12)
        missed = (Q - Pd) @ phi
        assert_allclose(np.linalg.norm(eta) ** 2,
                        np.linalg.norm(Pp @ phi) ** 2 + np.linalg.norm(Pm @ phi) ** 2
                        + np.linalg.norm(missed) ** 2, atol=1e-12)
        assert (np.linalg.norm(missed) > 1e-3) == (n_plus + n_minus < N - 1)
    # the example of the text: y = e_0 on C^3, P_+ onto e_1, P_- = 0
    e = np.eye(3)
    phi = (e[1] + e[2]) / np.sqrt(2)
    Pp = projection(e[1])
    Q3 = np.eye(3) - projection(e[0])
    assert_allclose(np.linalg.norm((Q3 - Pp) @ phi) ** 2, 0.5, atol=1e-12)


def test_zero_mean_grading_is_not_the_ground():
    y, others = ground_and_complement()
    Pp, Pm = frame(others, 2, 2)
    G = Pp - Pm
    balanced = unit(others[:, 0] + others[:, 2])                 # half on each aspect
    assert_allclose(np.vdot(balanced, G @ balanced), 0, atol=1e-12)
    assert_allclose(np.linalg.norm((Pp + Pm) @ balanced), 1, atol=1e-12)
    assert_allclose(np.vdot(y, G @ y), 0, atol=1e-12)
    assert_allclose(np.linalg.norm((Pp + Pm) @ y), 0, atol=1e-12)


def test_unsharp_effects_vanish_on_the_ground_too():
    y, others = ground_and_complement()
    Q = np.eye(N) - projection(y)
    for _ in range(5):
        A, B = random_complex(N, N), random_complex(N, N)
        Ep, Em = Q @ A @ dag(A) @ Q, Q @ B @ dag(B) @ Q
        scale = np.linalg.eigvalsh(Ep + Em).max()
        Ep, Em = Ep / scale, Em / scale
        for E in (Ep, Em):
            assert np.linalg.eigvalsh(E).min() > -1e-12
        assert np.linalg.eigvalsh(Q - Ep - Em).min() > -1e-10
        assert_allclose(Ep @ y, 0, atol=1e-10)
        assert_allclose(Em @ y, 0, atol=1e-10)
        # the effects are unsharp: not projections
        assert np.linalg.norm(Ep @ Ep - Ep) > 1e-3


def test_spectral_construction_and_its_converse():
    y, others = ground_and_complement()
    Q = np.eye(N) - projection(y)
    M = random_complex(N, N)
    D = Q @ (M + dag(M)) @ Q                                    # self-adjoint with D y = 0
    assert_allclose(D @ y, 0, atol=1e-12)
    values, vectors = np.linalg.eigh(D)
    pos = vectors[:, values > 1e-10]
    neg = vectors[:, values < -1e-10]
    Pp, Pm = pos @ dag(pos), neg @ dag(neg)
    assert_allclose(Pp @ Pm, 0, atol=1e-12)
    kernel = vectors[:, np.abs(values) <= 1e-10]
    assert_allclose(Pp + Pm, np.eye(N) - kernel @ dag(kernel), atol=1e-10)
    assert np.linalg.eigvalsh(Q - Pp - Pm).min() > -1e-10
    sign = vectors @ np.diag(np.sign(np.where(np.abs(values) > 1e-10, values, 0))) @ dag(vectors)
    assert_allclose(sign, Pp - Pm, atol=1e-10)
    # converse: every sharp frame is the spectral frame of its own grading
    Pp2, Pm2 = frame(others, 1, 3)
    G = Pp2 - Pm2
    values, vectors = np.linalg.eigh(G)
    assert_allclose(vectors[:, values > 0.5] @ dag(vectors[:, values > 0.5]), Pp2, atol=1e-10)
    assert_allclose(vectors[:, values < -0.5] @ dag(vectors[:, values < -0.5]), Pm2, atol=1e-10)
    assert_allclose(G @ y, 0, atol=1e-12)


def test_normalized_component_depends_on_state_and_frame():
    y, others = ground_and_complement()
    Q = np.eye(N) - projection(y)
    Pp, Pm = frame(others, 2, 2)
    Pd = Pp + Pm
    phi = unit(random_complex(N))
    component = Pp @ phi
    probability = np.linalg.norm(component) ** 2
    conditional = component / np.linalg.norm(component)
    # the Lueders operation of the outcome + of {P_+, P_-, I - P_d}
    rho = projection(phi)
    after = Pp @ rho @ Pp
    assert_allclose(np.trace(after).real, probability, atol=1e-12)
    assert_allclose(after / probability, projection(conditional), atol=1e-12)
    assert_allclose(sum(np.trace(P @ rho).real for P in (Pp, Pm, np.eye(N) - Pd)), 1, atol=1e-12)
    # another state, or another frame, gives another normalized component
    other_state = unit(random_complex(N))
    v = Pp @ other_state
    assert abs(abs(np.vdot(unit(v), conditional)) - 1) > 1e-3
    rotation = np.linalg.qr(random_complex(N - 1, N - 1))[0]
    rotated = others @ rotation
    Pp_rot = rotated[:, :2] @ dag(rotated[:, :2])
    w = Pp_rot @ phi
    assert abs(abs(np.vdot(unit(w), conditional)) - 1) > 1e-3
    # at the ground there is no component to normalize
    assert_allclose(Pp @ y, 0, atol=1e-12)
    assert rank(Pp) == 2
