"""Signed Schur reduction with bounded coupling and self-adjoint diagonal blocks."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.linalg import block_diag

from lib.fock import dag
from lib.linalg import rank, response_weights


LABELS = {
    "TAUM": ["ass:block-generator", "prop:block-branches", "def:effective-operator",
             "thm:Schur-factorization", "prop:oriented-self-energy", "eq:Leff-plus",
             "thm:self-energy-Herglotz", "cor:effective-imaginary-sign", "cor:static-self-energy"],
    "UM": ["eq:signed-block-generator", "eq:signed-schur-pencil", "prop:signed-schur",
           "eq:signed-schur-factorization", "eq:signed-block-inverse", "eq:signed-oriented-excursion",
           "eq:signed-self-energy-imaginary", "eq:signed-effective-imaginary",
           "eq:signed-static-return", "eq:signed-corner-imaginary", "cor:signed-full-corner",
           "eq:passive-signed-witness"],
}


def example_blocks():
    # Rectangular coupling with two distinct channels and one uncoupled retained direction.
    retained = np.array([[1, 1j, 0], [-1j, 2, 0.5], [0, 0.5, -1]])
    complement = np.array([[2, 1 + 1j], [1 - 1j, 4]])
    coupling = np.array([[1, 1j, 0], [1j, 0.5, 0]])
    return retained, complement, coupling


def generator(retained, complement, coupling, sign):
    return np.block([[retained, dag(coupling)], [sign * coupling, complement]])


def resolvent(operator, z):
    return np.linalg.inv(operator - z * np.eye(len(operator)))


def imaginary_part(operator):
    return (operator - dag(operator)) / (2j)


def test_factorization_and_four_inverse_blocks():
    retained, complement, coupling = example_blocks()
    ip, iq = np.eye(3), np.eye(2)
    zero = np.zeros_like(coupling)
    grading = block_diag(ip, -iq)
    for sign in (1, -1):
        full = generator(retained, complement, coupling, sign)
        adjoint_metric = block_diag(ip, sign * iq)
        assert_allclose(adjoint_metric @ dag(full) @ adjoint_metric, full)
        if sign == -1:
            assert not np.allclose(dag(full), full)
            assert_allclose(grading @ dag(full) @ grading, full)
        for z in (-3.0, 0.7 + 0.9j, 1.2 - 1.1j):
            r = resolvent(complement, z)
            pencil = retained - z * ip - sign * dag(coupling) @ r @ coupling
            upper = np.block([[ip, dag(coupling) @ r], [zero, iq]])
            lower = np.block([[ip, dag(zero)], [sign * r @ coupling, iq]])
            diagonal = block_diag(pencil, complement - z * iq)
            assert_allclose(upper @ diagonal @ lower, full - z * np.eye(5), atol=2e-14)
            t = np.linalg.inv(pencil)
            direct_inverse = resolvent(full, z)
            reconstructed = np.block([
                [t, -t @ dag(coupling) @ r],
                [-sign * r @ coupling @ t, r + sign * r @ coupling @ t @ dag(coupling) @ r],
            ])
            assert_allclose(reconstructed, direct_inverse, atol=2e-14)


def test_triangular_factors_on_the_unbounded_diagonal_domain():
    # In the infinite model both diagonal generators are N e_j = j e_j (j>=1).
    # C=|v><v| with v_j=1/j is bounded. The left factor may send smooth
    # right-hand data outside Dom N; R_Q C in the right factor restores the domain.
    z = 2j
    for count in (4, 16, 32):
        levels = np.arange(1, count + 1)
        n = np.diag(levels)
        e = np.eye(count)[0]
        v = 1 / levels
        c = np.outer(v, v)
        r = resolvent(n, z)
        data = np.r_[np.zeros(count), e]
        rough_input = dag(c) @ r @ e
        assert_allclose(np.vdot(n @ rough_input, n @ rough_input), count / abs(1 - z)**2)
        assert_allclose(np.linalg.norm(n @ c @ e)**2, count)
        assert_allclose(np.linalg.norm(n @ r @ c @ e)**2, np.sum(1 / (levels**2 + 4)))
        assert np.linalg.norm(n @ r @ c @ e)**2 < np.pi**2 / 6
        for sign in (1, -1):
            solution = resolvent(generator(n, n, c, sign), z) @ data
            p, q = np.split(solution, 2)
            assert_allclose(q, r @ (e - sign * c @ p), atol=1e-14)
            assert_allclose(p, -r @ dag(c) @ q, atol=1e-14)
            # ||C|| < pi²/6 and ||N R_Q(2i)|| <= 1. Bounded perturbation
            # gives a uniform solution norm, and the block equations then give
            # graph-norm bounds even though C e_1 lies outside Dom N.
            bound_c = np.pi**2 / 6
            assert np.linalg.norm(solution) <= 1 / (z.imag - bound_c) + 1e-13
            assert np.linalg.norm(n @ p) <= bound_c * np.linalg.norm(q) + 1e-13
            assert np.linalg.norm(n @ q) <= 1 + bound_c * np.linalg.norm(p) + 1e-13


def test_pencil_zeros_lift_and_elimination_can_fail_at_a_regular_point():
    # The balanced pair has lambda = delta/2 on both diagonal entries and
    # coupling -delta/2. At delta=2 the Hilbert levels are 0,2; flipping the
    # lower coupling alone gives a signed comparator with poles at 1 +/- i.
    # A second choice is the zero-diagonal example in the signed reduction.
    for center, coupling in ((1.0, -1.0), (0.0, 2.0)):
        for sign, roots in ((1, (-1, 1)), (-1, (-1j, 1j))):
            full = np.array([[center, coupling], [sign * coupling, center]])
            for root in roots:
                z = center + root * abs(coupling)
                pencil = center - z - sign * coupling**2 / (center - z)
                eigenvector = np.array([1, -sign * coupling / (center - z)])
                assert_allclose(pencil, 0, atol=1e-15)
                assert_allclose(full @ eigenvector, z * eigenvector)
            # The eliminated resolvent has a pole here, while the full inverse exists.
            assert_allclose(resolvent(full, center)[0, 0], 0, atol=1e-15)


def test_oriented_return_imaginary_forms_and_kernel():
    retained, complement, coupling = example_blocks()
    z = 1.7 + 0.6j
    r = resolvent(complement, z)
    sigma = dag(coupling) @ r @ coupling
    differential = np.zeros((5, 5), dtype=complex)
    differential[3:, :3] = coupling
    forward, backward = -differential, dag(differential)
    embedded = block_diag(np.zeros((3, 3)), r)
    assert_allclose(-backward @ embedded @ forward, block_diag(sigma, np.zeros((2, 2))))
    atoms = response_weights(complement, coupling)
    assert_allclose(sigma, sum(weight / (energy - z) for energy, weight in atoms), atol=2e-14)
    assert_allclose(dag(sigma), dag(coupling) @ resolvent(complement, z.conjugate()) @ coupling)
    imaginary = imaginary_part(sigma)
    assert_allclose(imaginary, z.imag * dag(r @ coupling) @ (r @ coupling), atol=1e-14)
    assert rank(imaginary) == rank(coupling) == 2
    assert_allclose(imaginary[:, 2], 0, atol=1e-14)
    for sign in (1, -1):
        effective = retained - sign * sigma
        assert_allclose(imaginary_part(effective), -sign * imaginary, atol=1e-14)
        assert_allclose(imaginary_part(effective - z * np.eye(3)),
                        -sign * imaginary - z.imag * np.eye(3), atol=1e-14)


def test_below_threshold_static_shift_and_norm_bound():
    retained, complement, coupling = example_blocks()
    bottom = np.linalg.eigvalsh(complement)[0]
    for energy in (-3.0, 0.3, bottom - 0.04):
        r = resolvent(complement, energy)
        sigma = dag(coupling) @ r @ coupling
        eigenvalues = np.linalg.eigvalsh(sigma)
        assert_allclose(sigma, dag(sigma), atol=1e-14)
        assert_allclose(eigenvalues[0], 0, atol=1e-14)
        assert eigenvalues[1] > 0
        assert eigenvalues[-1] <= np.linalg.norm(coupling, 2)**2 / (bottom - energy) + 1e-12
        for sign in (1, -1):
            effective = retained - sign * sigma
            for h in (np.array([1, 1j, 0]), np.array([0, 0, 1])):
                shift = np.vdot(h, (effective - retained) @ h)
                assert_allclose(shift, -sign * np.vdot(coupling @ h, r @ coupling @ h), atol=1e-13)
    # For H_Q e_n=n e_n, n>=1, C=I and energy zero, every nonzero vector has
    # a positive return form, yet no positive lower bound works at all cutoffs.
    for count in (4, 16, 64):
        levels = np.arange(1, count + 1)
        inverse_levels = resolvent(np.diag(levels), 0)
        assert_allclose(np.linalg.eigvalsh(inverse_levels)[[0, -1]], [1 / count, 1])


def test_krein_corner_sign_and_real_spectrum_witness():
    retained, complement, coupling = example_blocks()
    z = 1.5 + 0.9j
    r = resolvent(complement, z)
    threshold = np.linalg.norm(r @ coupling, 2)
    for ratio in (0.0, 0.5, 1.0, 1.8):
        c = coupling * ratio / threshold
        gram = dag(r @ c) @ (r @ c)
        for sign in (1, -1):
            t = resolvent(generator(retained, complement, c, sign), z)[:3, :3]
            imaginary = imaginary_part(t)
            assert_allclose(imaginary, z.imag * dag(t) @ (np.eye(3) + sign * gram) @ t, atol=1e-13)
            least = np.linalg.eigvalsh(imaginary)[0]
            if sign == 1 or ratio < 1:
                assert least > 1e-8
            elif ratio == 1:
                assert_allclose(least, 0, atol=1e-13)
            else:
                assert least < -1e-8
    # Section 9's signed witness has real eigenvalues although its reciprocal
    # corner has a negative spectral weight and a negative imaginary part here.
    z = 2 + 0.5j
    full = np.array([[0, np.sqrt(2)], [-np.sqrt(2), 3]])
    assert_allclose(np.sort(np.linalg.eigvals(full)), [1, 2])
    assert_allclose(abs(np.sqrt(2) / (3 - z))**2, 8 / 5)
    corner = resolvent(full, z)[0, 0]
    assert_allclose(corner, 2 / (1 - z) - 1 / (2 - z))
    assert_allclose(corner.imag, -6 / 5)

    # An additional retained level at 2 cancels that negative weight for this
    # scalar port, even though preparation through the port couples to Q.
    retained = np.diag([0.0, 2.0])
    complement = np.array([[3.0]])
    c = np.array([[np.sqrt(2), 0.0]])
    port = np.ones(2) / np.sqrt(2)
    signed = generator(retained, complement, c, -1)
    assert_allclose(c @ port, [1.0])
    for z in (2 + 0.5j, 1 + 0.01j, 3 - 0.8j, -2 + 0.4j):
        retained_corner = resolvent(signed, z)[:2, :2]
        assert_allclose(np.vdot(port, retained_corner @ port), 1 / (1 - z), atol=1e-13)
        if z == 2 + 0.5j:
            assert_allclose(np.linalg.eigvalsh(imaginary_part(retained_corner)), [-6 / 5, 2])


def test_resolvent_growth_near_complementary_spectral_mass():
    # At the complementary level of the §9 witness, C sends a unit input
    # to a vector of squared norm 2. The atom gives a pole of order one.
    complement = np.array([[3.0]])
    coupled_input = np.array([np.sqrt(2)])
    # For the atomless comparison take H_Q f(s) = s f(s) and C1(s) = 1
    # in L²([-1, 1], ds/2). Even with no eigenvector at 0, the norm grows.
    for eta in (0.25, 0.125, 0.0625, 0.03125):
        at_atom = resolvent(complement, 3 + 1j * eta) @ coupled_input
        assert_allclose(eta * np.linalg.norm(at_atom), np.sqrt(2))

        def squared_multiplier(s):
            return 0.5 / (s * s + eta * eta)

        squared_norm = quad(squared_multiplier, -1, 1)[0]
        assert_allclose(eta * squared_norm, np.arctan(1 / eta), rtol=1e-12)
        # The interval [-eta, eta] has mass eta. On it the squared
        # resolvent multiplier is at least 1/(2 eta²).
        assert squared_norm >= 1 / (2 * eta) > 1
