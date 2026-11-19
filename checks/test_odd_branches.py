"""Signed orientations, adjoints, natural domains, and polar factors of a closed channel."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import block_diag, expm

from lib.fock import dag


LABELS = {
    "TAUM": ["eq:odd", "lem:odd-adjoint", "thm:classification", "prop:canonical-sector",
             "eq:balanced-branches", "thm:AV-algebra", "thm:dichotomy", "cor:common-magnitude",
             "prop:cycle-positivity", "prop:Omega", "thm:polar-factorization",
             "prop:normalized-cycle", "app:domains"],
    "LG": ["def:odd-extensions", "prop:imported-algebra", "thm:nilpotent-isometry-obstruction",
           "prop:typed-isometry", "cor:typed-odd-separation", "prop:isometric-asymmetry",
           "prop:no-fixed-vector"],
    "UM": ["eq:carrier-grading", "prop:odd-block-adjoints", "eq:oriented-extensions",
           "eq:oriented-adjoints", "eq:oriented-products", "eq:oriented-brackets", "prop:odd-channel-branches",
           "eq:weighted-channel-branches", "eq:branch-squares", "eq:branch-return-form",
           "eq:branch-reference-image", "eq:balanced-channel-branches", "eq:odd-scalar-normal-form",
           "eq:branch-polar-factorization", "eq:inverse-phase-domain", "eq:support-normalized-pair",
           "eq:support-normalized-products", "eq:full-polar-rotation", "eq:typed-odd-polar",
           "prop:typed-odd-isometry", "eq:typed-isometric-range", "eq:active-laplacian-gap",
           "eq:graded-trace-general"],
}


def differential(channel):
    n_out, n_in = channel.shape
    d = np.zeros((n_in + n_out, n_in + n_out), dtype=complex)
    d[n_in:, :n_in] = channel
    return d


def test_signed_orientations_and_classification():
    # Recover the lower channel from a general odd, Krein-self-adjoint matrix.
    rng = np.random.default_rng(784)
    n_in, n_out = 3, 5
    size = n_in + n_out
    J = block_diag(np.eye(n_in), -np.eye(n_out))
    M = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    odd = (M - J @ M @ J) / 2
    K = (odd + J @ dag(odd) @ J) / 2
    A = -K[n_in:, :n_in]
    forward = -differential(A)
    backward = dag(differential(A))
    D = backward - forward
    delta = block_diag(dag(A) @ A, A @ dag(A))
    assert_allclose(K, forward + backward, atol=1e-13)
    assert_allclose(dag(forward), -backward)
    assert_allclose(dag(backward), -forward)
    assert_allclose(J @ dag(forward) @ J, backward)
    assert_allclose(J @ dag(backward) @ J, forward)
    assert_allclose(forward @ forward, 0)
    assert_allclose(backward @ backward, 0)
    assert_allclose(backward @ forward, block_diag(-dag(A) @ A, np.zeros((n_out, n_out))), atol=1e-13)
    assert_allclose(forward @ backward, block_diag(np.zeros((n_in, n_in)), -A @ dag(A)), atol=1e-13)
    assert_allclose(forward @ backward - backward @ forward, J @ delta, atol=1e-13)
    assert_allclose(K, J @ D)
    for branch, sign in ((D, 1), (K, -1)):
        assert_allclose(dag(branch), sign * branch)
        assert_allclose(J @ dag(branch) @ J, -sign * branch)
        assert_allclose(branch @ branch, sign * delta, atol=1e-13)
        assert_allclose(dag(branch) @ branch, delta, atol=1e-13)


def test_scalar_normal_form_and_balanced_branches():
    sigma_x = np.array([[0, 1], [1, 0]])
    i_sigma_y = np.array([[0, 1], [-1, 0]])
    J = np.diag([1, -1])
    for a in (0, -2, 2j, 1 - 3j):
        signed = np.array([[0, a], [-np.conj(a), 0]])
        gauge = np.diag(np.exp(0.5j * np.angle(a) * np.array([-1, 1])))
        assert_allclose(dag(gauge) @ J @ gauge, J, atol=1e-14)
        assert_allclose(gauge @ signed @ dag(gauge), abs(a) * i_sigma_y, atol=1e-14)
    # A self-adjoint channel can have both signs. Its modulus is not the original channel.
    channel = np.array([[1, 1j], [-1j, -1]])
    d = differential(channel)
    assert_allclose(d + dag(d), np.kron(sigma_x, channel))
    assert_allclose(dag(d) - d, np.kron(i_sigma_y, channel))
    assert_allclose(dag(channel) @ channel, 2 * np.eye(2))


def test_weighted_and_support_normalized_branches():
    # Two active modes of unequal magnitude; both sectors also have a kernel.
    # Fourier bases give genuinely complex singular vectors on the two different carriers.
    input_modes = np.fft.fft(np.eye(3), norm="ortho")[:, [1, 2]]
    output_modes = np.fft.fft(np.eye(4), norm="ortho")[:, [1, 3]]
    weights = np.diag([0.5, 3.0])
    A = output_modes @ weights @ dag(input_modes)
    phase = output_modes @ dag(input_modes)
    input_support = input_modes @ dag(input_modes)
    output_support = output_modes @ dag(output_modes)
    support = block_diag(input_support, output_support)
    magnitude = block_diag(input_modes @ weights @ dag(input_modes),
                           output_modes @ weights @ dag(output_modes))
    d = differential(A)
    signed_forward = -differential(phase)
    signed_backward = dag(differential(phase))
    reflection = signed_backward - signed_forward
    rotation = signed_backward + signed_forward
    assert_allclose(-signed_backward @ signed_forward,
                    block_diag(input_support, np.zeros((4, 4))), atol=1e-13)
    assert_allclose(-signed_forward @ signed_backward,
                    block_diag(np.zeros((3, 3)), output_support), atol=1e-13)
    for branch, polar, sign in ((d + dag(d), reflection, 1), (dag(d) - d, rotation, -1)):
        assert_allclose(branch, polar @ magnitude, atol=1e-13)
        assert_allclose(branch, magnitude @ polar, atol=1e-13)
        assert_allclose(branch @ np.linalg.pinv(magnitude), polar, atol=1e-13)
        assert_allclose(polar @ polar, sign * support, atol=1e-13)
        assert_allclose(polar @ (np.eye(7) - support), 0, atol=1e-13)
    active_basis = block_diag(input_modes, output_modes)
    active_return = dag(active_basis) @ (magnitude @ magnitude) @ active_basis
    assert_allclose(np.linalg.eigvalsh(active_return), [0.25, 0.25, 9, 9], atol=1e-13)
    # The graded trace formula applies even when neither block is a density or heat operator.
    input_block = (1 + 1j) * input_support
    output_block = -2j * output_support
    J = block_diag(np.eye(3), -np.eye(4))
    assert_allclose(np.trace(J @ block_diag(input_block, output_block)),
                    np.trace(input_block) - np.trace(output_block))
    t = 0.73
    assert_allclose(expm(t * rotation), np.eye(7) - support + np.cos(t) * support + np.sin(t) * rotation,
                    atol=1e-13)

    # The oscillator's weighted rotation runs at a different speed in each mode.
    # Its phase is the ordinary shift, giving speed one on the same planes.
    levels = 4
    shift = np.eye(levels + 1)[:, 1:]
    creation = shift @ np.diag(np.sqrt(np.arange(1, levels + 1)))
    weighted_arrow = differential(creation)
    polar_arrow = differential(shift)
    weighted_evolution = expm(t * (dag(weighted_arrow) - weighted_arrow))
    polar_evolution = expm(t * (dag(polar_arrow) - polar_arrow))
    basis = np.eye(2 * levels + 1)
    for n in (1, 3):
        q = basis[n]
        p = basis[levels + n + 1]
        angle = t * np.sqrt(n + 1)
        assert_allclose(weighted_evolution @ q, np.cos(angle) * q - np.sin(angle) * p, atol=1e-13)
        assert_allclose(polar_evolution @ q, np.cos(t) * q - np.sin(t) * p, atol=1e-13)


def test_first_order_and_ordered_return_domains():
    # The off-diagonal expression X(f,g)=(g',f') anticommutes formally with
    # J. A pair of affine functions with f(0)=g(0)=1 lies in the stated H¹
    # domain, whereas its J-image fails the boundary condition.
    trace = np.array([1, 1])
    slopes = np.array([2.0, -1j])
    J = np.diag([1, -1])
    assert_allclose(J @ slopes[::-1], -(J @ slopes)[::-1])
    assert trace[0] == trace[1]
    transformed_trace = J @ trace
    assert transformed_trace[0] != transformed_trace[1]

    # A e_n = 2**n e_n, n >= 0. The exact partial sums distinguish Hilbert norm,
    # first-order graph energy, and existence of a Hilbert-space return vector.
    for count in (4, 8, 16):
        singular = 2.0 ** np.arange(count)
        A = np.diag(singular)
        d = differential(A)
        x = singular ** -2
        y = singular ** -1
        zero = np.zeros(count)
        reference = np.r_[x, zero]
        image = -d @ reference
        assert_allclose(image, np.r_[zero, -y])
        assert_allclose(np.vdot(image, image).real, (1 - 4.0 ** -count) / (1 - 0.25))
        assert_allclose(np.vdot(image, np.r_[zero, -image[count:]]), -np.vdot(image, image))
        assert_allclose(np.vdot(image, image), np.vdot(x, dag(A) @ A @ x))
        assert_allclose(d @ image, 0)
        # Each cutoff adds one unit to the squared return norm: x is in Dom A but
        # outside Dom A*A. A second forward step is defined even in the infinite model.
        assert_allclose(dag(d) @ image, np.r_[-np.ones(count), zero])
        assert_allclose(np.linalg.norm(dag(d) @ image) ** 2, count)
        # A spectator y outside Dom A* is allowed in Dom d*d; the other ordered
        # product first applies d*, whose squared norm already diverges for this y.
        spectator = np.r_[zero, y]
        assert_allclose(dag(d) @ (d @ spectator), 0)
        assert_allclose(np.linalg.norm(dag(d) @ spectator) ** 2, count)
        # The roles of the two sectors reverse for dd*: the inactive input
        # component is unrestricted; the reverse orientation squares to zero.
        assert_allclose(d @ (dag(d) @ np.r_[y, zero]), 0)
        assert_allclose(dag(d) @ (dag(d) @ np.r_[zero, x]), 0)


def test_unbounded_inverse_and_extension_of_the_polar_phases():
    # The manuscript's A e_n = e_n/(n+1) has dense range in the active output.
    # The inverse applied to z_n = 1/(n+1) produces a nonsummable constant sequence.
    for count in (4, 8, 16):
        singular = 1 / np.arange(1, count + 1)
        d = differential(np.diag(singular))
        inverse = np.tile(1 / singular, 2)
        target = np.r_[singular, np.zeros(count)]
        preimage = inverse * target
        assert_allclose(np.vdot(preimage, preimage), count)
        phase_image = np.r_[np.zeros(count), singular]
        for branch, sign in ((d + dag(d), 1), (dag(d) - d, -1)):
            assert_allclose(branch @ preimage, sign * phase_image)
            assert_allclose(np.linalg.norm(phase_image), np.linalg.norm(target))
        # Integral bounds on the tail of sum 1/n² apply for every cutoff.
        # The truncated targets and phase images both converge in Hilbert norm.
        squared_tail = np.pi ** 2 / 6 - np.vdot(phase_image, phase_image).real
        assert 1 / (count + 1) <= squared_tail <= 1 / count


def test_typed_isometry_and_global_nilpotency():
    # An isometry into a proper output subspace: the image omits the first site.
    embedding = np.eye(4)[:, 1:]
    d = differential(embedding)
    assert_allclose(dag(embedding) @ embedding, np.eye(3))
    assert_allclose(embedding @ dag(embedding), np.diag([0, 1, 1, 1]))
    v = np.array([0, 1, 1j]) / np.sqrt(2)
    input_vector = np.r_[v, np.zeros(4)]
    output_vector = np.r_[np.zeros(3), embedding @ v]
    assert_allclose(-d @ input_vector, -output_vector)
    assert_allclose(np.linalg.norm(output_vector), np.linalg.norm(input_vector))
    for arrow in (-d, dag(d)):
        assert_allclose(arrow @ arrow, 0)
        assert np.linalg.matrix_rank(arrow - np.eye(7)) == 7
        assert not np.allclose(dag(arrow) @ arrow, np.eye(7))
    assert_allclose(d @ output_vector, 0)
