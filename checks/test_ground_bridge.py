"""Selected ground lines and the partial isometry in the polar decomposition."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import polar

LABELS = {
    "LG": ["thm:ground-channel-polar", "thm:strong-ground-equivalences"],
    "UM": ["prop:ground-polar-compatibility"],
}


def test_raw_ground_bridge_and_polar_ground_measurement():
    channel = np.array([[0.5, 0], [0, 1], [0, 0]])
    orientation, weight = polar(channel)
    # The same channel with two choices of ground lines. Both satisfy the raw bridge.
    cases = (
        (np.array([1, 1]) / np.sqrt(2), np.array([1, 2, 0]) / np.sqrt(5),
         np.sqrt(5 / 8), 0.9, False),
        (np.array([1, 0]), np.array([1, 0, 0]), 0.5, 1.0, True),
    )
    for g, h, gain, probability, strong in cases:
        assert_allclose(channel @ g, gain * h)
        input_ground = np.outer(g, g)
        output_ground = np.outer(h, h)
        line_orientation = channel @ input_ground / gain
        assert_allclose(line_orientation, np.outer(h, g))

        polar_preimage = orientation.T @ h
        assert_allclose(polar_preimage, weight @ g / gain)
        detector = orientation.T @ output_ground @ orientation
        assert_allclose(detector, np.outer(polar_preimage, polar_preimage))
        assert_allclose(g @ detector @ g, probability)
        for actual, expected in (
            (channel.T @ channel @ g, gain**2 * g),
            (weight @ g, gain * g),
            (orientation @ g, h),
            (orientation @ input_ground, line_orientation),
            (channel.T @ h, gain * g),
        ):
            assert np.allclose(actual, expected) == strong


def test_normalizing_polar_image_can_hide_a_kernel_component():
    # A projection is its own positive factor and its canonical polar partial isometry.
    orientation = np.diag([1, 0])
    g = np.array([1, 1]) / np.sqrt(2)
    h = np.array([1, 0])
    image = orientation @ g
    assert_allclose(np.linalg.norm(image) ** 2, 0.5)
    assert_allclose(image / np.linalg.norm(image), h)
    assert not np.allclose(orientation.T @ h, g)
