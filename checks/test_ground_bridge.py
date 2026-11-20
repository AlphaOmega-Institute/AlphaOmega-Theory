"""Selected ground lines and the partial isometry in the polar decomposition."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import polar

LABELS = {
    "LG": ["thm:ground-channel-polar", "thm:strong-ground-equivalences",
           "prop:normalized-volition-inverse", "eq:raw-volition-vector-return"],
    "UM": ["prop:ground-polar-compatibility", "eq:ground-return-energy",
           "subsec:volition-returns", "subsec:volition-actualization",
           "eq:weighted-volition", "eq:jug-physical-return"],
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


def test_path_ground_measurement_depends_on_cell_sizes():
    def path_hamiltonian(site_count):
        differences = np.diff(np.eye(site_count), axis=0)
        return differences.T @ differences / 2

    def constant_ground(hamiltonian):
        energies, vectors = np.linalg.eigh(hamiltonian)
        expected = np.ones(len(energies)) / np.sqrt(len(energies))
        assert_allclose(energies[0], 0, atol=1e-14)
        assert energies[1] > 0
        assert_allclose(np.outer(vectors[:, 0], vectors[:, 0]), np.outer(expected, expected))
        return expected

    coarse_hamiltonian = path_hamiltonian(2)
    coarse_ground = constant_ground(coarse_hamiltonian)
    cases = (((1, 4), 0.9), ((2, 2), 1.0))
    for cell_sizes, expected_probability in cases:
        fine_ground = constant_ground(path_hamiltonian(sum(cell_sizes)))
        channel = np.repeat(np.eye(2), cell_sizes, axis=0) / 2
        gain = np.sqrt(sum(cell_sizes) / 8)
        assert_allclose(channel @ coarse_ground, gain * fine_ground)

        orientation, _ = polar(channel)
        partner = orientation.T @ fine_ground
        assert_allclose(partner, np.sqrt(cell_sizes) / np.sqrt(sum(cell_sizes)))
        assert_allclose(abs(coarse_ground @ partner) ** 2, expected_probability)
        assert_allclose(partner @ coarse_hamiltonian @ partner, 1 - expected_probability, atol=1e-14)
        if cell_sizes == (1, 4):
            weighted = channel.T @ fine_ground
            assert_allclose(weighted, np.array([1, 4]) / (2 * np.sqrt(5)))
            assert_allclose(weighted @ weighted, 17 / 20)
            assert_allclose(abs(coarse_ground @ weighted) ** 2 / (weighted @ weighted), 25 / 34)

    # On the path with unequal cells, the two judgement operations of §16
    # change the support weight even though their conditional returns agree.
    copying_map = np.repeat(np.eye(2), [1, 4], axis=0) / 2
    orientation, _ = polar(copying_map)
    e1, e2, e3, e4, _ = np.eye(5)
    f, h = (e1 - e2) / np.sqrt(2), (e3 - e4) / np.sqrt(2)
    assert_allclose(orientation.T @ f, np.array([1, -0.5]) / np.sqrt(2), atol=1e-14)
    assert_allclose(copying_map.T @ f, np.array([1, -1]) / (2 * np.sqrt(2)), atol=1e-14)
    assert_allclose(orientation.T @ h, 0, atol=1e-14)
    assert_allclose(copying_map.T @ h, 0, atol=1e-14)
    cases = (
        ((f + h) / np.sqrt(2), 5 / 16),
        ((f + 2 * h) / np.sqrt(5), 1 / 8),
        ((2 * f + h) / np.sqrt(5), 1 / 2),
    )
    for appearance, expected_weight in cases:
        partner = orientation.T @ appearance
        weighted = copying_map.T @ appearance
        assert_allclose(partner @ partner, expected_weight)
        assert_allclose(partner / np.linalg.norm(partner), np.array([2, -1]) / np.sqrt(5))
        assert_allclose(weighted / np.linalg.norm(weighted), np.array([1, -1]) / np.sqrt(2))
