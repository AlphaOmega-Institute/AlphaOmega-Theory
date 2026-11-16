"""Generation and readout in the five-site path example."""
import numpy as np
from numpy.testing import assert_allclose

LABELS = {
    "LG": ["thm:concrete-GNS-appearance", "prop:covariance-content-inner-product"],
    "UM": ["prop:cyclic-determination-content", "eq:path-readout-states"],
}


def test_cell_and_site_determinations_on_the_path():
    ground = np.ones(5) / np.sqrt(5)
    site_indicators = np.eye(5)
    cell_indicators = np.column_stack((site_indicators[:, 0], site_indicators[:, 1:].sum(axis=1)))
    copying_map = cell_indicators / 2
    direct_projection = copying_map @ np.linalg.pinv(copying_map)

    for indicators, dimension in ((cell_indicators, 2), (site_indicators, 5)):
        # Each column is a diagonal multiplication operator applied to the ground.
        orbit = ground[:, None] * indicators
        means = ground @ orbit
        centered_orbit = orbit - np.outer(ground, means)
        covariance = orbit.T @ orbit - np.outer(means, means)
        assert np.linalg.matrix_rank(orbit) == dimension
        assert np.linalg.matrix_rank(centered_orbit) == dimension - 1
        assert_allclose(centered_orbit.T @ centered_orbit, covariance)

    cell_orbit = ground[:, None] * cell_indicators
    assert_allclose(cell_orbit @ np.linalg.pinv(cell_orbit), direct_projection, atol=1e-14)
    prepared = np.array([0, 1, 1, 1, 1]) / 2
    defect = np.array([0, 1, 1, -1, -1]) / 2
    assert_allclose(abs(prepared) ** 2, abs(defect) ** 2)
    assert_allclose(prepared @ direct_projection @ prepared, 1)
    assert_allclose(defect @ direct_projection @ defect, 0, atol=1e-14)
