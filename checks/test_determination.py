"""Generation and readout in the five-site path example."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm

LABELS = {
    "LG": ["thm:concrete-GNS-appearance", "prop:covariance-content-inner-product"],
    "UM": ["prop:cyclic-determination-content", "eq:path-readout-states", "prop:cyclic-common-excess"],
}

PREPARED = np.array([0, 1, 1, 1, 1]) / 2
DEFECT = np.array([0, 1, 1, -1, -1]) / 2


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
    assert_allclose(direct_projection @ ground, ground)
    # Projecting the cyclic ground onto the second site takes it outside the cell span.
    second_site = site_indicators[:, 1] * ground
    outside = second_site - direct_projection @ second_site
    assert_allclose(np.vdot(outside, outside), 3 / 20)
    assert_allclose(abs(PREPARED) ** 2, abs(DEFECT) ** 2)
    assert_allclose(PREPARED @ direct_projection @ PREPARED, 1)
    assert_allclose(DEFECT @ direct_projection @ DEFECT, 0, atol=1e-14)


def test_free_wait_reveals_the_sign_at_site_three():
    differences = np.diff(np.eye(5), axis=0)
    hamiltonian = differences.T @ differences / 2
    site_three = np.diag([0, 0, 1, 0, 0])

    # The Heisenberg effect has second derivative -[H,[H,P]] at zero.
    commutator = hamiltonian @ site_three - site_three @ hamiltonian
    second_derivative = commutator @ hamiltonian - hamiltonian @ commutator
    assert_allclose(PREPARED @ second_derivative @ PREPARED, 1 / 8)
    assert_allclose(DEFECT @ second_derivative @ DEFECT, -1 / 8)

    tau = 0.01
    wait = expm(-1j * tau * hamiltonian)
    probabilities = abs(wait[2] @ np.column_stack((PREPARED, DEFECT))) ** 2
    assert_allclose((probabilities[0] - probabilities[1]) / tau**2, 1 / 8, rtol=1e-4)
