"""Finite blocks of the grounded shift family and finite input extensions."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm

from lib.linalg import rank

LABELS = {
    "UM": ["eq:index-shift-family", "eq:index-extension-counts"],
    "LG": ["eq:unilateral-shift-example", "eq:shift-index-minus-one"],
}


def test_shift_defect_can_rotate_while_the_ground_stays_fixed():
    # All remaining columns of the infinite shift pair the two tails unitarily.
    g0, g1, g2 = np.eye(3)
    shift = np.column_stack((g1, g2))
    generator = np.outer(g2, g0) - np.outer(g0, g2)
    for angle in (1e-6, np.pi / 2):
        rotated = expm(angle * generator) @ shift
        assert_allclose(rotated[:, 0], g1)
        assert_allclose(rotated.T @ rotated, np.eye(2), atol=1e-14)
        assert rank(rotated - shift) == 1
        assert_allclose(np.linalg.norm(rotated - shift, 2), 2 * np.sin(angle / 2))

        defect = np.cos(angle) * g0 + np.sin(angle) * g2
        assert_allclose(np.eye(3) - rotated @ rotated.T, np.outer(defect, defect), atol=1e-14)
        # Any common defect must be on the old line g0, which has nonzero
        # overlap with the rotated image of f1 for both choices of angle.
        assert_allclose(g0 @ rotated[:, 1] / np.sin(angle), -1)


def test_input_extension_counts_independent_new_output_directions():
    g0, g1, g2, g3 = np.eye(4)
    channel = np.column_stack((g0, 2 * g0, g1))
    # The original kernel has dimension one and the defect spans g2 and g3.
    old_defect_coordinates = np.vstack((g2, g3))
    cases = (
        ((g0 + g2,), 1),
        ((g0 + g2, g1 + 2 * g2), 1),
        ((g0 + g2, g1 + g3), 2),
    )
    for columns, projected_rank in cases:
        added = np.column_stack(columns)
        assert rank(old_defect_coordinates @ added) == projected_rank
        extended = np.column_stack((channel, added))
        image_dimension = rank(extended)
        assert extended.shape[0] - image_dimension == 2 - projected_rank
        assert extended.shape[1] - image_dimension == 1 + len(columns) - projected_rank
