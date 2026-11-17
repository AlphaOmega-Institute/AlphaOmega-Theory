"""Finite blocks of the grounded shift family and finite input extensions."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm

from lib.linalg import rank

LABELS = {
    "UM": ["eq:index-shift-family", "eq:index-shift-energy",
           "eq:index-channel-coherence", "eq:index-extension-counts"],
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


def test_excitation_energy_selects_a_common_shift_defect():
    g0, g1, g2 = np.eye(3)
    shift = np.column_stack((g1, g2))
    input_energies = np.array([0.0, 1.0])
    h_input = np.diag(input_energies)
    angle = 0.37
    rotated = np.column_stack((g1, np.cos(angle) * g2 - np.sin(angle) * g0))
    missing_line = np.outer(g0, g0)
    for energy, expected_defect in (
        (1.0, np.zeros((3, 3))),
        (0.5, missing_line),
        (7.3, missing_line),
    ):
        output_energies = np.array([energy, 0.0, 1.0])
        h_output = np.diag(output_energies)
        assert_allclose(h_output @ shift, shift @ h_input)
        residual = h_output @ rotated - rotated @ h_input
        assert_allclose(np.linalg.norm(residual, 2), abs((energy - 1) * np.sin(angle)))

        # Matrix unit |g_i><f_j| intertwines precisely at equal energies.
        # Row i is in the common defect when every such entry vanishes.
        matching_energies = output_energies[:, None] == input_energies[None, :]
        common_defect = np.diag(~matching_energies.any(axis=1))
        assert_allclose(common_defect, expected_defect)


def test_covariant_channel_saturates_two_factor_coherence_bound():
    g0, g1, g2 = np.eye(3)
    f0, f1 = np.eye(2)
    h_input = np.diag([0, 1])
    epsilon = 0.5
    h_output = np.diag([epsilon, 0, 1])
    input_coherence = np.outer(f0, f1)
    for a in (0, 0.2, 0.5, 0.9, 1):
        for p in (0, 0.1, 0.3, 0.8, 1):
            keep = np.column_stack((np.sqrt(1 - a) * g1, np.sqrt(1 - p) * g2))
            from_ground = np.sqrt(a) * np.outer(g0, f0)
            from_excitation = np.sqrt(p) * np.outer(g0, f1)
            kraus = (keep, from_ground, from_excitation)
            assert_allclose(sum(k.T @ k for k in kraus), np.eye(2))
            for k, shift in zip(kraus, (0, epsilon, epsilon - 1)):
                assert_allclose(h_output @ k - k @ h_input, shift * k)
            for state, weight in ((f0, a), (f1, p)):
                assert_allclose(sum(abs(g0 @ k @ state) ** 2 for k in kraus), weight)
            coherence = sum(k @ input_coherence @ k.T for k in kraus)
            bound = np.sqrt((1 - a) * (1 - p))
            assert_allclose(coherence, bound * np.outer(g1, g2))
            assert_allclose(np.linalg.norm(coherence, ord="nuc"), bound)


def test_line_coherence_requires_resonance_and_can_require_moving_the_ground():
    g0, g1, g2 = np.eye(3)
    f0, f1 = np.eye(2)
    input_energies = np.array([0, 1])
    # Use a positive amplitude for every allowed transition. Terms from
    # different Kraus operators then cannot cancel on a matrix unit.
    for epsilon, preserve_ground, expect_line_coherence in (
        (0.5, False, False), (0.5, True, False),
        (2, False, True), (2, True, False), (1, True, True),
    ):
        ground_prob = [0, 1, 0] if preserve_ground else [0.2, 0.3, 0.5]
        amplitudes = np.sqrt(np.column_stack((ground_prob, [0.4, 0.2, 0.4])))
        shifts = np.array([epsilon, 0, 1])[:, None] - input_energies[None, :]
        kraus = [amplitudes * (shifts == shift) for shift in np.unique(shifts)]
        assert_allclose(sum(k.T @ k for k in kraus), np.eye(2))
        line_coherences = []
        for left in (f0, f1):
            for right in (f0, f1):
                out = sum(np.outer(k @ left, k @ right) for k in kraus)
                line_coherences.extend(out[0, 1:])
        assert bool(np.any(line_coherences)) == expect_line_coherence
        ordinary_coherence = sum((g1 @ k @ f0) * (g2 @ k @ f1) for k in kraus)
        assert ordinary_coherence > 0

    # At epsilon=2 an isometry shifting both energies by +1 transfers all
    # the original coherence. It takes the input ground to an excitation.
    resonant = np.column_stack((g2, g0))
    assert_allclose(resonant.T @ resonant, np.eye(2))
    assert_allclose(np.diag([2, 0, 1]) @ resonant - resonant @ np.diag(input_energies), resonant)
    assert_allclose(resonant @ np.outer(f0, f1) @ resonant.T, np.outer(g2, g0))


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
