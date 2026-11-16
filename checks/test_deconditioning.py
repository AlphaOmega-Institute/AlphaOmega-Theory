"""Two-moment deconditioning on the unequal-cell path."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import null_space, polar

LABELS = {
    "LG": ["thm:canonical-deconditioning-semigroup", "thm:deconditioning-contracts-excess"],
    "UM": ["prop:two-moment-deconditioning", "prop:ground-polar-compatibility"],
}


def path_pairing():
    copying_map = np.repeat(np.eye(2), [1, 4], axis=0) / 2
    orientation, _ = polar(copying_map)
    return orientation


def damping_operators(target, retention):
    ground = np.outer(target, target.conj())
    no_jump = ground + np.sqrt(retention) * (np.eye(len(target)) - ground)
    excited_basis = null_space(target.conj()[None, :])
    jumps = [np.sqrt(1 - retention) * np.outer(target, state.conj())
             for state in excited_basis.T]
    return [no_jump, *jumps]


def evolve(operators, state):
    return sum(op @ state @ op.conj().T for op in operators)


def test_path_deconditioning_selects_the_ground_partner():
    pairing = path_pairing()
    fine_ground = np.ones(5) / np.sqrt(5)
    coarse_ground = np.ones(2) / np.sqrt(2)
    partner = pairing.T @ fine_ground
    partner_state = np.outer(partner, partner)
    assert_allclose(abs(coarse_ground @ partner)**2, 9 / 10)

    fine_step = damping_operators(fine_ground, 0.5)
    partner_step = damping_operators(partner, 0.5)
    # Matrix units span all coarse density matrices, including both coherences.
    for left in np.eye(2):
        for right in np.eye(2):
            unit = np.outer(left, right)
            fine_output = evolve(fine_step, pairing @ unit @ pairing.T)
            assert_allclose(fine_output, pairing @ evolve(partner_step, unit) @ pairing.T,
                            atol=1e-14)

    assert_allclose(evolve(partner_step, partner_state), partner_state)
    coarse_ground_step = damping_operators(coarse_ground, 0.5)
    assert not np.allclose(evolve(coarse_ground_step, partner_state), partner_state)
    excitation = np.array([2, -1]) / np.sqrt(5)
    relaxed = evolve(partner_step, np.outer(excitation, excitation))
    assert_allclose(np.linalg.eigvalsh(relaxed), [0.5, 0.5])


def test_target_position_controls_excess_weight():
    pairing = path_pairing()
    excess = np.eye(5) - pairing @ pairing.T
    fine_ground = np.ones(5) / np.sqrt(5)
    outside = np.array([0, 1, 1, -1, -1]) / 2

    # Pull back the effect itself so the check includes all cross-space coherences.
    for retention in (0.25, 0.5, 0.8):
        operators = damping_operators(fine_ground, retention)
        effect = sum(op.conj().T @ excess @ op for op in operators)
        assert_allclose(effect, retention * excess, atol=1e-14)

    # The top eigenvalue on the paired space maximizes leakage over all paired inputs.
    oblique_target = np.sqrt(2 / 3) * fine_ground + 1j * outside / np.sqrt(3)
    operators = damping_operators(oblique_target, 0.5)
    effect = sum(op.conj().T @ excess @ op for op in operators)
    paired_effect = pairing.T @ effect @ pairing
    assert_allclose(np.linalg.eigvalsh(paired_effect)[-1], 1 / 6)
