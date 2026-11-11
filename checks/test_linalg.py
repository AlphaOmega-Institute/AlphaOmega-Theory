"""Regression checks for numerical spectral grouping, independent of a manuscript label."""
import numpy as np

from lib.linalg import eigenspace, response_weights, spectral_blocks


def test_nearby_spectral_blocks_are_disjoint():
    """Resolved nearby eigenvalues must not share eigenvectors or double-count port mass."""
    for offset, gap in ((1.0, 1e-6), (1.0, 4e-9), (1000.0, 1e-6)):
        H = np.diag([offset, offset + gap])
        blocks = spectral_blocks(H)
        assert len(blocks) == 2
        assert [E.shape[1] for _, E in blocks] == [1, 1]
        projectors = [E @ E.conj().T for _, E in blocks]
        assert np.allclose(sum(projectors), np.eye(2))
        assert np.allclose(projectors[0] @ projectors[1], 0)
        for Gamma in (np.eye(2), np.array([[1.0], [1.0]])):
            weights = [R for _, R in response_weights(H, Gamma)]
            assert np.allclose(sum(weights), Gamma.conj().T @ Gamma)
        assert eigenspace(H, offset).shape[1] == 1

    # Numerically coincident eigenvalues still share one block at the declared precision.
    H = np.diag([1.0, 1.0 + 2e-11, 2.0])
    blocks = spectral_blocks(H)
    assert [E.shape[1] for _, E in blocks] == [2, 1]
    assert np.allclose(sum(E @ E.conj().T for _, E in blocks), np.eye(3))
