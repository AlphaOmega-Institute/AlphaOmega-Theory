"""Regression checks for numerical spectral grouping, independent of a manuscript label."""
import numpy as np

from lib.linalg import eigenspace, minimal_space, rank, response_weights, spectral_blocks


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


def _disjoint_and_exhaustive(blocks, dim):
    projectors = [E @ E.conj().T for _, E in blocks]
    assert np.allclose(sum(projectors), np.eye(dim))
    for i in range(len(projectors)):
        for k in range(i + 1, len(projectors)):
            assert np.allclose(projectors[i] @ projectors[k], 0)


def test_unresolved_pair_across_a_rounding_boundary_shares_a_block():
    """Two eigenvalues 2e-13 apart are unresolved at the declared tolerance wherever they
    sit on the number line; in particular not only when they round to the same decimal."""
    H = np.diag([0.5000000005 - 1e-13, 0.5000000005 + 1e-13, 2.0])
    Gamma = np.array([[1.0], [1.0], [0.0]]) / np.sqrt(2)
    blocks = spectral_blocks(H)
    assert [E.shape[1] for _, E in blocks] == [2, 1]
    _disjoint_and_exhaustive(blocks, 3)
    weights = [R for _, R in response_weights(H, Gamma)]
    assert np.allclose(sum(weights), Gamma.conj().T @ Gamma)
    assert sum(rank(R) for R in weights) == 1
    assert minimal_space(H, Gamma).shape[1] == 1


def test_groups_never_span_more_than_the_tolerance():
    """Neighbouring eigenvalues each within the tolerance must not chain into a wide group:
    every group's total spread stays within the declared tolerance."""
    tol = 1e-9
    for values in ([0.0, 0.7e-9, 1.4e-9], list(0.6e-9 * np.arange(10)), [3.0, 3.0 + 0.9e-9, 3.0 + 1.8e-9, 5.0]):
        H = np.diag(values)
        blocks = spectral_blocks(H, tol=tol)
        _disjoint_and_exhaustive(blocks, len(values))
        w = np.array(values)
        for _, E in blocks:
            members = w[np.abs(E).sum(axis=1) > 0.5]
            assert members.max() - members.min() <= tol
        assert len(blocks) >= 2
    assert [E.shape[1] for _, E in spectral_blocks(np.diag([0.0, 0.7e-9, 1.4e-9]))] == [2, 1]


def test_blocks_partition_random_clustered_spectra():
    """Disjoint, exhaustive blocks and weights that sum to Gamma^* Gamma for random
    Hermitian matrices with tight eigenvalue clusters."""
    rng = np.random.default_rng(3)
    for _ in range(20):
        centres = rng.uniform(-2, 2, size=4)
        values = np.concatenate([c + rng.uniform(-3e-10, 3e-10, size=rng.integers(1, 4)) for c in centres])
        Q, _ = np.linalg.qr(rng.normal(size=(len(values),) * 2) + 1j * rng.normal(size=(len(values),) * 2))
        H = Q @ np.diag(values) @ Q.conj().T
        H = (H + H.conj().T) / 2
        Gamma = rng.normal(size=(len(values), 2)) + 1j * rng.normal(size=(len(values), 2))
        _disjoint_and_exhaustive(spectral_blocks(H), len(values))
        weights = [R for _, R in response_weights(H, Gamma)]
        assert np.allclose(sum(weights), Gamma.conj().T @ Gamma)


def test_merging_distinct_energies_is_an_approximation_by_design():
    """Two genuinely distinct energies closer than the tolerance are merged.  The exact
    response of a port coupled to both needs two directions; the grouped count reports one.
    Lowering the tolerance below their separation resolves them."""
    H = np.diag([1.0, 1.0 + 5e-10])
    Gamma = np.array([[1.0], [1.0]]) / np.sqrt(2)
    assert minimal_space(H, Gamma).shape[1] == 1
    assert [E.shape[1] for _, E in spectral_blocks(H, tol=1e-10)] == [1, 1]
