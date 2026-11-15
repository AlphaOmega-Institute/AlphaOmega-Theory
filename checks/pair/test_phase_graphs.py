"""Checks for phase graphs on the active carrier and for their evolution in the pair model."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import block_diag, expm

from lib.fock import dag
from pair.model import PairModel


LABELS = {
    "UM": ["eq:polar-graph-projector", "prop:phase-readout",
           "eq:graph-evolution-phase", "eq:pair-graph-departure"],
}


def graph_map(orientation, phase):
    return np.vstack((np.eye(orientation.shape[1]), np.exp(1j * phase) * orientation)) / np.sqrt(2)


def cross_observable(upper):
    n_in, n_out = upper.shape
    return np.block([[np.zeros((n_in, n_in)), upper],
                     [dag(upper), np.zeros((n_out, n_out))]])


# Embed a complex input plane in three output coordinates; the third direction
# perpendicular to the image tests the distinction between active and ambient.
ORIENTATION = np.array([[1 / np.sqrt(2), 0], [1j / np.sqrt(2), 0], [0, 1j]])
OUTSIDE = np.array([1, -1j, 0]) / np.sqrt(2)


def test_phase_graph_projectors_and_sector_population():
    v = ORIENTATION
    active = block_diag(np.eye(2), v @ dag(v))
    grading = block_diag(np.eye(2), -np.eye(3))
    reflection = cross_observable(dag(v))
    quadrature = cross_observable(-1j * dag(v))
    psi = np.array([3, 4j]) / 5
    effect = np.array([[0.6, 0.2j], [-0.2j, 0.3]])
    paired = block_diag(effect, v @ effect @ dag(v))

    for phase in (0, 0.72, np.pi / 2, np.pi):
        graph = graph_map(v, phase)
        projection = (active + np.cos(phase) * reflection + np.sin(phase) * quadrature) / 2
        assert_allclose(graph @ dag(graph), projection, atol=1e-12)
        assert_allclose(dag(graph) @ grading @ graph, 0, atol=1e-12)
        assert_allclose(dag(graph) @ graph, np.eye(2), atol=1e-12)
        state = graph @ psi
        for measurement, target in ((reflection, np.cos(phase)), (quadrature, np.sin(phase))):
            assert_allclose(np.vdot(state, measurement @ state), target, atol=1e-12)

        for population in (0, 0.3, 0.5, 1):
            state = np.concatenate((np.sqrt(population) * psi,
                                    np.sqrt(1 - population) * np.exp(1j * phase) * (v @ psi)))
            assert_allclose(np.vdot(state, paired @ state), np.vdot(psi, effect @ psi), atol=1e-12)
            assert_allclose(np.vdot(state, grading @ state), 2 * population - 1, atol=1e-12)


def test_phase_readout_depends_on_active_cross_block_and_input():
    v = ORIENTATION
    first, second = np.eye(2)
    imaginary = 1j * np.outer(second, second) @ dag(v)
    outside = np.outer(second, OUTSIDE.conj())
    assert_allclose(np.linalg.norm(outside), 1)
    psi = (first + 1j * second) / np.sqrt(2)

    for phase in (0, np.pi / 2, np.pi, 3 * np.pi / 2):
        graph = graph_map(v, phase)
        # The ambient cross block is nonzero but vanishes on every phase graph.
        assert_allclose(dag(graph) @ cross_observable(outside) @ graph, 0, atol=1e-12)
        # The active cross block is visible for psi and silent for first.
        compressed = dag(graph) @ cross_observable(imaginary) @ graph
        assert_allclose(np.vdot(first, compressed @ first), 0, atol=1e-12)
        assert_allclose(np.vdot(psi, compressed @ psi), -np.sin(phase) / 2, atol=1e-12)


def test_pair_phase_graph_evolution_and_departure():
    model = PairModel(theta=np.pi / 8, j=1.3)
    original_pairs = np.column_stack((model.C, model.D))
    v = dag(model.F2) @ model.A
    phase = 0.73
    t = 0.81
    minimal = model.exchange_hamiltonian(model.j)
    restricted = dag(original_pairs) @ minimal @ original_pairs
    output_evolution = expm(-1j * t * (dag(model.F2) @ minimal @ model.F2))

    for shift in (0.0, 0.2 * model.j):
        input_evolution = expm(-1j * t * (restricted + shift * np.eye(2)))
        evolved_graph = block_diag(input_evolution, output_evolution) @ graph_map(v, phase)
        assert_allclose(evolved_graph, graph_map(v, phase + shift * t) @ input_evolution, atol=1e-12)

    transfer_time = np.pi / model.j
    input_evolution = expm(-1j * transfer_time * (dag(original_pairs) @ model.H @ original_pairs))
    output_evolution = expm(-1j * transfer_time * (dag(model.F2) @ model.H @ model.F2))
    evolved = block_diag(input_evolution, output_evolution) @ graph_map(v, phase) @ model.ePsi
    # Output is orthogonal to the old final support, so its term in every
    # original graph projection is zero. The remaining input term has norm 1/2.
    assert_allclose(dag(v) @ evolved[2:], 0, atol=1e-12)
    for target_phase in (0, phase, np.pi / 2, np.pi):
        overlap = dag(graph_map(v, target_phase)) @ evolved
        assert_allclose(np.vdot(overlap, overlap), 0.25, atol=1e-12)
