"""Random-model tests of the three shell estimates and of the constants in them (hbar = 1).

The committed checks of UM prop:shell-accuracy use two-level systems and one truncated
trap.  Here the hypotheses are sampled directly: H0 has a degenerate lowest level of
dimension six at energy zero and every other level at or above the gap; V is a random
Hermitian matrix whose spectrum lies in an interval of width nu placed anywhere, with nu
ranging from small values to 1.5 gaps.  Seeds are fixed, so the run is reproducible.
"""
import numpy as np

LABELS = {
    "UM": ["prop:shell-accuracy", "eq:shell-ritz", "eq:shell-leak", "eq:shell-drift"],
}

GAP = 1.7
SHELL, OTHER = 6, 10


def _model(rng, width):
    n = SHELL + OTHER
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    unitary, _ = np.linalg.qr(raw)
    position = rng.beta(0.4, 0.4, size=n)                      # pushed toward both ends of the interval
    offset = width * rng.uniform(-1, 1)
    interaction = (unitary * (offset - width * position)) @ unitary.conj().T
    upper = GAP * (1 + rng.exponential(0.3, size=OTHER) * (rng.random(OTHER) < 0.5))
    return np.diag(np.r_[np.zeros(SHELL), upper]).astype(complex), interaction


def test_three_estimates_on_random_models():
    rng = np.random.default_rng(1618)
    closest = {"energy": 0.0, "leakage": 0.0, "drift": 0.0}
    for _ in range(600):
        width = rng.uniform(0.02, 1.5) * GAP
        h0, v = _model(rng, width)
        h = h0 + v
        ritz, ritz_vectors = np.linalg.eigh(h[:SHELL, :SHELL])
        exact, vectors = np.linalg.eigh(h)
        errors = ritz - exact[:SHELL]
        assert errors.min() >= -1e-12
        assert errors.max() <= width**2 / (4 * GAP) + 1e-12
        closest["energy"] = max(closest["energy"], errors.max() / (width**2 / (4 * GAP)))
        for t in rng.uniform(0, 200, size=12):
            evolved = (vectors * np.exp(-1j * exact * t)) @ vectors.conj().T[:, :SHELL]
            leakage = np.linalg.norm(evolved[SHELL:], 2) ** 2
            assert leakage <= min(1.0, (width / GAP) ** 2) + 1e-12
            closest["leakage"] = max(closest["leakage"], leakage / min(1.0, (width / GAP) ** 2))
            shell_evolution = (ritz_vectors * np.exp(-1j * ritz * t)) @ ritz_vectors.conj().T
            drift = np.linalg.norm(evolved[:SHELL] - shell_evolution, 2)
            assert drift <= width**2 * t / (2 * GAP) + 1e-12
            closest["drift"] = max(closest["drift"], drift / (width**2 * t / (2 * GAP)))
    # The samples are not vacuous: each estimate is approached to within a factor of four.
    assert closest["energy"] > 0.5 and closest["leakage"] > 0.5 and closest["drift"] > 0.25


def test_energy_bound_holds_beyond_the_gap():
    rng = np.random.default_rng(577)
    for _ in range(300):
        width = rng.uniform(1.0, 3.0) * GAP
        h0, v = _model(rng, width)
        h = h0 + v
        errors = np.linalg.eigvalsh(h[:SHELL, :SHELL]) - np.linalg.eigvalsh(h)[:SHELL]
        assert errors.min() >= -1e-12 and errors.max() <= width**2 / (4 * GAP) + 1e-12


def test_both_constants_are_attained_by_two_level_reflections():
    for ratio in (0.05, 0.3, 0.9):
        width = ratio * GAP
        for cosine, kind in ((width / (2 * GAP), "energy"), (width / GAP, "leakage")):
            sine = np.sqrt(1 - cosine**2)
            v = (width / 2) * np.array([[cosine, sine], [sine, -cosine]])
            h = np.diag([0.0, GAP]) + v
            values, vectors = np.linalg.eigh(h)
            if kind == "energy":
                assert np.isclose(h[0, 0] - values[0], width**2 / (4 * GAP), atol=1e-13)
            else:
                gap = values[1] - values[0]
                amplitude = (vectors * np.exp(-1j * values * np.pi / gap)) @ vectors.conj().T[:, 0]
                assert np.isclose(abs(amplitude[1]) ** 2, (width / GAP) ** 2, atol=1e-12)


def test_one_dimensional_shell_gives_an_isolated_branch():
    rng = np.random.default_rng(5)
    for _ in range(800):
        n = 9
        width = rng.uniform(0.02, 0.98) * GAP
        raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        unitary, _ = np.linalg.qr(raw)
        v = -(unitary * (width * rng.beta(0.4, 0.4, size=n))) @ unitary.conj().T      # attraction: -width <= V <= 0
        upper = GAP * (1 + rng.exponential(0.3, size=n - 1) * (rng.random(n - 1) < 0.5))
        values = np.linalg.eigvalsh(np.diag(np.r_[0.0, upper]) + v)
        assert values[0] <= 1e-12
        assert values[1] >= GAP - width - 1e-12
