"""Check each method against an independent reference before trusting it.

This is not a unit-test suite; it is the short verification a careful reader would
want: the effective-resistance machinery agrees with networkx, the energy distance
matches its definition, and a diagram's distance to itself is zero.

    python experiments/verify_methods.py
"""
import networkx as nx
import numpy as np
from scipy.spatial.distance import cdist, pdist

import config  # noqa: F401  (puts the package on the path)
from topoperturb import distances, homology, twosample


def check_effective_resistance():
    """The Laplacian-pseudoinverse resistance equals networkx's resistance_distance."""
    rng = np.random.default_rng(0)
    A = distances.knn_graph(rng.normal(size=(40, 5)))
    graph = nx.from_numpy_array(A)
    deg = A.sum(1)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        Lp = np.linalg.pinv(np.diag(deg) - A)
        diag = np.diag(Lp)
        R = diag[:, None] + diag[None, :] - 2 * Lp      # resistance before the paper's correction
        error = max(abs(R[i, j] - nx.resistance_distance(graph, i, j))
                    for i, j in [(0, 1), (2, 7), (5, 20)])
    assert error < 1e-6, error
    print(f"effective resistance agrees with networkx (max error {error:.1e})")


def check_energy_distance():
    """Energy distance matches 2E|A-B| - E|A-A'| - E|B-B'| and is non-negative."""
    rng = np.random.default_rng(1)
    A, B = rng.normal(size=(60, 4)), rng.normal(size=(60, 4)) + 0.5
    reference = 2 * cdist(A, B).mean() - pdist(A).mean() - pdist(B).mean()
    assert abs(twosample.energy_distance(A, B) - reference) < 1e-9
    assert twosample.energy_distance(A, B) > 0
    print("energy distance matches its definition and is positive for differing clouds")


def check_diagram_distance():
    """A persistence diagram is zero distance from itself (ripser + persim)."""
    rng = np.random.default_rng(2)
    diagram = homology.summarize(rng.normal(size=(30, 3)), ("euclidean",))["euclidean"]["H1"]
    assert homology.diagram_distance(diagram, diagram) == 0.0
    print("diagram distance is zero between identical diagrams")


def check_permutation_calibration():
    """The permutation test on identical clouds rarely rejects (sanity on the null)."""
    rng = np.random.default_rng(3)
    energy_only = lambda A, B: {"energy": twosample.energy_distance(A, B)}
    rejects = 0
    for _ in range(40):
        A, B = rng.normal(size=(60, 5)), rng.normal(size=(60, 5))
        p, _ = twosample.permutation_test(energy_only, A, B, rng=rng)["energy"]
        rejects += p < 0.05
    assert rejects <= 6, rejects         # ~5% expected, allow slack on 40 draws
    print(f"permutation test on the null rejected {rejects}/40 times (≈5% expected)")


if __name__ == "__main__":
    check_effective_resistance()
    check_energy_distance()
    check_diagram_distance()
    check_permutation_calibration()
    print("all method checks passed")
