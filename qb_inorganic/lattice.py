import re
import numpy as np

def parse_cell(text: str) -> np.ndarray:
    """Extract a,b,c,alpha,beta,gamma from CIF text."""
    cell = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    keys = [
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",
    ]
    for i, k in enumerate(keys):
        m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", text)
        if m:
            cell[i] = float(m.group(1))
    return np.array(cell, float)


def build_metric(cell: np.ndarray) -> np.ndarray:
    """Build 3×3 lattice matrix from cell parameters."""
    a, b, c, al, be, ga = cell
    al, be, ga = np.radians([al, be, ga])

    v = np.sqrt(
        1
        - np.cos(al) ** 2
        - np.cos(be) ** 2
        - np.cos(ga) ** 2
        + 2 * np.cos(al) * np.cos(be) * np.cos(ga)
    )

    return np.array(
        [
            [a, b * np.cos(ga), c * np.cos(be)],
            [0.0, b * np.sin(ga), c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)],
            [0.0, 0.0, c * v / np.sin(ga)],
        ]
    )
