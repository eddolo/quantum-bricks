import statistics
import numpy as np
from .constants import COVALENT_RADII

def build_bonds_inorg(atoms, metric: np.ndarray, scale: float = 1.15):
    """
    Build simple bond list inside the unit cell using covalent radii.
    """
    if not atoms:
        return []

    fracs = np.array([a["frac"] for a in atoms])
    carts = fracs @ metric

    n = len(atoms)
    bonds = []
    for i in range(n):
        for j in range(i + 1, n):
            ri = COVALENT_RADII.get(atoms[i]["elem"], 1.5)
            rj = COVALENT_RADII.get(atoms[j]["elem"], 1.5)
            d = float(np.linalg.norm(carts[i] - carts[j]))
            if 0.4 < d < scale * (ri + rj):
                bonds.append((i, j, d))
    return bonds


def build_bonds_periodic(atoms, metric, scale=1.15):
    """
    FULL PERIODIC BOND BUILDER (3×3×3 SUPERCELL)
    Like Materials Project / pymatgen.
    """
    if not atoms:
        return []

    fracs = np.array([a["frac"] for a in atoms])
    cartes = fracs @ metric
    bonds = []

    shifts = np.array(
        [[i, j, k] for i in (-1, 0, 1)
                   for j in (-1, 0, 1)
                   for k in (-1, 0, 1)]
    )

    n = len(atoms)

    for i in range(n):
        ei = atoms[i]["elem"]
        ri = COVALENT_RADII.get(ei, 1.5)

        for j in range(i + 1, n):
            ej = atoms[j]["elem"]
            rj = COVALENT_RADII.get(ej, 1.5)

            cutoff = scale * (ri + rj)

            base_j = fracs[j]
            for s in shifts:
                f_j = base_j + s
                cart_j = f_j @ metric
                d = float(np.linalg.norm(cartes[i] - cart_j))

                if 0.4 < d < cutoff:
                    bonds.append((i, j, d))
                    break

    return bonds


def inorganic_gap_from_bonds(atoms, bonds):
    """
    Base inorganic band gap estimator (covalent part):

    Eg_base ≈ K / d_nn^2, where d_nn is average nearest-neighbour bond length.
    """
    if not atoms or not bonds:
        return 0.0, 0.0, None

    ds = [d for _, _, d in bonds]
    ds_sorted = sorted(ds)
    half = max(1, len(ds_sorted) // 2)
    d_nn = statistics.mean(ds_sorted[:half])

    K = 6.2  # eV·Å²
    Eg_base = K / (d_nn ** 2)

    disp = 2.0 / d_nn  # heuristic "dispersion-like" proxy

    return Eg_base, disp, d_nn
