import math
import statistics
import numpy as np

from .constants import (
    ELECTRONEGATIVITY,
    TYPICAL_VALENCE,
    ANION_LIKE,
    CATION_LIKE,
)

def analyze_local_env(atoms, bonds, metric: np.ndarray):
    """
    Analyze cation–anion bonds, coordination, and angle disorder.

    Returns:
        avg_delta_ch, avg_valence, coord_type, rms_angle_deg
    """
    if not atoms or not bonds:
        return None, None, None, None

    fracs = np.array([a["frac"] for a in atoms])
    carts = fracs @ metric
    n = len(atoms)

    elems = [a["elem"] for a in atoms]

    adj = [[] for _ in range(n)]
    for i, j, d in bonds:
        adj[i].append((j, d))
        adj[j].append((i, d))

    delta_ch_vals = []
    valences = []
    cation_centers = []

    for i, j, d in bonds:
        ei = elems[i]
        ej = elems[j]

        chi_i = ELECTRONEGATIVITY.get(ei)
        chi_j = ELECTRONEGATIVITY.get(ej)
        if chi_i is None or chi_j is None:
            continue

        if ei in CATION_LIKE and ej in ANION_LIKE:
            cation = i
            anion = j
        elif ej in CATION_LIKE and ei in ANION_LIKE:
            cation = j
            anion = i
        else:
            continue

        delta = abs(chi_i - chi_j)
        delta_ch_vals.append(delta)

        val = TYPICAL_VALENCE.get(elems[cation])
        if val is not None:
            valences.append(val)

        if cation not in cation_centers:
            cation_centers.append(cation)

    if not delta_ch_vals:
        return None, None, None, None

    avg_delta_ch = float(statistics.mean(delta_ch_vals)) if delta_ch_vals else None
    avg_valence = float(statistics.mean(valences)) if valences else None

    coord_counts = []
    angle_devs = []

    for c in cation_centers:
        neighs = [j for (j, d) in adj[c] if elems[j] in ANION_LIKE]
        if len(neighs) < 2:
            continue

        coord_counts.append(len(neighs))

        rc = carts[c]
        vecs = [carts[j] - rc for j in neighs]

        local_angles = []
        for i in range(len(vecs)):
            for j in range(i + 1, len(vecs)):
                v1 = vecs[i]
                v2 = vecs[j]
                n1 = np.linalg.norm(v1)
                n2 = np.linalg.norm(v2)
                if n1 < 1e-6 or n2 < 1e-6:
                    continue
                cosang = np.dot(v1, v2) / (n1 * n2)
                cosang = max(-1.0, min(1.0, cosang))
                ang = math.degrees(math.acos(cosang))
                local_angles.append(ang)

        if not local_angles:
            continue

        if len(neighs) <= 4:
            ideal = 109.5
        else:
            ideal = 90.0

        devs = [(ang - ideal) for ang in local_angles]
        for d in devs:
            angle_devs.append(d)

    if coord_counts:
        avg_coord = statistics.mean(coord_counts)
        if avg_coord < 4.5:
            coord_type = "tetrahedral-like"
        elif avg_coord < 5.5:
            coord_type = "mixed"
        else:
            coord_type = "octahedral-like"
    else:
        coord_type = None

    if angle_devs:
        rms_angle = math.sqrt(statistics.mean([d * d for d in angle_devs]))
    else:
        rms_angle = None

    return avg_delta_ch, avg_valence, coord_type, rms_angle
