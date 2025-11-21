import math
from .family_classifier import normalized_stoichiometry

IONIC_RADIUS_A = {
    "Ca": 1.34,
    "Sr": 1.44,
    "Ba": 1.61,
    "La": 1.36,
    "K":  1.64,
    "Na": 1.39,
    "Pb": 1.49,
}

IONIC_RADIUS_B = {
    "Ti": 0.605,
    "Ta": 0.64,
    "Nb": 0.64,
    "Mn": 0.645,
    "Al": 0.535,
    "Ga": 0.62,
}

IONIC_RADIUS_X_OXIDE = {
    "O": 1.40,
}

IONIC_RADIUS_X_HALIDE = {
    "Cl": 1.81,
    "Br": 1.96,
    "I":  2.20,
}


def compute_perovskite_metrics(atoms, family):
    """
    Compute Goldschmidt tolerance factor t and octahedral factor μ
    for ABO3 (oxides) or ABX3 (halide) perovskites.
    """
    stoich = normalized_stoichiometry(atoms)
    elems_set = set(stoich.keys())

    A_site = {
        "Ba","Sr","Ca","K","Na","Pb","La","Ce","Pr","Nd","Sm","Eu",
        "Gd","Dy","Y","Bi"
    }
    B_site = {
        "Ti","Zr","Hf","V","Nb","Ta","Mn","Fe","Co","Ni","Cr","Sc",
        "Al","Ga"
    }

    A_elems = [e for e in elems_set if e in A_site]
    B_elems = [e for e in elems_set if e in B_site]

    if not A_elems or not B_elems:
        return None, None

    if family == "ABO3_PEROVSKITE":
        X_elems = ["O"] if "O" in elems_set else []
        X_radii_table = IONIC_RADIUS_X_OXIDE
    elif family == "ABX3_HALIDE":
        halides = {"Cl","Br","I"}
        X_elems = [e for e in elems_set if e in halides]
        X_radii_table = IONIC_RADIUS_X_HALIDE
    else:
        return None, None

    if not X_elems:
        return None, None

    def avg_radius(elems, table):
        vals = [table[e] for e in elems if e in table]
        if not vals:
            return None
        return float(sum(vals) / len(vals))

    r_A = avg_radius(A_elems, IONIC_RADIUS_A)
    r_B = avg_radius(B_elems, IONIC_RADIUS_B)
    r_X = avg_radius(X_elems, X_radii_table)

    if r_A is None or r_B is None or r_X is None:
        return None, None

    t = (r_A + r_X) / (math.sqrt(2.0) * (r_B + r_X))
    mu = r_B / r_X

    return float(t), float(mu)
