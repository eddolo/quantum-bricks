from collections import Counter
from math import gcd

# ----------------------------
# Extended family maps
# ----------------------------

III_V = {
    ("Ga","As"):  "GAAS",
    ("In","P"):   "INP",
    ("In","As"):  "INAS",
    ("Al","As"):  "ALAS",
    ("Ga","P"):   "GAP",
}

II_VI_EXT = {
    ("Zn","S"):   "ZNS",
    ("Zn","Se"):  "ZNSE",
    ("Zn","Te"):  "ZNTE",
    ("Cd","Te"):  "CDTE",
}

OXIDE_SPECIAL = {
    ("Fe","O"):   "FEO",
    ("Co","O"):   "COO",
    ("Ni","O"):   "NIO",
    ("Cu","O"):   "CUO",
    ("Mn","O"):   "MNO",
}

CHALCOGENIDES = {
    ("Sn","S"): "SNS",
    ("Sn","Se"): "SNSE",
    ("Sn","Te"): "SNTE",
    ("Pb","S"): "PBS",
    ("Pb","Se"): "PBSE",
    ("Pb","Te"): "PBTE",
}

HALIDE_SALTS = {
    ("Na","Cl"): "NACL",
    ("K","Cl"):  "KCL",
    ("Cs","Cl"): "CSCL",
    ("Li","F"):  "LIF",
    ("Na","F"):  "NAF",
}

# ----------------------------
# Helpers
# ----------------------------

def normalized_stoichiometry(atoms):
    counts = Counter(a["elem"] for a in atoms)
    if not counts:
        return {}
    g = 0
    for c in counts.values():
        g = c if g == 0 else gcd(g, c)
    if g <= 1:
        return dict(counts)
    return {e: c // g for e, c in counts.items()}


def formula_from_atoms(atoms):
    sto = normalized_stoichiometry(atoms)
    order = sorted(sto.items(), key=lambda x: x[0])
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in order)


# ----------------------------
# Main classification
# ----------------------------

def classify_family(atoms):
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

    # --- ABO3 perovskites ---
    if "O" in elems_set:
        if any(e in A_site for e in elems_set) and any(e in B_site for e in elems_set):
            return "ABO3_PEROVSKITE"

    # --- ABX3 halide perovskites ---
    if any(e in {"Cl","Br","I"} for e in elems_set):
        A_found = any(e in A_site or e in {"Cs","Rb","Tl"} for e in elems_set)
        B_found = any(e in {"Pb","Sn"} for e in elems_set)
        if A_found and B_found:
            return "ABX3_HALIDE"

    # --- Elemental ---
    if elems_set == {"C"}:
        return "DIAMOND"
    if elems_set == {"Si"}:
        return "SI"
    if elems_set == {"Ge"}:
        return "GE"

    # --- Simple binaries already in engine ---
    if elems_set == {"Zn","O"}:
        return "ZNO"
    if elems_set == {"Ti","O"}:
        return "TIO2"
    if elems_set == {"Mg","O"}:
        return "MGO"

    if elems_set == {"Al","O"}:
        a = stoich.get("Al", 0)
        o = stoich.get("O", 0)
        if 3 * a == 2 * o:
            return "AL2O3"

    if elems_set == {"Al","N"}:
        return "ALN"
    if elems_set == {"Ga","N"}:
        return "GAN"

    if elems_set == {"Cd","S"}:
        return "CDS"
    if elems_set == {"Cd","Se"}:
        return "CDSE"

    if elems_set == {"Mo","S"} or elems_set == {"Mo","S","H"}:
        return "MOS2"

    # --- Ruddlesden–Popper n=1 ---
    if stoich.get("O", 0) == 4:
        if any(e in A_site for e in elems_set) and any(e in B_site for e in elems_set):
            return "RP_LAYERED"

    # ----------------------------------------
    # NEW FAMILIES
    # ----------------------------------------

    # III–V
    for (a, b), fam in III_V.items():
        if elems_set == {a, b}:
            return fam

    # II–VI extended
    for (a, b), fam in II_VI_EXT.items():
        if elems_set == {a, b}:
            return fam

    # Transition-metal monoxides
    for (a, b), fam in OXIDE_SPECIAL.items():
        if elems_set == {a, b}:
            return fam

    # SnX/PbX chalcogenides
    for (a, b), fam in CHALCOGENIDES.items():
        if elems_set == {a, b}:
            return fam

    # Simple halide salts
    for (a, b), fam in HALIDE_SALTS.items():
        if elems_set == {a, b}:
            return fam

    return "GENERIC"
