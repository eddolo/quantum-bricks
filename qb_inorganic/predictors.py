def canonical_formula(stoich: dict, family: str) -> str:
    """Return canonical chemical formula string."""

    # Special perovskite canonicalization
    if family in {"ABO3_PEROVSKITE", "ABX3_HALIDE"}:
        A = None
        B = None
        X = []

        for elem, count in stoich.items():
            if elem in {"O", "F", "Cl", "Br", "I"}:
                X.append((elem, count))
            else:
                if A is None:
                    A = (elem, count)
                else:
                    B = (elem, count)

        X_str = "".join(f"{e}{(c if c > 1 else '')}" for e, c in X)
        return f"{A[0]}{B[0]}{X_str}"

    metals = []
    nonmetals = []

    for elem, count in stoich.items():
        if elem in {"C","H","N","O","F","Cl","Br","I","S","Se","P"}:
            nonmetals.append((elem, count))
        else:
            metals.append((elem, count))

    metals.sort()
    nonmetals.sort()

    items = metals + nonmetals
    return "".join(f"{e}{c if c > 1 else ''}" for e, c in items)


def estimate_fermi_level(Eg: float, family: str, ionicity: float):
    """Return relative Fermi offset from mid-gap."""
    if Eg <= 0:
        return 0.0

    Ef = 0.0

    if family in {"ZNO","GAN","ALN","CDS"}:
        Ef += 0.10 * (ionicity if ionicity is not None else 0.0)

    if family in {"TIO2","MGO","AL2O3"}:
        Ef -= 0.05 * ((ionicity if ionicity is not None else 0.0) + 1.0)

    if "PEROVSKITE" in family:
        Ef += 0.05 * ((ionicity if ionicity is not None else 0.0) - 1.0)

    Ef = max(-0.25, min(0.25, Ef))
    return Ef


def estimate_dielectric(Eg: float, family: str):
    """Return (eps_inf, eps_static)."""
    if Eg <= 0:
        return (1.0, 1.0)

    try:
        n = (108.0 / Eg) ** 0.25
    except Exception:
        n = 1.5

    eps_inf = n**2

    if family in {"ABO3_PEROVSKITE"}:
        eps_static = eps_inf * 4.0
    elif family in {"ZNO","GAN"}:
        eps_static = eps_inf * 2.0
    else:
        eps_static = eps_inf * 1.5

    return (eps_inf, eps_static)


def estimate_effective_masses(Eg, family, coord_type, angle_disorder, ionicity):
    """Return (m_e, m_h) approximate effective masses (in units of m0)."""

    base = max(0.1, min(2.0, 0.25 + Eg * 0.15))

    if coord_type == "octahedral-like":
        base *= 1.3
    if coord_type == "tetrahedral-like":
        base *= 0.9

    if angle_disorder is not None:
        base *= 1.0 + (angle_disorder / 120.0)

    ion = ionicity if ionicity is not None else 0.0

    m_e = base * (1.0 - 0.2 * ion)
    m_h = base * (1.0 + 0.3 * ion)

    m_e = max(0.05, min(5.0, m_e))
    m_h = max(0.05, min(5.0, m_h))

    return m_e, m_h
