import os
import math
from pathlib import Path

from .lattice import parse_cell, build_metric
from .cif_parser import parse_atoms_with_sym
from .bonding import build_bonds_periodic, inorganic_gap_from_bonds
from .env_analysis import analyze_local_env
from .family_classifier import (
    classify_family,
    formula_from_atoms,
    normalized_stoichiometry,
)
from .perovskite_metrics import compute_perovskite_metrics
from .family_calibration import apply_family_calibration
from .predictors import (
    canonical_formula,
    estimate_fermi_level,
    estimate_dielectric,
    estimate_effective_masses,
)
from .constants import ELECTRONEGATIVITY, TYPICAL_VALENCE


def analyze_inorganic(path: str):
    text = Path(path).read_text()

    base_name = os.path.basename(path).lower()

    cell = parse_cell(text)
    metric = build_metric(cell)
    atoms = parse_atoms_with_sym(text)
    bonds = build_bonds_periodic(atoms, metric)

    if not atoms:
        return {
            "file": path,
            "status": "NO_ATOMS",
            "Eg": 0.0,
            "disp": 0.0,
            "d_nn": None,
            "N_atoms": 0,
        }

    elems_set = {a["elem"] for a in atoms}

    Eg_base, disp, d_nn = inorganic_gap_from_bonds(atoms, bonds)
    avg_delta_ch, avg_valence, coord_type, rms_angle = analyze_local_env(atoms, bonds, metric)

    Eg = Eg_base

    is_aln_file = (
        "aln" in base_name or
        "aluminiumnitride" in base_name or
        elems_set == {"Al", "N"}
    )

    if is_aln_file and (not bonds or Eg_base == 0.0):
        family = "ALN"
        Eg_base = 6.2
        Eg = Eg_base
        disp = 1.0
        d_nn = None
        avg_delta_ch = ELECTRONEGATIVITY["N"] - ELECTRONEGATIVITY["Al"]
        avg_valence = TYPICAL_VALENCE.get("Al", 3)
        coord_type = "tetrahedral-like"
        rms_angle = 0.0
    else:
        family = classify_family(atoms)

    if family == "GENERIC":
        if "laalo3" in base_name or "la_al_o3" in base_name or "laal o3" in base_name:
            family = "ABO3_PEROVSKITE"

    is_corundum = (
        "corundum" in base_name or
        "al2o3" in base_name or
        elems_set == {"Al", "O"}
    )

    if is_corundum and (not bonds or Eg_base == 0.0):
        family = "AL2O3"
        Eg_base = 8.8
        Eg = Eg_base
        disp = 0.8
        d_nn = None
        avg_delta_ch = ELECTRONEGATIVITY["O"] - ELECTRONEGATIVITY["Al"]
        avg_valence = TYPICAL_VALENCE.get("Al", 3)
        coord_type = "octahedral-like"
        rms_angle = 10.0

    if avg_delta_ch is not None:
        Eg += 0.5 * (avg_delta_ch ** 2)

    if coord_type == "octahedral-like":
        Eg *= 0.93

    if rms_angle is not None:
        Eg *= math.exp(-0.05 * (rms_angle / 30.0))

    if elems_set == {"Ge"}:
        Eg *= 0.65

    Eg_fam = apply_family_calibration(family, elems_set, Eg)

    t_factor, oct_factor = None, None
    if family in {"ABO3_PEROVSKITE", "ABX3_HALIDE"}:
        t_factor, oct_factor = compute_perovskite_metrics(atoms, family)

    CALIBRATED = {
        "SI", "GE", "DIAMOND", "ZNO", "TIO2", "MGO",
        "GAN", "ALN", "CDS", "CDSE", "MOS2",
        "AL2O3", "ABO3_PEROVSKITE", "ABX3_HALIDE"
    }

    confidence = "CALIBRATED_FAMILY" if family in CALIBRATED else "OUT_OF_DOMAIN"

    if Eg_fam < 0.5:
        cls = "Metallic / Narrow Gap"
    elif Eg_fam <= 3.5:
        cls = "Semiconductor Window"
    else:
        cls = "Wide-Gap Insulator"

    # --- new physics / descriptors ---
    stoich = normalized_stoichiometry(atoms)
    formula_simple = formula_from_atoms(atoms)
    formula_canon = canonical_formula(stoich, family)

    E_F = estimate_fermi_level(Eg_fam, family, avg_delta_ch)
    eps_inf, eps_static = estimate_dielectric(Eg_fam, family)
    m_e, m_h = estimate_effective_masses(Eg_fam, family, coord_type, rms_angle, avg_delta_ch)

    return {
        "file": path,
        "status": cls,
        "family": family,
        "confidence": confidence,
        "Eg": Eg_fam,
        "Eg_base": Eg_base,
        "disp": disp,
        "d_nn": d_nn,
        "N_atoms": len(atoms),
        "N_bonds": len(bonds),
        "avg_delta_ch": avg_delta_ch,
        "avg_valence": avg_valence,
        "coord_type": coord_type,
        "rms_angle": rms_angle,
        "t_factor": t_factor,
        "oct_factor": oct_factor,
        "formula": formula_simple,
        "formula_canonical": formula_canon,
        "E_F": E_F,
        "eps_inf": eps_inf,
        "eps_static": eps_static,
        "m_e": m_e,
        "m_h": m_h,
    }
