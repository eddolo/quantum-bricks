import sys
import os
import re
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.optimize import minimize

from quantum_bricks import parse_cif, TightBindingSolver, ATOMIC_PARAMS

# ==========================================
#   CONFIGURATION
# ==========================================

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "S": 1.80, "Cl": 1.75, "Si": 2.10, "P": 1.80,
    "Se": 1.90, "Br": 1.85, "Ge": 2.11, "Te": 2.06, "As": 1.85
}

ALLOWED_SWAPS = {
    "H": ["F", "Cl"], "F": ["H", "Cl"], "Cl": ["F", "H"],
    "O": ["S", "Se"], "S": ["O", "Se", "Te"], "Se": ["S", "Te"], "Te": ["Se", "S"],
    "N": ["P", "As"], "P": ["N", "As"], "As": ["P", "N"],
    "C": ["Si", "Ge"], "Si": ["C", "Ge"], "Ge": ["Si", "C"],
}

# ==========================================
#   SEMICONDUCTOR SCORING
# ==========================================

def calculate_viability_score(disp, gap, strain_pct):
    """
    Combined semiconductor score:
    - GapScore: best for Eg in [1, 3] eV, smoothly decays outside
    - MobilityScore: based on dispersion, saturates at ~0.5 eV
    - Strain penalty: same idea as v40, but applied after physics
    """

    # ------------------
    # 1) Gap quality
    # ------------------
    try:
        Eg = float(gap)
    except (TypeError, ValueError):
        return 0.0

    if Eg < 0.5:
        gap_score = 0.0
    elif Eg < 1.0:
        # ramp from 0 to 0.5 between 0.5 and 1.0 eV
        gap_score = 0.5 * (Eg - 0.5) / 0.5
    elif Eg <= 3.0:
        # ideal semiconductor window
        gap_score = 1.0
    elif Eg <= 5.0:
        # smooth decay from 1 to 0 between 3 and 5 eV
        gap_score = 1.0 - (Eg - 3.0) / 2.0
    else:
        gap_score = 0.0

    gap_score = max(0.0, min(1.0, gap_score))

    # ------------------
    # 2) Mobility score from dispersion
    # ------------------
    try:
        d = float(disp)
    except (TypeError, ValueError):
        d = 0.0

    if d <= 0:
        mobility_score = 0.0
    else:
        mobility_score = min(1.0, d / 0.5)  # saturate at ~0.5 eV

    # ------------------
    # 3) Strain penalty
    # ------------------
    try:
        s = abs(float(strain_pct))
    except (TypeError, ValueError):
        s = 0.0

    if s <= 3:
        penalty = 0.0
    elif s <= 10:
        penalty = (s - 3) * 0.05
    else:
        penalty = 0.35 + (s - 10) * 0.10

    penalty = min(0.90, penalty)  # cap

    score = gap_score * mobility_score * (1.0 - penalty)
    return max(score, 0.0)

# ==========================================
#   HELPER FUNCTIONS
# ==========================================

def get_chemical_formula(atoms):
    """Generates the Hill system chemical formula (C first, H second, then alphabetical)."""
    counts = {}
    for a in atoms:
        el = a['elem']
        counts[el] = counts.get(el, 0) + 1

    def hill_sort(e):
        if e == 'C':
            return (0, e)
        if e == 'H':
            return (1, e)
        return (2, e)

    elements = sorted(counts.keys(), key=hill_sort)
    formula = ""
    for el in elements:
        c = counts[el]
        formula += f"{el}" + (str(c) if c > 1 else "")
    return formula


def generate_mutant_text(cif_text, e_from, e_to):
    """
    Replace element labels in CIF text for mutation, keeping numbers/indices.
    Two patterns:
    - e_from followed by digits/parentheses
    - standalone element with spaces around it
    """
    mutant = re.sub(rf"\b{e_from}(?=[0-9\(])", e_to, cif_text)
    mutant = re.sub(rf"(\s){e_from}(\s)", rf"\1{e_to}\2", mutant)
    return mutant


def apply_new_cell_params(cif_text, params):
    """
    Replace _cell_length_[abc] and _cell_angle_[alpha|beta|gamma] with new values.
    """
    a, b, c, al, be, ga = params

    def replace(match):
        tag = match.group(1)
        if "_length_a" in tag:
            val = a
        elif "_length_b" in tag:
            val = b
        elif "_length_c" in tag:
            val = c
        elif "_angle_alpha" in tag:
            val = al
        elif "_angle_beta" in tag:
            val = be
        elif "_angle_gamma" in tag:
            val = ga
        else:
            return match.group(0)
        return f"{tag}   {val:.4f}"

    return re.sub(
        r"(_cell_length_[abc]|_cell_angle_[a-z]+)\s+([0-9.]+(?:\([0-9]+\))?)",
        replace,
        cif_text
    )


def _parse_atoms_from_text(text):
    atoms = []
    lines = text.splitlines()
    in_loop = False
    idxs = {}
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("loop_"):
            in_loop = False
            idxs = {}
            continue
        if s.startswith("_atom_site_"):
            in_loop = True
            idxs[s.split()[0]] = len(idxs)
            continue
        if in_loop and not s.startswith("_"):
            cols = s.split()
            try:
                l_idx = idxs.get("_atom_site_label") or idxs.get("_atom_site_type_symbol")
                if l_idx is not None:
                    el = "".join([c for c in cols[l_idx] if c.isalpha()])
                    if el not in ATOMIC_PARAMS:
                        el = el[:1]
                    if el in ATOMIC_PARAMS:
                        fx, fy, fz = [
                            float(re.sub(r"\([^)]*\)", "", cols[idxs[k]]))
                            for k in [
                                "_atom_site_fract_x",
                                "_atom_site_fract_y",
                                "_atom_site_fract_z",
                            ]
                        ]
                        atoms.append({"elem": el, "frac": np.array([fx, fy, fz])})
            except Exception:
                pass
    return atoms


def build_lattice_matrix(params):
    a, b, c, al, be, ga = params
    al_r, be_r, ga_r = np.radians([al, be, ga])
    try:
        v = np.sqrt(
            1
            - np.cos(al_r) ** 2
            - np.cos(be_r) ** 2
            - np.cos(ga_r) ** 2
            + 2 * np.cos(al_r) * np.cos(be_r) * np.cos(ga_r)
        )
        return np.array(
            [
                [a, b * np.cos(ga_r), c * np.cos(be_r)],
                [
                    0,
                    b * np.sin(ga_r),
                    c * (np.cos(al_r) - np.cos(be_r) * np.cos(ga_r))
                    / np.sin(ga_r),
                ],
                [0, 0, c * v / np.sin(ga_r)],
            ]
        )
    except Exception:
        return np.eye(3) * 10.0


def collision_cost(lengths, angles, atoms):
    a, b, c = lengths
    al, be, ga = angles
    lat = build_lattice_matrix([a, b, c, al, be, ga])
    fracs = np.array([a_['frac'] for a_ in atoms])
    cart = fracs @ lat
    overlap = 0.0
    offsets = [
        np.array([0, 0, 0]),
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1]),
    ]
    for i in range(len(atoms)):
        for j in range(i, len(atoms)):
            limit = (
                VDW_RADII.get(atoms[i]['elem'], 1.5)
                + VDW_RADII.get(atoms[j]['elem'], 1.5)
            ) * 0.70
            for off in offsets:
                if i == j and np.linalg.norm(off) < 0.1:
                    continue
                d = np.linalg.norm(cart[j] + (off @ lat) - cart[i])
                if 0.1 < d < limit:
                    overlap += (limit - d) ** 2
    return overlap


def relax_geometry(cif_text):
    """
    Very light L-BFGS relaxation of cell lengths only.
    """
    cell = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    keys = [
        "_cell_length_a",
        "_cell_length_b",
        "_cell_length_c",
        "_cell_angle_alpha",
        "_cell_angle_beta",
        "_cell_angle_gamma",
    ]
    for i, k in enumerate(keys):
        m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", cif_text)
        if m:
            cell[i] = float(re.sub(r"\(.*\)", "", m.group(1)))
    atoms = _parse_atoms_from_text(cif_text)
    if not atoms:
        return cell
    res = minimize(
        collision_cost,
        cell[:3],
        args=(cell[3:], atoms),
        method="L-BFGS-B",
        bounds=[(l * 0.8, l * 1.4) for l in cell[:3]],
        options={"maxiter": 30},
    )
    return list(res.x) + cell[3:]


def print_live_leaderboard(log_data, mode):
    """
    Refreshes terminal leaderboard for current run.
    Shows top candidates by Score for the selected mode.
    """
    if not log_data:
        return
    df = pd.DataFrame(log_data)
    if "Score" not in df.columns:
        return

    df_sorted = df.sort_values(by="Score", ascending=False)

    os.system('cls' if os.name == 'nt' else 'clear')
    print(f"--- ⚗️  ALCHEMIST v41 ({mode.upper()} MODE) ---")
    print(f"Candidates: {len(log_data)}")
    print("=" * 110)
    print(
        f"{'TYPE':<10} | {'FORMULA':<16} | {'VARIANT':<12} | "
        f"{'DISP(eV)':<9} | {'GAP(eV)':<8} | {'SCORE':<7} | {'STRAIN%':<7}"
    )
    print("=" * 110)

    for _, row in df_sorted.head(8).iterrows():
        strain = row.get("Strain_Pct", "")
        if isinstance(strain, (int, float)):
            strain_str = f"{strain:+d}"
        elif strain is None:
            strain_str = ""
        else:
            strain_str = str(strain)

        print(
            f"{row['Type']:<10} | "
            f"{str(row['Formula'])[:16]:<16} | "
            f"{str(row.get('Variant',''))[:12]:<12} | "
            f"{row['Disp']:<9.3f} | "
            f"{row['Gap']:<8.3f} | "
            f"{row['Score']:<7.3f} | "
            f"{strain_str:<7}"
        )
    print("=" * 110)


# ==========================================
#   MAIN ALCHEMY ENGINE (v41.0)
# ==========================================

def run_alchemy(filepath, mode):
    """
    Main entry for a single CIF file and a single mode.
    Modes:
      - fast  : SCREENED only
      - relax : RELAXED only
      - brute : BRUTE only
      - both  : RELAXED + BRUTE
    """

    path_obj = Path(filepath)
    base_name = path_obj.stem

    # Folder: alchemist/<cif_basename>/
    root_dir = Path("alchemist") / base_name
    root_dir.mkdir(parents=True, exist_ok=True)

    # Mode + timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = root_dir / f"{mode}_{timestamp}"
    run_dir.mkdir(exist_ok=True)

    print(f"\n=== ALCHEMIST v41 | File: {filepath} | Mode: {mode} ===")
    print(f"Output dir: {run_dir}")

    original_text = path_obj.read_text()
    full_log = []     # all trials
    summary_log = {}  # one row per mutation

    # ---------------------------------------------------
    #   BASELINE: ORIGINAL (depends on mode)
    # ---------------------------------------------------
    # Parse original atoms for formula
    try:
        lat_orig, atoms_orig = parse_cif(filepath)
        orig_formula = get_chemical_formula(atoms_orig) if atoms_orig else "?"
    except Exception:
        lat_orig, atoms_orig, orig_formula = None, [], "?"

    # Helper to store summary row if not exists
    def ensure_summary(mtype, base_type):
        if mtype not in summary_log:
            summary_log[mtype] = {
                "Type": base_type,
                "Formula": orig_formula if mtype == "Original" else "?",
            }

    # FAST: only SCREENED original
    if mode == "fast":
        try:
            fname = run_dir / "Original_fast.cif"
            fname.write_text(original_text)
            lat, atoms = parse_cif(str(fname))
            formula = get_chemical_formula(atoms) if atoms else orig_formula
            solver = TightBindingSolver(lat, atoms)
            gap, disp, _, _ = solver.solve()
            score = calculate_viability_score(disp, gap, 0)

            full_log.append({
                "Filename": fname.name,
                "Type": "Original",
                "Formula": formula,
                "Variant": "screened",
                "Disp": float(disp),
                "Gap": float(gap),
                "Score": float(score),
                "Strain_Pct": 0,
            })

            ensure_summary("Original", "Original")
            summary_log["Original"].update({
                "Formula": formula,
                "Disp_Screened": disp,
                "Gap_Screened": gap,
                "Score_Screened": score,
            })
            print_live_leaderboard(full_log, mode)
        except Exception as e:
            print(f"[Original FAST] Error: {e}")

    # RELAX: only RELAXED original
    if mode == "relax" or mode == "both":
        try:
            relaxed_params = relax_geometry(original_text)
            relaxed_text = apply_new_cell_params(original_text, relaxed_params)
            fname = run_dir / f"Original_relaxed_{mode}.cif"
            fname.write_text(relaxed_text)

            lat, atoms = parse_cif(str(fname))
            formula = get_chemical_formula(atoms) if atoms else orig_formula
            solver = TightBindingSolver(lat, atoms)
            gap_r, disp_r, _, _ = solver.solve()
            score_r = calculate_viability_score(disp_r, gap_r, 0)

            full_log.append({
                "Filename": fname.name,
                "Type": "Original",
                "Formula": formula,
                "Variant": "relaxed",
                "Disp": float(disp_r),
                "Gap": float(gap_r),
                "Score": float(score_r),
                "Strain_Pct": 0,
            })

            ensure_summary("Original", "Original")
            summary_log["Original"].update({
                "Formula": formula,
                "Disp_Relaxed": disp_r,
                "Gap_Relaxed": gap_r,
                "Score_Relaxed": score_r,
            })
            print_live_leaderboard(full_log, mode)
        except Exception as e:
            print(f"[Original RELAX] Error: {e}")

    # BRUTE: only BRUTE original (or with BOTH)
    if mode == "brute" or mode == "both":
        try:
            # Extract original cell
            cell = [10., 10., 10., 90., 90., 90.]
            keys = [
                "_cell_length_a", "_cell_length_b", "_cell_length_c",
                "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"
            ]
            for i, k in enumerate(keys):
                m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", original_text)
                if m:
                    cell[i] = float(re.sub(r"\(.*\)", "", m.group(1)))

            best_score = 0.0
            best_disp = 0.0
            best_gap = 0.0
            best_strain = 0

            for i in range(-5, 16):  # -10% to +30% in steps of 2%
                exp = 1.0 + (i * 0.02)
                exp_pct = int(round((exp - 1.0) * 100))
                brute_params = [p * exp for p in cell[:3]] + cell[3:]
                brute_text = apply_new_cell_params(original_text, brute_params)
                fname = run_dir / f"Original_Exp{exp_pct}_{mode}.cif"
                fname.write_text(brute_text)

                try:
                    lat, atoms = parse_cif(str(fname))
                    formula = get_chemical_formula(atoms) if atoms else orig_formula
                    solver = TightBindingSolver(lat, atoms)
                    gap_b, disp_b, _, _ = solver.solve()
                    score_b = calculate_viability_score(disp_b, gap_b, exp_pct)

                    full_log.append({
                        "Filename": fname.name,
                        "Type": "Original",
                        "Formula": formula,
                        "Variant": f"Exp{exp_pct}",
                        "Disp": float(disp_b),
                        "Gap": float(gap_b),
                        "Score": float(score_b),
                        "Strain_Pct": int(exp_pct),
                    })
                    print_live_leaderboard(full_log, mode)

                    if score_b > best_score:
                        best_score = score_b
                        best_disp = disp_b
                        best_gap = gap_b
                        best_strain = exp_pct
                except Exception:
                    continue

            ensure_summary("Original", "Original")
            summary_log["Original"].update({
                "Formula": orig_formula,
                "Disp_Brute": best_disp,
                "Gap_Brute": best_gap,
                "Score_Brute": best_score,
                "Strain_Brute": best_strain,
            })
        except Exception as e:
            print(f"[Original BRUTE] Error: {e}")

    # ---------------------------------------------------
    #   MUTATION LOOP
    # ---------------------------------------------------
    present = set()
    for el in ALLOWED_SWAPS:
        # detect if element appears in labels or as standalone
        if re.search(rf"\b{el}(?=[0-9\(])", original_text) or re.search(rf"(\s){el}(\s)", original_text):
            present.add(el)

    # Pre-calc original cell for brute
    cell_orig = [10., 10., 10., 90., 90., 90.]
    keys = [
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"
    ]
    for i, k in enumerate(keys):
        m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", original_text)
        if m:
            cell_orig[i] = float(re.sub(r"\(.*\)", "", m.group(1)))

    for e_from in present:
        for e_to in ALLOWED_SWAPS.get(e_from, []):
            mtype = f"{e_from}->{e_to}"
            ensure_summary(mtype, mtype)

            # FAST: SCREENED only
            if mode == "fast":
                raw_mutant = generate_mutant_text(original_text, e_from, e_to)
                fname = run_dir / f"{e_from}-{e_to}_Mutant_fast.cif"
                fname.write_text(raw_mutant)
                try:
                    lat, atoms = parse_cif(str(fname))
                    formula = get_chemical_formula(atoms) if atoms else "?"
                    solver = TightBindingSolver(lat, atoms)
                    gap, disp, _, _ = solver.solve()
                    score = calculate_viability_score(disp, gap, 0)

                    full_log.append({
                        "Filename": fname.name,
                        "Type": mtype,
                        "Formula": formula,
                        "Variant": "screened",
                        "Disp": float(disp),
                        "Gap": float(gap),
                        "Score": float(score),
                        "Strain_Pct": 0,
                    })
                    summary_log[mtype].update({
                        "Formula": formula,
                        "Disp_Screened": disp,
                        "Gap_Screened": gap,
                        "Score_Screened": score,
                    })
                    print_live_leaderboard(full_log, mode)
                except Exception:
                    continue

            # RELAX: RELAXED only (or BOTH)
            if mode == "relax" or mode == "both":
                raw_mutant = generate_mutant_text(original_text, e_from, e_to)
                relaxed_params = relax_geometry(raw_mutant)
                relaxed_text = apply_new_cell_params(raw_mutant, relaxed_params)
                fname = run_dir / f"{e_from}-{e_to}_Mutant_relaxed_{mode}.cif"
                fname.write_text(relaxed_text)
                try:
                    lat, atoms = parse_cif(str(fname))
                    formula = get_chemical_formula(atoms) if atoms else "?"
                    solver = TightBindingSolver(lat, atoms)
                    gap_r, disp_r, _, _ = solver.solve()
                    score_r = calculate_viability_score(disp_r, gap_r, 0)

                    full_log.append({
                        "Filename": fname.name,
                        "Type": mtype,
                        "Formula": formula,
                        "Variant": "relaxed",
                        "Disp": float(disp_r),
                        "Gap": float(gap_r),
                        "Score": float(score_r),
                        "Strain_Pct": 0,
                    })
                    summary_log[mtype].update({
                        "Formula": formula,
                        "Disp_Relaxed": disp_r,
                        "Gap_Relaxed": gap_r,
                        "Score_Relaxed": score_r,
                    })
                    print_live_leaderboard(full_log, mode)
                except Exception:
                    continue

            # BRUTE: BRUTE only (or BOTH)
            if mode == "brute" or mode == "both":
                raw_mutant = generate_mutant_text(original_text, e_from, e_to)
                best_score = 0.0
                best_disp = 0.0
                best_gap = 0.0
                best_strain = 0
                best_formula = "?"

                for i in range(-5, 16):  # -10% to +30% in steps of 2%
                    exp = 1.0 + (i * 0.02)
                    exp_pct = int(round((exp - 1.0) * 100))
                    brute_params = [p * exp for p in cell_orig[:3]] + cell_orig[3:]
                    brute_text = apply_new_cell_params(raw_mutant, brute_params)
                    fname = run_dir / f"{e_from}-{e_to}_Mutant_Exp{exp_pct}_{mode}.cif"
                    fname.write_text(brute_text)

                    try:
                        lat, atoms = parse_cif(str(fname))
                        formula = get_chemical_formula(atoms) if atoms else "?"
                        solver = TightBindingSolver(lat, atoms)
                        gap_b, disp_b, _, _ = solver.solve()
                        score_b = calculate_viability_score(disp_b, gap_b, exp_pct)

                        full_log.append({
                            "Filename": fname.name,
                            "Type": mtype,
                            "Formula": formula,
                            "Variant": f"Exp{exp_pct}",
                            "Disp": float(disp_b),
                            "Gap": float(gap_b),
                            "Score": float(score_b),
                            "Strain_Pct": int(exp_pct),
                        })
                        print_live_leaderboard(full_log, mode)

                        if score_b > best_score:
                            best_score = score_b
                            best_disp = disp_b
                            best_gap = gap_b
                            best_strain = exp_pct
                            best_formula = formula
                    except Exception:
                        continue

                summary_log[mtype].update({
                    "Formula": best_formula,
                    "Disp_Brute": best_disp,
                    "Gap_Brute": best_gap,
                    "Score_Brute": best_score,
                    "Strain_Brute": best_strain,
                })

    # ---------------------------------------------------
    #   EXPORT LOGS
    # ---------------------------------------------------
    if full_log:
        df_log = pd.DataFrame(full_log)
        df_log.to_csv(run_dir / f"leaderboard_full_log_{mode}.csv", index=False)

    if summary_log:
        summary_df = pd.DataFrame.from_dict(summary_log, orient='index')
        summary_df.index.name = "Mutation"

        # Choose columns and sort key depending on mode
        if mode == "fast":
            cols = ["Type", "Formula", "Disp_Screened", "Gap_Screened", "Score_Screened"]
            sort_col = "Score_Screened"
        elif mode == "relax":
            cols = ["Type", "Formula", "Disp_Relaxed", "Gap_Relaxed", "Score_Relaxed"]
            sort_col = "Score_Relaxed"
        elif mode == "brute":
            cols = ["Type", "Formula", "Disp_Brute", "Gap_Brute", "Score_Brute", "Strain_Brute"]
            sort_col = "Score_Brute"
        else:  # both
            cols = [
                "Type", "Formula",
                "Disp_Relaxed", "Gap_Relaxed", "Score_Relaxed",
                "Disp_Brute", "Gap_Brute", "Score_Brute", "Strain_Brute"
            ]
            sort_col = "Score_Brute"

        # Keep only columns that exist
        cols = [c for c in cols if c in summary_df.columns]
        summary_df = summary_df[cols]

        if sort_col in summary_df.columns:
            summary_df = summary_df.sort_values(by=sort_col, ascending=False)

        csv_path = run_dir / f"summary_{mode}.csv"
        summary_df.to_csv(csv_path)

        # Human-readable TXT summary
        txt_path = run_dir / f"summary_{mode}.txt"
        with txt_path.open("w", encoding="utf-8") as f:
            f.write(f"ALCHEMIST v41 SUMMARY | Mode: {mode} | File: {filepath}\n")
            f.write("=" * 80 + "\n\n")
            for idx, row in summary_df.iterrows():
                f.write(f"Mutation : {idx}\n")
                f.write(f"Type     : {row.get('Type', '')}\n")
                f.write(f"Formula  : {row.get('Formula', '')}\n\n")
                if mode == "fast":
                    f.write("  SCREENED\n")
                    f.write(f"    Disp  : {row.get('Disp_Screened', 0):.4f} eV\n")
                    f.write(f"    Gap   : {row.get('Gap_Screened', 0):.4f} eV\n")
                    f.write(f"    Score : {row.get('Score_Screened', 0):.4f}\n")
                elif mode == "relax":
                    f.write("  RELAXED\n")
                    f.write(f"    Disp  : {row.get('Disp_Relaxed', 0):.4f} eV\n")
                    f.write(f"    Gap   : {row.get('Gap_Relaxed', 0):.4f} eV\n")
                    f.write(f"    Score : {row.get('Score_Relaxed', 0):.4f}\n")
                elif mode == "brute":
                    f.write("  BRUTE (best)\n")
                    f.write(f"    Disp   : {row.get('Disp_Brute', 0):.4f} eV\n")
                    f.write(f"    Gap    : {row.get('Gap_Brute', 0):.4f} eV\n")
                    f.write(f"    Score  : {row.get('Score_Brute', 0):.4f}\n")
                    f.write(f"    Strain : {row.get('Strain_Brute', 0)} %\n")
                else:  # both
                    f.write("  RELAXED\n")
                    f.write(f"    Disp  : {row.get('Disp_Relaxed', 0):.4f} eV\n")
                    f.write(f"    Gap   : {row.get('Gap_Relaxed', 0):.4f} eV\n")
                    f.write(f"    Score : {row.get('Score_Relaxed', 0):.4f}\n\n")
                    f.write("  BRUTE (best)\n")
                    f.write(f"    Disp   : {row.get('Disp_Brute', 0):.4f} eV\n")
                    f.write(f"    Gap    : {row.get('Gap_Brute', 0):.4f} eV\n")
                    f.write(f"    Score  : {row.get('Score_Brute', 0):.4f}\n")
                    f.write(f"    Strain : {row.get('Strain_Brute', 0)} %\n")

                best_score = 0.0
                if mode == "fast":
                    best_score = row.get("Score_Screened", 0.0)
                elif mode == "relax":
                    best_score = row.get("Score_Relaxed", 0.0)
                elif mode == "brute":
                    best_score = row.get("Score_Brute", 0.0)
                else:
                    best_score = max(
                        row.get("Score_Relaxed", 0.0),
                        row.get("Score_Brute", 0.0),
                    )
                f.write(f"\n  BEST SCORE: {best_score:.4f}\n")
                f.write("-" * 80 + "\n\n")

        print("\n" + "=" * 80)
        print(f"FINAL SUMMARY ({mode}) written to: {csv_path}")
        print(f"TXT summary: {txt_path}")
        print("=" * 80)


# ==========================================
#   CLI ENTRY POINT
# ==========================================

if __name__ == "__main__":
    args = sys.argv[1:]
    mode = "brute"  # default

    if "--both" in args:
        mode = "both"
        args.remove("--both")
    elif "--relax" in args:
        mode = "relax"
        args.remove("--relax")
    elif "--fast" in args:
        mode = "fast"
        args.remove("--fast")
    elif "--brute" in args:
        # explicit, but same as default
        args.remove("--brute")

    files = args
    if not files:
        print("Usage: python alchemist.py <cif_file_or_folder> [--fast|--brute|--relax|--both]")
        sys.exit(0)

    for f in files:
        p = Path(f)
        if p.is_dir():
            for name in sorted(os.listdir(p)):
                if name.lower().endswith(".cif"):
                    run_alchemy(str(p / name), mode)
        else:
            run_alchemy(str(p), mode)
