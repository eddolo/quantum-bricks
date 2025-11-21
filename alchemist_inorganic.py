"""
AC Inorganics — Alchemist for Inorganic CIFs
v1.0 — compatible with analyze_inorganic()
"""

import re
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from qb_inorganic.analyze import analyze_inorganic

# ---------------------------------------------
# Allowed inorganic swaps (minimal safe set)
# ---------------------------------------------

ALLOWED_SWAPS_INORG = {
    "Ga": ["Al", "In"],
    "In": ["Ga", "Al"],
    "Al": ["Ga"],

    "N": ["P", "As"],
    "P": ["N", "As"],
    "As": ["N", "P"],

    "O": ["S", "Se", "Te"],
    "S": ["O", "Se", "Te"],
    "Se": ["S", "Te"],
    "Te": ["S", "Se"],

    "Mg": ["Zn", "Ca"],
    "Zn": ["Mg", "Cd"],
    "Cd": ["Zn"],

    "Sr": ["Ba", "Ca"],
    "Ba": ["Sr", "Ca"],
    "Ca": ["Sr", "Ba"],

    "Ti": ["Zr", "Hf"],
    "Zr": ["Ti", "Hf"],
    "Hf": ["Zr", "Ti"],
}

# ---------------------------------------------
# Score for inorganic semiconductor quality
# ---------------------------------------------

def score_inorganic(result):
    Eg = result.get("Eg", 0.0)
    disp = result.get("disp", 0.0)
    family = result.get("family", "")
    strain = abs(result.get("strain_pct", 0))

    # ----- gap score -----
    if Eg < 0.5:
        g = 0.0
    elif 0.5 <= Eg < 1.0:
        g = 0.3
    elif 1.0 <= Eg <= 3.0:
        g = 1.0
    elif 3.0 < Eg <= 5.0:
        g = 1.0 - (Eg - 3) / 2.0
    else:
        g = 0.0

    g = max(0.0, min(1.0, g))

    # ------ mobility ------
    m = min(1.0, disp / 1.0)  # normalize around 1.0 Å⁻¹

    # ------ family bonus ------
    if family in ["GAN", "ALN", "INP", "GAAS", "CDS", "CDSE", "MOS2"]:
        fb = 1.2
    elif "PEROVSKITE" in family:
        fb = 1.1
    else:
        fb = 1.0

    # ------ strain penalty ------
    if strain <= 3:
        sp = 1.0
    elif strain <= 10:
        sp = 1.0 - (strain - 3) * 0.05
    else:
        sp = max(0.1, 1.0 - 0.35 - (strain - 10) * 0.10)

    return max(0.0, g * m * fb * sp)

# ---------------------------------------------
# Apply mutation: text replacement
# ---------------------------------------------

def mutate_text(cif_text, e_from, e_to):
    out = re.sub(rf"\b{e_from}\b", e_to, cif_text)
    out = re.sub(rf"{e_from}(?=[0-9\(])", e_to, out)
    return out

# ---------------------------------------------
# Brutal strain (isotropic)
# ---------------------------------------------

def strain_cell(cif_text, pct):
    f = 1 + pct / 100
    def repl(m):
        tag = m.group(1)
        val = float(m.group(2))
        return f"{tag}   {val * f:.4f}"
    return re.sub(
        r"(_cell_length_[abc])\s+([0-9.]+)",
        repl,
        cif_text
    )

# ---------------------------------------------
# MAIN ENGINE
# ---------------------------------------------

def run_alchemy_inorganic(filepath, mode="both"):
    """
    Modes:
      fast  → screened only
      brute → strain only
      both  → screened + strain
    """

    p = Path(filepath)
    base = p.stem

    root = Path("ac_inorganic") / base
    root.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = root / f"{mode}_{timestamp}"
    run_dir.mkdir()

    original_text = p.read_text()

    full_log = []
    summary = {}

    def add_summary(key, res):
        summary[key] = res

    # ---- baseline ----
    res0 = analyze_inorganic(filepath)
    res0["mutation"] = "Original"
    res0["strain_pct"] = 0
    res0["score"] = score_inorganic(res0)
    full_log.append(res0)
    add_summary("Original", res0)

    # ---- detect present elements ----
    present = set()
    for el in ALLOWED_SWAPS_INORG:
        # match elements with numeric/underscore suffixes: O1, O2, Ba1, Ti2, O(1)
        pattern = rf"\b{el}(?=[^A-Za-z])"
        if re.search(pattern, original_text):
            present.add(el)

    # ---------------------------------------------
    #  LOOP: mutate each allowed element
    # ---------------------------------------------

    for e_from in present:
        for e_to in ALLOWED_SWAPS_INORG[e_from]:

            def safe_name(name):
                return re.sub(r'[^A-Za-z0-9_\-]', '_', name)

            mkey = safe_name(f"{e_from}_to_{e_to}")

            # ========== SCREENED ==========
            if mode in ("fast", "both"):
                mut_text = mutate_text(original_text, e_from, e_to)
                f = run_dir / f"{mkey}_screened.cif"
                f.write_text(mut_text)

                res = analyze_inorganic(str(f))
                res["mutation"] = mkey
                res["variant"] = "screened"
                res["strain_pct"] = 0
                res["score"] = score_inorganic(res)

                full_log.append(res)
                add_summary(mkey, res)

            # ========== BRUTE (strain sweep) ==========
            if mode in ("brute", "both"):
                best = None
                for pct in range(-8, 20, 2):  # -8% to +18%
                    mut_text = mutate_text(original_text, e_from, e_to)
                    strained = strain_cell(mut_text, pct)

                    f = run_dir / f"{mkey}_strain{pct}.cif"
                    f.write_text(strained)

                    res = analyze_inorganic(str(f))
                    res["mutation"] = mkey
                    res["variant"] = f"strain{pct}"
                    res["strain_pct"] = pct
                    res["score"] = score_inorganic(res)

                    full_log.append(res)

                    if best is None or res["score"] > best["score"]:
                        best = res

                # store best strain for summary
                if best:
                    add_summary(mkey, best)

    # ---------------------------------------------
    # EXPORT
    # ---------------------------------------------

    df_full = pd.DataFrame(full_log)
    df_full.to_csv(run_dir / f"full_log_{mode}.csv", index=False)

    df_summary = pd.DataFrame(summary.values())
    df_summary = df_summary.sort_values(by="score", ascending=False)
    df_summary.to_csv(run_dir / f"summary_{mode}.csv", index=False)

    return run_dir
