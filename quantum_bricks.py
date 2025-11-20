import re
from pathlib import Path
import numpy as np
from scipy.linalg import eigh  # kept for compatibility, though we use np.linalg.eigvalsh

# ===========================================================
#   PHYSICS ENGINE v2.4-SAFE-FAST
#   (Same physics as v2.3, light precomputation for speed)
# ===========================================================

ENGINE_VERSION = "2.4-SAFE-FAST"

# -----------------------------------------------------------
#   THE BRICKS: ATOMIC PHYSICS + COULOMB REPULSION (U)
#   (Solid State Parameters - Tuned for Screening)
# -----------------------------------------------------------
ATOMIC_PARAMS = {
    "H":  {"Es": -13.6, "Ep": 0.0,   "U": 3.5}, 
    "C":  {"Es": -19.4, "Ep": -10.7, "U": 3.4}, 
    "N":  {"Es": -23.0, "Ep": -12.5, "U": 3.8}, 
    "O":  {"Es": -29.1, "Ep": -14.1, "U": 4.2}, 
    "S":  {"Es": -20.2, "Ep": -11.6, "U": 2.8}, 
    "F":  {"Es": -36.0, "Ep": -17.0, "U": 4.8},
    "Cl": {"Es": -26.0, "Ep": -14.0, "U": 3.5},
    "Br": {"Es": -23.0, "Ep": -12.0, "U": 3.0},
    "I":  {"Es": -20.0, "Ep": -10.0, "U": 2.5},
    "Au": {"Es": -9.2,  "Ep": -6.0,  "U": 0.5}, 
    "Si": {"Es": -14.7, "Ep": -7.8,  "U": 1.2},
    "Cu": {"Es": -11.4, "Ep": -6.0,  "U": 1.5},
    "Zn": {"Es": -9.39, "Ep": -5.0,  "U": 1.5},
    "Se": {"Es": -18.5, "Ep": -9.8,  "U": 2.3},
    "Te": {"Es": -17.5, "Ep": -9.0,  "U": 2.0},
    "P":  {"Es": -18.6, "Ep": -10.5, "U": 4.0}, 
    "Ge": {"Es": -15.6, "Ep": -7.5,  "U": 1.8}, 
    "As": {"Es": -17.6, "Ep": -9.1,  "U": 3.0}, 
}

# Harrison constants
ETA_SS_SIGMA = -1.40
ETA_SP_SIGMA = 1.84
ETA_PP_SIGMA = 3.24
ETA_PP_PI    = -0.81
CONST_K      = 7.62

# ===========================================================
#   CIF PARSING & SYMMETRY  (same logic as v2.3)
# ===========================================================

def apply_symmetry(atoms, ops):
    if not ops:
        return atoms

    new_atoms = []
    generated_fracs = []

    for op in ops:
        op = op.lower().replace(" ", "")
        for atom in atoms:
            x, y, z = atom['frac']
            ctx = {'x': x, 'y': y, 'z': z}

            try:
                parts = op.split(',')
                nx = eval(parts[0], {"__builtins__": None}, ctx)
                ny = eval(parts[1], {"__builtins__": None}, ctx)
                nz = eval(parts[2], {"__builtins__": None}, ctx)

                nf = np.array([nx % 1.0, ny % 1.0, nz % 1.0])

                if not any(np.linalg.norm(nf - ef) < 0.01 for ef in generated_fracs):
                    generated_fracs.append(nf)
                    new_atoms.append({"elem": atom["elem"], "frac": nf})

            except:
                pass

    return new_atoms if new_atoms else atoms


def parse_cif(path: str):
    text = Path(path).read_text()

    # --- Cell parameters ---
    cell = [10., 10., 10., 90., 90., 90.]
    keys = [
        "_cell_length_a","_cell_length_b","_cell_length_c",
        "_cell_angle_alpha","_cell_angle_beta","_cell_angle_gamma"
    ]
    for i, k in enumerate(keys):
        m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", text)
        if m:
            cell[i] = float(m.group(1))

    a, b, c, al, be, ga = cell
    al_r, be_r, ga_r = np.radians([al, be, ga])
    v = np.sqrt(
        1 - np.cos(al_r)**2 - np.cos(be_r)**2 - np.cos(ga_r)**2 +
        2*np.cos(al_r)*np.cos(be_r)*np.cos(ga_r)
    )

    lattice = np.array([
        [a, b*np.cos(ga_r),          c*np.cos(be_r)],
        [0, b*np.sin(ga_r),          c*(np.cos(al_r)-np.cos(be_r)*np.cos(ga_r))/np.sin(ga_r)],
        [0, 0,                       c*v/np.sin(ga_r)]
    ])

    atoms = []
    sym_ops = []

    lines = text.splitlines()
    in_sym = False
    in_atoms = False
    idxs = {}

    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue

        if s.startswith("loop_"):
            in_sym = False
            in_atoms = False
            idxs = {}
            continue

        if s.startswith("_symmetry_equiv_pos_as_xyz"):
            in_sym = True
            continue

        if s.startswith("_atom_site_"):
            in_atoms = True
            idxs[s.split()[0]] = len(idxs)
            continue

        if in_sym:
            op = s.replace("'", "").replace('"', "")
            if op:
                sym_ops.append(op)

        elif in_atoms:
            if s.startswith("_"):
                continue

            cols = s.split()
            try:
                l_idx = idxs.get("_atom_site_label") or idxs.get("_atom_site_type_symbol")
                x_idx = idxs.get("_atom_site_fract_x")
                y_idx = idxs.get("_atom_site_fract_y")
                z_idx = idxs.get("_atom_site_fract_z")

                if l_idx is None or x_idx is None:
                    continue

                label = cols[l_idx]
                raw = "".join([c for c in label if c.isalpha()])
                elem = raw

                if elem not in ATOMIC_PARAMS:
                    if any(c.isdigit() for c in label):
                        if elem[:1] in ATOMIC_PARAMS:
                            elem = elem[:1]
                        elif elem[:2] in ATOMIC_PARAMS:
                            elem = elem[:2]

                if elem not in ATOMIC_PARAMS:
                    continue

                fx = float(re.sub(r"\([^)]*\)", "", cols[x_idx]))
                fy = float(re.sub(r"\([^)]*\)", "", cols[y_idx]))
                fz = float(re.sub(r"\([^)]*\)", "", cols[z_idx]))

                atoms.append({"elem": elem, "frac": np.array([fx, fy, fz])})
            except:
                pass

    full_atoms = apply_symmetry(atoms, sym_ops)
    return lattice, full_atoms

# ===========================================================
#   TIGHT-BINDING SOLVER (FAST precomputation, same physics)
# ===========================================================

class TightBindingSolver:
    def __init__(self, lattice, atoms):
        self.lattice = lattice
        self.atoms = atoms
        self.n_atoms = len(atoms)

        # Precompute fractional and Cartesian coords once
        if self.n_atoms > 0:
            self.frac = np.array([a['frac'] for a in atoms])
            self.cart = self.frac @ self.lattice
        else:
            self.frac = np.zeros((0, 3))
            self.cart = np.zeros((0, 3))

        # Precompute offsets and their shifts
        self.offsets = np.array(
            [[x, y, z] for x in [-1, 0, 1]
                      for y in [-1, 0, 1]
                      for z in [-1, 0, 1]],
            dtype=float
        )  # shape (27, 3)

        self.offset_shifts = self.offsets @ self.lattice  # (27, 3)

        # Orbitals & Hubbard U
        self.orbitals = []
        self.atom_orb_map = []
        total_U = 0.0
        idx = 0

        for i, atom in enumerate(atoms):
            el = atom['elem']
            params = ATOMIC_PARAMS[el]
            total_U += params['U']

            # s orbital
            self.orbitals.append({'atom': i, 'l': 0, 'E': params['Es']})
            idx += 1

            # p orbitals, if present
            if params['Ep'] != 0.0:
                E_p = params['Ep']
                for axis in range(3):
                    self.orbitals.append({'atom': i, 'l': 1, 'axis': axis, 'E': E_p})
                idx += 3

            N = 4 if params['Ep'] != 0.0 else 1
            self.atom_orb_map.append((idx - N, idx))

        self.n_orb = len(self.orbitals)
        self.avg_U = (total_U / self.n_atoms) if self.n_atoms > 0 else 0.0

    def get_interaction(self, o1, o2, d_vec):
        d = np.linalg.norm(d_vec)
        if d < 0.1:
            return 0.0

        l, m, n = d_vec / d
        base_V = CONST_K / (d ** 2)

        if o1['l'] == 0 and o2['l'] == 0:
            return ETA_SS_SIGMA * base_V

        if o1['l'] == 0 and o2['l'] == 1:
            return [l, m, n][o2['axis']] * ETA_SP_SIGMA * base_V

        if o1['l'] == 1 and o2['l'] == 0:
            return -[l, m, n][o1['axis']] * ETA_SP_SIGMA * base_V

        if o1['l'] == 1 and o2['l'] == 1:
            a1, a2 = o1['axis'], o2['axis']
            c1, c2 = [l, m, n][a1], [l, m, n][a2]
            if a1 == a2:
                return ((c1 ** 2) * ETA_PP_SIGMA + (1 - c1 ** 2) * ETA_PP_PI) * base_V
            else:
                return (c1 * c2 * (ETA_PP_SIGMA - ETA_PP_PI)) * base_V

        return 0.0

    def build_H(self, k_vec):
        H = np.zeros((self.n_orb, self.n_orb), dtype=complex)

        # Loop structure and term ordering kept identical;
        # only shift computation is precomputed.
        for i in range(self.n_atoms):
            si, ei = self.atom_orb_map[i]

            # On-site terms
            for oi in range(si, ei):
                H[oi, oi] = self.orbitals[oi]['E']

            ri = self.cart[i]

            for j in range(self.n_atoms):
                sj, ej = self.atom_orb_map[j]
                rj = self.cart[j]

                for off_idx in range(len(self.offsets)):
                    shift = self.offset_shifts[off_idx]
                    off_vec = self.offsets[off_idx]

                    d_vec = rj + shift - ri
                    d = np.linalg.norm(d_vec)

                    if d < 0.1 or d > 4.5:
                        continue

                    # same phase definition as before
                    phase = np.exp(2j * np.pi * np.dot(k_vec, off_vec))

                    for oi in range(si, ei):
                        for oj in range(sj, ej):
                            val = self.get_interaction(self.orbitals[oi], self.orbitals[oj], d_vec)
                            H[oi, oj] += val * phase

        return H

    def solve(self):
        # Gamma point
        vals_g = np.linalg.eigvalsh(self.build_H(np.array([0.0, 0.0, 0.0])))

        # electron count (same logic)
        n_e = 0
        for a in self.atoms:
            el = a['elem']
            if el in ["H","Au","Ag","Cu","K","Na","Li"]: n_e += 1
            elif el in ["Zn","Cd","Mg","Ca"]: n_e += 2
            elif el in ["B","Al","Ga"]: n_e += 3
            elif el in ["C","Si","Ge","Sn","Pb"]: n_e += 4
            elif el in ["N","P","As","Sb"]: n_e += 5
            elif el in ["O","S","Se","Te"]: n_e += 6
            elif el in ["F","Cl","Br","I"]: n_e += 7

        if n_e % 2 != 0:
            return 0.0, 0.0, True, "Radical"

        homo = n_e // 2 - 1
        if homo >= len(vals_g) - 1:
            return 0.0, 0.0, False, "Metal"

        E_v_g = vals_g[homo]
        E_c_g = vals_g[homo + 1]

        min_gap_tb = E_c_g - E_v_g
        max_disp = 0.0
        axis = "None"

        k_scans = {"X": [0.5, 0.0, 0.0], "Y": [0.0, 0.5, 0.0], "Z": [0.0, 0.0, 0.5]}
        for label, kv in k_scans.items():
            vals_k = np.linalg.eigvalsh(self.build_H(np.array(kv)))
            gap_k = vals_k[homo + 1] - vals_k[homo]
            if gap_k < min_gap_tb:
                min_gap_tb = gap_k

            disp = abs(E_v_g - vals_k[homo])
            if disp > max_disp:
                max_disp = disp
                axis = label

        # Hubbard correction (unchanged)
        screening = min(0.8, max_disp * 1.5)
        eff_U = self.avg_U * (1.0 - screening)
        final_gap = min_gap_tb + eff_U

        return final_gap, max_disp, False, axis

# ===========================================================
#   OPTIONAL UTILITIES (SAFE)
# ===========================================================

def quick_gap(path):
    lat, atoms = parse_cif(path)
    solver = TightBindingSolver(lat, atoms)
    return solver.solve()

def export_json_results(path, outfile="tb_results.json"):
    import json
    gap, disp, _, axis = quick_gap(path)
    json.dump(
        {"gap_eV": float(gap), "max_dispersion": float(disp), "axis": axis},
        open(outfile, "w"),
        indent=2,
    )

# ===========================================================
#   CLI WRAPPER
# ===========================================================

def run(path):
    print(f"\n--- PHYSICS ENGINE v{ENGINE_VERSION} ---")
    print(f"File: {path}")
    try:
        lat, atoms = parse_cif(path)
        if not atoms:
            print("No atoms parsed.")
            return

        solver = TightBindingSolver(lat, atoms)
        gap, disp, rad, axis = solver.solve()

        if rad:
            print("! SYSTEM IS A RADICAL !   Gap = 0.000 eV")
        else:
            print(f"> Band Gap:       {gap:.3f} eV (Optical)")
            print(f"> Max Dispersion: {disp:.3f} eV")
            print(f"> Best Axis:      {axis}")
            if disp > 0.25:
                print("> Status:         HIGH MOBILITY (Speed King)")
            elif disp > 0.15:
                print("> Status:         Medium Mobility")
            else:
                print("> Status:         Low Mobility")

            # ----------------------------------------------------
            #   SAVE RESULT FILE
            # ----------------------------------------------------
            from pathlib import Path
            outpath = Path(path)
            outfile = outpath.with_suffix(".bandgap.txt")

            with open(outfile, "w") as f:
                f.write(f"File: {outpath.name}\n")
                f.write(f"Engine: v{ENGINE_VERSION}\n")
                f.write(f"Band Gap (eV): {gap:.6f}\n")
                f.write(f"Max Dispersion (eV): {disp:.6f}\n")
                f.write(f"Best Axis: {axis}\n")
                mobility = (
                    "HIGH" if disp > 0.25
                    else "MEDIUM" if disp > 0.15
                    else "LOW"
                )
                f.write(f"Mobility: {mobility}\n")
    except Exception as e:
        print(f"ERR: {e}")

if __name__ == "__main__":
    import sys, os
    paths = sys.argv[1:]

    for p in paths:
        if os.path.isdir(p):
            # Run on all .cif files in the folder
            for f in sorted(os.listdir(p)):
                if f.lower().endswith(".cif"):
                    run(os.path.join(p, f))
        else:
            run(p)
