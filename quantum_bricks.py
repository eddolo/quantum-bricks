import re
import sys
import os
import copy
import numpy as np
from pathlib import Path
from scipy.linalg import eigh
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize

# ===========================================================
#   QUANTUM BRICKS v6.0 - CLEAN TB SIMULATOR
#   - Orthogonal Tight Binding (Solves Hc = Ec)
#   - Harrison-like distance scaling for hoppings
#   - Dual Analysis (Raw vs. Relaxed)
#   - Optional Auto-Relaxation (simple force field)
#   - Transparency reports
# ===========================================================

ENGINE_VERSION = "6.0-TB-ORTHO"

# --- TB PARAMETERS (still toy, but consistent) ---
ATOMIC_PARAMS = {
    "H":  {"Es": -13.6, "Ep": 0.0},
    "C":  {"Es": -21.4, "Ep": -11.4},
    "S":  {"Es": -20.2, "Ep": -11.6},
    "O":  {"Es": -32.3, "Ep": -14.8},
    "N":  {"Es": -26.0, "Ep": -13.4},
    "F":  {"Es": -40.0, "Ep": -18.1},
    "Cl": {"Es": -26.3, "Ep": -14.2},
    "Br": {"Es": -22.0, "Ep": -12.3},
    "Au": {"Es": -10.9, "Ep": -6.5},
}
DEFAULT_PARAM = {"Es": -15.0, "Ep": -8.0}

# Harrison-like prefactor (eV·Å²) for σ/π hopping
HARRISON_K = 7.62

# Cutoffs
MAX_HOP_DIST = 4.5  # Å
MIN_DIST = 0.1      # avoid self / numerical junk


class QuantumBricksEngine:
    """
    Orthogonal tight-binding engine with s + p orbitals.
    Periodic via 3×3×3 image summation, k-space scan at Γ + X/Y/Z.
    """

    def __init__(self, lattice, atoms, relax=False, dimensionality="3D"):
        self.atoms = copy.deepcopy(atoms)
        self.lattice = lattice
        self.n_atoms = len(self.atoms)
        self.dimensionality = dimensionality
        self.was_relaxed = False
        self.relax_info = "No changes made."
        self.sanity_warning = False
        self.sanity_msg = "OK"

        # 1. Geometry sanity check
        self._check_sanity()

        # 2. Optional relaxation (very crude)
        if relax:
            self._relax_geometry()

        # 3. Precompute geometry
        self.frac = np.array([a["frac"] for a in self.atoms])
        self.cart = self.frac @ self.lattice

        # 4. Neighbor offsets for periodic images
        self.offsets = np.array(
            [[x, y, z] for x in [-1, 0, 1]
                       for y in [-1, 0, 1]
                       for z in [-1, 0, 1]],
            dtype=float,
        )
        self.offset_shifts = self.offsets @ self.lattice

        # 5. Basis setup
        self._setup_orbitals()

    # -------------------------------------------------------
    #  GEOMETRY & RELAXATION
    # -------------------------------------------------------
    def _check_sanity(self, threshold=0.5):
        """Detect atoms closer than threshold (Å)."""
        temp_frac = np.array([a["frac"] for a in self.atoms])
        temp_cart = temp_frac @ self.lattice
        dists = pdist(temp_cart)
        if len(dists) > 0 and np.any(dists < threshold):
            min_d = np.min(dists)
            self.sanity_warning = True
            self.sanity_msg = f"CRITICAL: Atoms overlap at {min_d:.3f} Å."
        else:
            self.sanity_warning = False
            self.sanity_msg = "OK"

    def _relax_geometry(self):
        """
        Simple molecular-like relaxation:
        - harmonic springs for "bonds"
        - repulsion at short distances
        NOTE: This ignores periodicity in a serious way, so treat as
        a crude fixer for badly-overlapping CIFs, not a real optimizer.
        """
        start_pos = np.array([a["frac"] for a in self.atoms]) @ self.lattice
        pos = start_pos.copy()
        elements = [a["elem"] for a in self.atoms]

        def get_target(e1, e2):
            pair = sorted([e1, e2])
            if pair == ["C", "C"]:
                return 1.42
            if pair == ["C", "S"]:
                return 1.76
            if pair == ["C", "H"]:
                return 1.09
            if pair == ["C", "N"]:
                return 1.35
            if pair == ["C", "O"]:
                return 1.25
            return 1.5

        dist_mat = squareform(pdist(pos))
        bonds = []
        for i in range(len(elements)):
            for j in range(i + 1, len(elements)):
                if dist_mat[i, j] < 1.9:
                    bonds.append((i, j))

        def energy(flat_p):
            curr = flat_p.reshape((-1, 3))
            E = 0.0
            # bond springs
            for (i, j) in bonds:
                d = np.linalg.norm(curr[i] - curr[j])
                t = get_target(elements[i], elements[j])
                E += 500.0 * (d - t) ** 2
            # soft repulsion for non-bonded close pairs
            for i in range(len(curr)):
                for j in range(i + 1, len(curr)):
                    if (i, j) not in bonds:
                        d = np.linalg.norm(curr[i] - curr[j])
                        if d < 2.4:
                            E += 50.0 * (2.4 - d) ** 2
            return E

        try:
            res = minimize(
                energy,
                pos.flatten(),
                method="L-BFGS-B",
                tol=0.1,
            )
            new_pos = res.x.reshape((-1, 3))
        except Exception:
            new_pos = pos

        displacement = np.linalg.norm(new_pos - start_pos, axis=1)
        max_shift = float(np.max(displacement))
        mean_shift = float(np.mean(displacement))

        if max_shift > 0.05:
            self.was_relaxed = True
            self.relax_info = (
                "Geometry Optimized.\n"
                f"      - Max atom shift: {max_shift:.3f} Å\n"
                f"      - Mean shift:     {mean_shift:.3f} Å"
            )
            inv_lat = np.linalg.inv(self.lattice)
            new_frac = new_pos @ inv_lat
            for i, a in enumerate(self.atoms):
                a["frac"] = new_frac[i]
        else:
            self.was_relaxed = False
            self.relax_info = "Geometry was already near-optimal."

    def save_cif(self, filename):
        """Export current geometry to a simple P1 CIF."""
        with open(filename, "w") as f:
            f.write(f"data_{Path(filename).stem}\n")
            f.write("_symmetry_space_group_name_H-M 'P1'\n")
            f.write(
                "loop_\n"
                "_atom_site_label\n"
                "_atom_site_type_symbol\n"
                "_atom_site_fract_x\n"
                "_atom_site_fract_y\n"
                "_atom_site_fract_z\n"
            )
            for i, a in enumerate(self.atoms):
                el = a["elem"]
                fr = a["frac"]
                f.write(
                    f"{el}{i+1} {el} "
                    f"{fr[0]:.6f} {fr[1]:.6f} {fr[2]:.6f}\n"
                )
        return filename

    # -------------------------------------------------------
    #  ORBITALS & HARRISON HOPPINGS
    # -------------------------------------------------------
    def _setup_orbitals(self):
        """
        Minimal basis: s + (optional) px,py,pz.
        self.atom_orb_map[i] = (start_index, end_index)
        """
        self.orbitals = []
        self.atom_orb_map = []

        for i, atom in enumerate(self.atoms):
            el = atom["elem"]
            params = ATOMIC_PARAMS.get(el, DEFAULT_PARAM)
            Es = params["Es"]
            Ep = params["Ep"]

            start_idx = len(self.orbitals)

            # s
            self.orbitals.append({"atom": i, "l": 0, "E": Es})

            # p
            if Ep != 0.0:
                for axis in range(3):
                    self.orbitals.append(
                        {"atom": i, "l": 1, "axis": axis, "E": Ep}
                    )

            end_idx = len(self.orbitals)
            self.atom_orb_map.append((start_idx, end_idx))

        self.n_orb = len(self.orbitals)

    def _calc_harrison(self, o1, o2, l, m, n, d):
        """
        Harrison-like two-center hopping:
        V ~ 1/d^2 with angular factors for s/p.
        Values are toy but structurally reasonable.
        """
        if d < MIN_DIST:
            return 0.0

        base_val = HARRISON_K / (d ** 2)

        # s-s sigma
        if o1["l"] == 0 and o2["l"] == 0:
            return -1.40 * base_val

        # s-p sigma
        if o1["l"] == 0 and o2["l"] == 1:
            vec = [l, m, n][o2["axis"]]
            return 1.84 * base_val * vec

        # p-s sigma
        if o1["l"] == 1 and o2["l"] == 0:
            vec = [l, m, n][o1["axis"]]
            return -1.84 * base_val * vec

        # p-p (σ and π mix)
        if o1["l"] == 1 and o2["l"] == 1:
            a1, a2 = o1["axis"], o2["axis"]
            c1, c2 = [l, m, n][a1], [l, m, n][a2]
            V_sigma = 3.24 * base_val
            V_pi = -0.81 * base_val

            if a1 == a2:
                # same axis
                return V_sigma * (c1 ** 2) + V_pi * (1.0 - c1 ** 2)
            else:
                # cross term
                return (V_sigma - V_pi) * c1 * c2

        return 0.0

    # -------------------------------------------------------
    #  HAMILTONIAN CONSTRUCTION (ORTHOGONAL)
    # -------------------------------------------------------
    def get_hamiltonian(self, k_vec, mask_axis=None):
        """
        Build orthogonal TB Hamiltonian H(k).
        S is identity, so we solve Hc = Ec.
        """
        H = np.zeros((self.n_orb, self.n_orb), dtype=complex)

        # on-site terms
        for i in range(self.n_atoms):
            si, ei = self.atom_orb_map[i]
            for oi in range(si, ei):
                H[oi, oi] = self.orbitals[oi]["E"]

        # off-diagonal hoppings + phase factors
        for i in range(self.n_atoms):
            si, ei = self.atom_orb_map[i]
            ri = self.cart[i]

            for j in range(self.n_atoms):
                sj, ej = self.atom_orb_map[j]
                rj = self.cart[j]

                for off_idx, off_vec in enumerate(self.offsets):
                    # dimensionality mask (1D "wire" mode)
                    if mask_axis:
                        if mask_axis == "X" and (off_vec[1] != 0 or off_vec[2] != 0):
                            continue
                        if mask_axis == "Y" and (off_vec[0] != 0 or off_vec[2] != 0):
                            continue
                        if mask_axis == "Z" and (off_vec[0] != 0 or off_vec[1] != 0):
                            continue

                    shift = self.offset_shifts[off_idx]
                    d_vec = rj + shift - ri
                    d = float(np.linalg.norm(d_vec))
                    if d < MIN_DIST or d > MAX_HOP_DIST:
                        continue

                    l, m, n = d_vec / d
                    phase = np.exp(2j * np.pi * np.dot(k_vec, off_vec))

                    for oi in range(si, ei):
                        o1 = self.orbitals[oi]
                        for oj in range(sj, ej):
                            if i == j and off_idx == 13:  # central cell self
                                continue
                            o2 = self.orbitals[oj]

                            V = self._calc_harrison(o1, o2, l, m, n, d)
                            if V == 0.0:
                                continue

                            H[oi, oj] += V * phase

        # enforce Hermiticity (small numerical noise)
        H = 0.5 * (H + H.conj().T)
        return H

    # -------------------------------------------------------
    #  SOLVER
    # -------------------------------------------------------
    def _count_electrons(self):
        """Very rough valence electron count."""
        n_e = 0
        for a in self.atoms:
            el = a["elem"]
            if el == "H":
                n_e += 1
            elif el in ["S", "O"]:
                n_e += 6
            elif el in ["C", "Si"]:
                n_e += 4
            elif el in ["N", "P"]:
                n_e += 5
            elif el in ["F", "Cl", "Br"]:
                n_e += 7
            elif el == "Au":
                n_e += 1
        return n_e

    def solve(self):
        """
        Solve TB Hamiltonian at Γ + X/Y/Z and estimate:
        - bandgap (min direct gap at sampled k)
        - valence-band dispersion (as a proxy for mobility)
        """
        n_e = self._count_electrons()
        if n_e <= 0:
            return 0.0, 0.0, "No electrons", "None"
        if n_e % 2 != 0:
            return 0.0, 0.0, "Radical system (odd electrons)", "None"

        homo_idx = n_e // 2 - 1

        def solve_k(k_vec, mask=None):
            H_k = self.get_hamiltonian(k_vec, mask_axis=mask)
            vals = eigh(H_k, eigvals_only=True)
            return np.real(vals)

        # Γ reference
        try:
            vals_g = solve_k([0.0, 0.0, 0.0])
        except Exception:
            return 0.0, 0.0, "Diagonalization error at Γ", "None"

        if homo_idx + 1 >= len(vals_g):
            return 0.0, 0.0, "Metal (no empty state above HOMO)", "None"

        E_v_gamma = vals_g[homo_idx]
        gap_gamma = vals_g[homo_idx + 1] - vals_g[homo_idx]

        # k-path scan
        k_scans = {"X": [0.5, 0.0, 0.0], "Y": [0.0, 0.5, 0.0], "Z": [0.0, 0.0, 0.5]}
        min_gap = float(gap_gamma)
        max_disp = 0.0
        best_axis = "Γ"

        for label, kv in k_scans.items():
            mask = label if self.dimensionality == "1D_Isolated" else None
            try:
                vals_k = solve_k(kv, mask=mask)
            except Exception:
                continue

            if homo_idx + 1 >= len(vals_k):
                continue

            E_v_k = vals_k[homo_idx]
            gap_k = vals_k[homo_idx + 1] - vals_k[homo_idx]

            if gap_k < min_gap:
                min_gap = float(gap_k)

            disp = abs(E_v_gamma - E_v_k)
            if disp > max_disp:
                max_disp = float(disp)
                best_axis = label

        # ensure non-negative
        min_gap = max(0.0, min_gap)

        return min_gap, max_disp, "OK", best_axis


# ===========================================================
#   CIF PARSER (simple P1-ish)
# ===========================================================
def parse_cif_simple(path):
    text = Path(path).read_text()

    # cell parameters
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
        m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", text)
        if m:
            cell[i] = float(m.group(1))

    a, b, c, al, be, ga = cell
    al_r, be_r, ga_r = np.radians([al, be, ga])

    omega = np.sqrt(
        1
        - np.cos(al_r) ** 2
        - np.cos(be_r) ** 2
        - np.cos(ga_r) ** 2
        + 2 * np.cos(al_r) * np.cos(be_r) * np.cos(ga_r)
    )

    lattice = np.array(
        [
            [a, b * np.cos(ga_r), c * np.cos(be_r)],
            [
                0.0,
                b * np.sin(ga_r),
                c * (np.cos(al_r) - np.cos(be_r) * np.cos(ga_r)) / np.sin(ga_r),
            ],
            [0.0, 0.0, c * omega / np.sin(ga_r)],
        ]
    )

    atoms = []
    lines = text.splitlines()
    in_loop = False

    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith("loop_"):
            in_loop = False
            continue
        if "_atom_site_fract_z" in s:
            in_loop = True
            continue

        if in_loop and len(s.split()) >= 5:
            cols = s.split()
            try:
                el_raw = cols[1]
                el = "".join(c for c in el_raw if c.isalpha())
                if not el:
                    el = "".join(c for c in cols[0] if c.isalpha())
                x = float(cols[2])
                y = float(cols[3])
                z = float(cols[4])
                atoms.append({"elem": el, "frac": np.array([x, y, z])})
            except Exception:
                pass

    return lattice, atoms

# ===========================================================
#   COMPATIBILITY LAYER FOR STREAMLIT APP
#   (keeps old API: parse_cif + TightBindingSolver)
# ===========================================================

def parse_cif(path: str):
    """
    Backwards-compatible wrapper for the old parse_cif() name.
    Returns (lattice, atoms) from a CIF file.
    """
    return parse_cif_simple(path)


class TightBindingSolver:
    """
    Thin wrapper around QuantumBricksEngine so existing Streamlit
    code can keep calling:

        solver = TightBindingSolver(lattice, atoms)
        gap, disp, radical, axis = solver.solve()
    """

    def __init__(self, lattice, atoms, dimensionality: str = "3D", relax: bool = False):
        self.engine = QuantumBricksEngine(
            lattice=lattice,
            atoms=atoms,
            relax=relax,
            dimensionality=dimensionality,
        )

    def solve(self):
        """
        Returns:
            gap (float): band gap in eV
            disp (float): valence-band dispersion (eV)
            radical (bool): True for radical/metal/problem cases
            axis (str): axis of highest dispersion ("Γ", "X", "Y", "Z", or "None")
        """
        gap, disp, status, axis = self.engine.solve()

        # Map text status -> boolean radical flag for the app
        radical = False
        if status != "OK":
            # treat these as radical/metal-type systems
            if any(
                key in status
                for key in [
                    "Radical",
                    "odd electrons",
                    "Metal",
                    "No electrons",
                ]
            ):
                radical = True
            else:
                # diagonalization errors etc → mark as "bad" too
                radical = True

        return float(gap), float(disp), radical, axis

# ===========================================================
#   MAIN
# ===========================================================
if __name__ == "__main__":
    paths = sys.argv[1:]
    if not paths:
        print("Usage: python quantum_bricks_tb_v6.py <file.cif> [more.cif...]")
        sys.exit(1)

    print("\n╔════════════════════════════════════════════════╗")
    print(f"║  QUANTUM BRICKS v{ENGINE_VERSION:<35}║")
    print("║  Orthogonal Tight-Binding TB Simulator         ║")
    print("╚════════════════════════════════════════════════╝")

    for p in paths:
        try:
            print(f"\n📂 TARGET: {Path(p).name}")
            lat, atoms = parse_cif_simple(p)
            if not atoms:
                print("   [!] Error: No atoms found in CIF.")
                continue

            # PASS 1: raw geometry
            print("   ► [Pass 1] Raw Geometry Analysis...")
            engine_raw = QuantumBricksEngine(lat, atoms, relax=False)
            sanity_icon = "✅" if not engine_raw.sanity_warning else "⚠️"
            print(f"      Geometry Integrity: {sanity_icon} {engine_raw.sanity_msg}")

            gap1, disp1, status1, axis1 = engine_raw.solve()
            if status1 != "OK":
                print(f"      Result: ❌ FAILED ({status1})")
            else:
                print(f"      Bandgap (min sampled): {gap1:.3f} eV")
                print(f"      VB Dispersion:         {disp1:.3f} eV ({axis1})")

            # PASS 2: relaxed geometry
            print("   ► [Pass 2] Auto-Relaxation Attempt...")
            engine_fix = QuantumBricksEngine(lat, atoms, relax=True)

            if engine_fix.was_relaxed:
                new_filename = f"{Path(p).stem}_RELAXED.cif"
                engine_fix.save_cif(new_filename)
                print(f"      {engine_fix.relax_info}")
                print(f"      💾 Exported fix to: {new_filename}")

                gap2, disp2, status2, axis2 = engine_fix.solve()
                if status2 != "OK":
                    print(f"      New Result: ❌ FAILED ({status2})")
                else:
                    print(f"      New Bandgap:          {gap2:.3f} eV")
                    print(f"      New VB Dispersion:    {disp2:.3f} eV ({axis2})")
                    if disp2 > 0.25:
                        print("      🏃 High-mobility-like (broad valence band)")
                    else:
                        print("      Standard / localized-like")
            else:
                print("      ✅ Geometry already near-optimal. Pass 1 result is fine.")

        except Exception as e:
            print(f"   [!] System Error: {e}")
