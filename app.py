import streamlit as st
import re
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from io import StringIO

# ==========================================
#  PHYSICS ENGINE CONFIGURATION
# ==========================================

st.set_page_config(
    page_title="Quantum Bricks",
    page_icon="⚛️",
    layout="wide"
)

# Solid State Atomic Parameters (Screened U)
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
}

# Harrison's Constants
ETA_SS_SIGMA = -1.40
ETA_SP_SIGMA = 1.84
ETA_PP_SIGMA = 3.24
ETA_PP_PI    = -0.81
CONST_K      = 7.62 

# ==========================================
#  CORE PHYSICS LOGIC (v2.2)
# ==========================================

def apply_symmetry(atoms, ops):
    if not ops: return atoms
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
                
                is_duplicate = False
                for ef in generated_fracs:
                    if np.linalg.norm(nf - ef) < 0.01: 
                        is_duplicate = True
                        break
                if not is_duplicate:
                    generated_fracs.append(nf)
                    new_atoms.append({"elem": atom["elem"], "frac": nf})
            except: pass
    if not new_atoms: return atoms
    return new_atoms

def parse_cif_content(cif_text):
    """Parses CIF content from a string."""
    cell_params = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    keys = ["_cell_length_a", "_cell_length_b", "_cell_length_c",
            "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]
    for i, k in enumerate(keys):
        m = re.search(rf"{re.escape(k)}\s+([0-9.]+)", cif_text)
        if m: cell_params[i] = float(m.group(1))
    
    a, b, c, al, be, ga = cell_params
    al, be, ga = np.radians([al, be, ga])
    v = np.sqrt(1 - np.cos(al)**2 - np.cos(be)**2 - np.cos(ga)**2 + 2*np.cos(al)*np.cos(be)*np.cos(ga))
    
    lattice = np.array([
        [a, b * np.cos(ga), c * np.cos(be)],
        [0.0, b * np.sin(ga), c * (np.cos(al) - np.cos(be)*np.cos(ga))/np.sin(ga)],
        [0.0, 0.0, c * v / np.sin(ga)]
    ])
    
    sym_ops = []
    atoms = []
    lines = cif_text.splitlines()
    in_sym_loop = False
    in_atom_loop = False
    loop_indices = {}
    
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"): continue
        
        if s.startswith("loop_"):
            in_sym_loop = False
            in_atom_loop = False
            loop_indices = {}
            continue
        if s.startswith("_symmetry_equiv_pos_as_xyz"):
            in_sym_loop = True
            continue
        if s.startswith("_atom_site_"):
            in_atom_loop = True
            parts = s.split()
            loop_indices[parts[0]] = len(loop_indices)
            continue
            
        if in_sym_loop:
            op = s.replace("'", "").replace('"', "")
            if op: sym_ops.append(op)
        elif in_atom_loop:
            if s.startswith("_"): continue
            cols = s.split()
            try:
                lab_idx = loop_indices.get("_atom_site_label")
                if lab_idx is None: lab_idx = loop_indices.get("_atom_site_type_symbol")
                x_idx = loop_indices.get("_atom_site_fract_x")
                y_idx = loop_indices.get("_atom_site_fract_y")
                z_idx = loop_indices.get("_atom_site_fract_z")
                
                if lab_idx is not None:
                    label = cols[lab_idx]
                    elem = "".join([c for c in label if c.isalpha()])
                    if elem not in ATOMIC_PARAMS: elem = elem[:1]
                    if elem not in ATOMIC_PARAMS: elem = elem[:2]
                    
                    if elem in ATOMIC_PARAMS:
                        fx = float(re.sub(r"\([^)]*\)", "", cols[x_idx]))
                        fy = float(re.sub(r"\([^)]*\)", "", cols[y_idx]))
                        fz = float(re.sub(r"\([^)]*\)", "", cols[z_idx]))
                        atoms.append({"elem": elem, "frac": np.array([fx, fy, fz])})
            except: pass

    full_atoms = apply_symmetry(atoms, sym_ops)
    return lattice, full_atoms

class TightBindingSolver:
    def __init__(self, lattice, atoms):
        self.lattice = lattice
        self.atoms = atoms
        self.n_atoms = len(atoms)
        self.orbitals = [] 
        self.atom_orb_map = []
        self.avg_U = 0.0
        idx = 0
        
        total_U = 0
        for atom in atoms: total_U += ATOMIC_PARAMS[atom['elem']]['U']
        if self.n_atoms > 0: self.avg_U = total_U / self.n_atoms

        for i, atom in enumerate(atoms):
            el = atom['elem']
            params = ATOMIC_PARAMS[el]
            self.orbitals.append({'atom': i, 'l': 0, 'E': params['Es']})
            idx += 1
            if params['Ep'] != 0.0:
                E_p = params['Ep']
                self.orbitals += [{'atom': i, 'l': 1, 'axis': x, 'E': E_p} for x in range(3)]
                idx += 3
            self.atom_orb_map.append((idx - (4 if params['Ep']!=0 else 1), idx))
            
        self.n_orb = len(self.orbitals)

    def get_interaction(self, o1, o2, d_vec):
        d = np.linalg.norm(d_vec)
        if d < 0.1: return 0.0
        l, m, n = d_vec / d
        base_V = CONST_K / (d**2)
        if o1['l']==0 and o2['l']==0: return ETA_SS_SIGMA * base_V
        if o1['l']==0 and o2['l']==1: return [l,m,n][o2['axis']] * ETA_SP_SIGMA * base_V
        if o1['l']==1 and o2['l']==0: return -[l,m,n][o1['axis']] * ETA_SP_SIGMA * base_V
        if o1['l']==1 and o2['l']==1:
            a1, a2 = o1['axis'], o2['axis']
            c1, c2 = [l,m,n][a1], [l,m,n][a2]
            if a1 == a2: return ((c1**2)*ETA_PP_SIGMA + (1-c1**2)*ETA_PP_PI) * base_V
            else: return (c1*c2*(ETA_PP_SIGMA - ETA_PP_PI)) * base_V
        return 0.0

    def build_H(self, k_vec):
        H = np.zeros((self.n_orb, self.n_orb), dtype=complex)
        offsets = [np.array([x,y,z]) for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1]]
        frac = np.array([a['frac'] for a in self.atoms])
        cart = frac @ self.lattice
        
        for i in range(self.n_atoms):
            start_i, end_i = self.atom_orb_map[i]
            for k in range(start_i, end_i): H[k,k] = self.orbitals[k]['E']
            
            for j in range(self.n_atoms):
                start_j, end_j = self.atom_orb_map[j]
                r_j_base = cart[j]
                for off in offsets:
                    shift = off @ self.lattice
                    d_vec = r_j_base + shift - cart[i]
                    d = np.linalg.norm(d_vec)
                    if d > 4.5 or d < 0.1: continue
                    phase = np.exp(2j * np.pi * np.dot(k_vec, off))
                    for oi in range(start_i, end_i):
                        for oj in range(start_j, end_j):
                            val = self.get_interaction(self.orbitals[oi], self.orbitals[oj], d_vec)
                            H[oi, oj] += val * phase
        return H

    def solve(self):
        # Gamma Point
        vals_g = np.linalg.eigvalsh(self.build_H(np.array([0.,0.,0.])))
        
        n_elec = 0
        for a in self.atoms:
            el = a['elem']
            if el in ["H","Au","Ag","Cu","K","Na","Li"]: n_elec+=1
            elif el in ["Zn","Cd","Mg","Ca"]: n_elec+=2
            elif el in ["B","Al","Ga"]: n_elec+=3
            elif el in ["C","Si","Ge","Sn","Pb"]: n_elec+=4
            elif el in ["N","P","As","Sb"]: n_elec+=5
            elif el in ["O","S","Se","Te"]: n_elec+=6
            elif el in ["F","Cl","Br","I"]: n_elec+=7
            
        if n_elec % 2 != 0: return 0.0, 0.0, True, "Radical"
        homo = (n_elec // 2) - 1
        if homo >= len(vals_g)-1: return 0.0, 0.0, False, "Metal"

        E_v_g = vals_g[homo]
        E_c_g = vals_g[homo+1]
        
        # Dispersion Scan
        k_scans = {"X":[0.5,0,0], "Y":[0,0.5,0], "Z":[0,0,0.5]}
        min_gap_tb = E_c_g - E_v_g
        max_disp = 0.0
        axis = "None"
        
        for label, kv in k_scans.items():
            vals_k = np.linalg.eigvalsh(self.build_H(np.array(kv)))
            gap_k = vals_k[homo+1] - vals_k[homo]
            if gap_k < min_gap_tb: min_gap_tb = gap_k
            disp = abs(E_v_g - vals_k[homo])
            if disp > max_disp:
                max_disp = disp
                axis = label
                
        screening = min(0.8, max_disp * 1.5)
        eff_U = self.avg_U * (1.0 - screening)
        final_gap = min_gap_tb + eff_U
        return final_gap, max_disp, False, axis

# ==========================================
#  STREAMLIT UI
# ==========================================

st.title("⚛️ Quantum Bricks")
st.subheader("Deterministic Physics Screener for Organic Semiconductors")

st.markdown("""
This tool uses **Tight-Binding Theory + Hubbard U Correction** to screen organic crystals instantly.
No AI. No black boxes. Just orbital physics.

**How to read results:**
*   **Dispersion > 0.4 eV:** 🏆 **High Mobility** (Speed King / OFET candidate).
*   **Dispersion < 0.1 eV:** 🧱 **Low Mobility** (Insulator / Molecular Crystal).
*   **Gap ~ 1.0 - 2.0 eV:** ⚡ **Semiconductor**.
*   **Gap = 0.0 eV:** ⚠️ **Radical/Metal** (or broken file).
""")

uploaded_files = st.file_uploader("Upload CIF files", accept_multiple_files=True, type="cif")

if uploaded_files:
    st.write(f"Processing {len(uploaded_files)} files...")
    
    results_data = []
    
    progress_bar = st.progress(0)
    
    for i, file in enumerate(uploaded_files):
        try:
            # Parse string content
            stringio = StringIO(file.getvalue().decode("utf-8"))
            cif_content = stringio.read()
            
            lattice, atoms = parse_cif_content(cif_content)
            
            if not atoms:
                results_data.append({
                    "Filename": file.name,
                    "Status": "Error (No Atoms)",
                    "Gap (eV)": 0.0,
                    "Dispersion (eV)": 0.0,
                    "Axis": "-"
                })
                continue
                
            solver = TightBindingSolver(lattice, atoms)
            gap, disp, rad, axis = solver.solve()
            
            status = "Semiconductor"
            if rad:
                status = "Radical/Metal (0 eV)"
                gap = 0.0
            elif disp > 0.4:
                status = "🏆 HIGH MOBILITY"
            elif disp < 0.1:
                status = "Low Mobility"
            
            results_data.append({
                "Filename": file.name,
                "Status": status,
                "Gap (eV)": round(gap, 3),
                "Dispersion (eV)": round(disp, 3),
                "Axis": axis,
                "Atoms": len(atoms)
            })
            
        except Exception as e:
            results_data.append({
                "Filename": file.name,
                "Status": f"Error: {str(e)}",
                "Gap (eV)": 0.0,
                "Dispersion (eV)": 0.0,
                "Axis": "-"
            })
            
        progress_bar.progress((i + 1) / len(uploaded_files))

    # Display Data
    df = pd.DataFrame(results_data)
    
    # Highlight high mobility rows
    def highlight_winners(row):
        if row['Dispersion (eV)'] > 0.4:
            return ['background-color: #d4edda'] * len(row)
        elif row['Gap (eV)'] == 0.0:
            return ['background-color: #f8d7da'] * len(row)
        else:
            return [''] * len(row)

    st.dataframe(df.style.apply(highlight_winners, axis=1), use_container_width=True)
    
    # CSV Download
    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button(
        "Download Results as CSV",
        csv,
        "quantum_bricks_results.csv",
        "text/csv",
        key='download-csv'
    )

else:
    st.info("Drag and drop CIF files to begin.")