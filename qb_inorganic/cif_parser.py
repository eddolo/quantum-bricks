import re
import numpy as np
from .constants import COVALENT_RADII

def apply_symmetry(atoms, ops):
    """Apply symmetry operations to asymmetric atoms, deduplicating by frac position."""
    if not ops:
        return atoms

    new_atoms = []
    generated = []

    for op in ops:
        op_clean = op.lower().replace(" ", "")
        parts = op_clean.split(",")
        if len(parts) != 3:
            continue

        for atom in atoms:
            x, y, z = atom["frac"]
            ctx = {"x": x, "y": y, "z": z}
            try:
                nx = eval(parts[0], {"__builtins__": None}, ctx)
                ny = eval(parts[1], {"__builtins__": None}, ctx)
                nz = eval(parts[2], {"__builtins__": None}, ctx)
            except Exception:
                continue

            f = np.array([nx % 1.0, ny % 1.0, nz % 1.0])
            if not any(np.linalg.norm(f - g) < 1e-3 for g in generated):
                generated.append(f)
                new_atoms.append({"elem": atom["elem"], "frac": f})

    return new_atoms or atoms


def parse_atoms_with_sym(text: str):
    """
    Parse atoms and symmetry from CIF.

    Supports:
    - _symmetry_equiv_pos_as_xyz
    - _space_group_symop_operation_xyz
    """
    lines = text.splitlines()
    atoms = []
    sym_ops = []

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

        # symmetry blocks (two common keys)
        if s.startswith("_symmetry_equiv_pos_as_xyz") or s.startswith("_space_group_symop_operation_xyz"):
            in_sym = True
            in_atoms = False
            continue

        if s.startswith("_atom_site_"):
            in_atoms = True
            in_sym = False
            idxs[s.split()[0]] = len(idxs)
            continue

        if in_sym:
            if s.startswith("_"):
                continue
            op = s.replace("'", "").replace('"', "")
            if op:
                sym_ops.append(op)

        elif in_atoms:
            if s.startswith("_"):
                continue

            cols = s.split()
            if not cols:
                continue

            def col(prefix):
                for key, idx in idxs.items():
                    if key.startswith(prefix):
                        return idx
                return None

            lab_c = col("_atom_site_label")
            if lab_c is None:
                lab_c = col("_atom_site_type_symbol")
            x_c = col("_atom_site_fract_x")
            y_c = col("_atom_site_fract_y")
            z_c = col("_atom_site_fract_z")

            if None in (lab_c, x_c, y_c, z_c):
                continue
            if len(cols) <= max(lab_c, x_c, y_c, z_c):
                continue

            label = cols[lab_c]
            elem = "".join(ch for ch in label if ch.isalpha()) or label

            if elem not in COVALENT_RADII:
                elem = elem[:2]
            if elem not in COVALENT_RADII:
                elem = elem[:1]
            if elem not in COVALENT_RADII:
                continue

            def clean(v):
                return float(re.sub(r"\([^)]*\)", "", v))

            try:
                fx = clean(cols[x_c])
                fy = clean(cols[y_c])
                fz = clean(cols[z_c])
            except Exception:
                continue

            atoms.append({"elem": elem, "frac": np.array([fx, fy, fz])})

    full_atoms = apply_symmetry(atoms, sym_ops)
    return full_atoms
