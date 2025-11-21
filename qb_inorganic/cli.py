from pathlib import Path
import os

from .analyze import analyze_inorganic
from .constants import ENGINE_VERSION


def save_summary_txt(res):
    """Write a summary .txt next to the CIF."""
    out_path = Path(res["file"]).with_suffix(".txt")
    lines = []
    lines.append(f"File: {res['file']}")
    lines.append(f"Formula (simple):    {res.get('formula','?')}")
    lines.append(f"Formula (canonical): {res.get('formula_canonical','?')}")
    lines.append(f"Family: {res['family']} ({res['confidence']})")
    lines.append(f"Base Gap: {res['Eg_base']:.3f} eV")
    lines.append(f"Corrected Gap: {res['Eg']:.3f} eV")
    lines.append(f"Status: {res['status']}")
    lines.append("")
    lines.append(f"N_atoms: {res['N_atoms']}")
    lines.append(f"N_bonds: {res['N_bonds']}")
    if res["d_nn"] is not None:
        lines.append(f"Mean NN bond: {res['d_nn']:.3f} Å")
    if res["avg_delta_ch"] is not None:
        lines.append(f"Δχ avg: {res['avg_delta_ch']:.3f}")
    if res["coord_type"]:
        lines.append(f"Coordination: {res['coord_type']}")
    if res["rms_angle"] is not None:
        lines.append(f"RMS angle dev: {res['rms_angle']:.2f}°")
    if res.get("t_factor") is not None:
        lines.append(f"Goldschmidt t: {res['t_factor']:.3f}")
    if res.get("oct_factor") is not None:
        lines.append(f"Octahedral μ: {res['oct_factor']:.3f}")
    if res.get("E_F") is not None:
        lines.append(f"Fermi offset (rel.): {res['E_F']:.3f} (± mid-gap)")
    if res.get("eps_inf") is not None and res.get("eps_static") is not None:
        lines.append(f"ε∞: {res['eps_inf']:.2f}   ε₀: {res['eps_static']:.2f}")
    if res.get("m_e") is not None and res.get("m_h") is not None:
        lines.append(f"m_e*: {res['m_e']:.3f} m0   m_h*: {res['m_h']:.3f} m0")

    out_path.write_text("\n".join(lines))


def run_on_path(p: str):
    from .analyze import analyze_inorganic

    print(f"\n--- INORGANIC ENGINE v{ENGINE_VERSION} ---")
    print(f"File: {p}")
    try:
        res = analyze_inorganic(p)
        if res.get("N_atoms", 0) == 0:
            print("No atoms parsed (or unsupported format).")
            return

        print(f"  Atoms (after symmetry): {res['N_atoms']}")
        print(f"  Formula (simple):       {res['formula']}")
        print(f"  Formula (canonical):    {res['formula_canonical']}")
        if res.get("family"):
            print(f"  Family prototype:       {res['family']}")
            print(f"  Domain confidence:      {res['confidence']}")
        if res["d_nn"] is not None:
            print(f"  Mean nearest-neighbour bond length: {res['d_nn']:.3f} Å")
        if res["avg_delta_ch"] is not None:
            print(f"  Avg cation–anion Δχ:    {res['avg_delta_ch']:.3f}")
        if res["coord_type"] is not None:
            print(f"  Coordination:           {res['coord_type']}")
        if res["rms_angle"] is not None:
            print(f"  RMS angle deviation:    {res['rms_angle']:.2f}°")
        if res.get("t_factor") is not None and res.get("oct_factor") is not None:
            print(f"  Goldschmidt t-factor:   {res['t_factor']:.3f}")
            print(f"  Octahedral factor μ:    {res['oct_factor']:.3f}")

        if res.get("E_F") is not None:
            print(f"  Fermi level offset:     {res['E_F']:.3f} (relative to mid-gap)")
        if res.get("eps_inf") is not None and res.get("eps_static") is not None:
            print(f"  ε∞, ε₀:                  {res['eps_inf']:.2f}, {res['eps_static']:.2f}")
        if res.get("m_e") is not None and res.get("m_h") is not None:
            print(f"  m_e*, m_h*:             {res['m_e']:.3f} m0, {res['m_h']:.3f} m0")

        print(f"> Base Band Gap (geom):   {res['Eg_base']:.3f} eV")
        print(f"> Corrected Band Gap:     {res['Eg']:.3f} eV")
        print(f"> 'Dispersion':           {res['disp']:.3f} (1/Å, heuristic)")
        print(f"> Status:                 {res['status']}")
        save_summary_txt(res)

    except Exception as e:
        print(f"ERR: {e}")


def main():
    import sys

    args = sys.argv[1:]
    if not args:
        print("Usage: python quantum_bricks_inorganic.py <file_or_folder> [...]")
        sys.exit(1)

    for arg in args:
        p = Path(arg)
        if p.is_dir():
            for f in sorted(p.iterdir()):
                if f.suffix.lower() == ".cif":
                    run_on_path(str(f))
        else:
            run_on_path(str(p))
