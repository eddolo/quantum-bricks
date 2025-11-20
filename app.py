import os
from pathlib import Path
from datetime import datetime

import streamlit as st
import numpy as np
import pandas as pd

from quantum_bricks import parse_cif, TightBindingSolver
from alchemist import run_alchemy

# =====================================================
#  BASIC CONFIG
# =====================================================

st.set_page_config(
    page_title="Quantum Bricks + Alchemist",
    page_icon="⚛️",
    layout="wide"
)

UPLOAD_ROOT = Path("uploads")
RESULTS_ROOT = Path("results")
UPLOAD_ROOT.mkdir(exist_ok=True)
RESULTS_ROOT.mkdir(exist_ok=True)


# =====================================================
#  SMALL HELPERS
# =====================================================

def save_uploaded_file(uploaded_file, subdir: str) -> Path:
    """Save Streamlit UploadedFile to disk and return its path."""
    target_dir = UPLOAD_ROOT / subdir
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / uploaded_file.name
    with open(path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return path


def append_history_csv(path: Path, df_new: pd.DataFrame):
    """Append df_new to CSV, creating it if needed."""
    if path.exists():
        df_old = pd.read_csv(path)
        df_all = pd.concat([df_old, df_new], ignore_index=True)
    else:
        df_all = df_new
    df_all.to_csv(path, index=False)


def highlight_qb(row):
    """Row highlighter for Quantum Bricks results."""
    disp = row.get("Dispersion (eV)", 0.0)
    gap = row.get("Gap (eV)", 0.0)
    base = [''] * len(row)
    if gap == 0.0:
        return ['background-color: #f8d7da'] * len(row)  # red-ish for radicals
    if disp > 0.4:
        return ['background-color: #d4edda'] * len(row)  # green for high mobility
    if disp < 0.1:
        return ['background-color: #f0f0f0'] * len(row)  # grey for low mobility
    return base


def highlight_alchemist(row):
    """Row highlighter for Alchemist results (by Score if present)."""
    score_cols = [c for c in row.index if "Score" in c]
    base = [''] * len(row)
    if not score_cols:
        return base
    score = max(row[c] for c in score_cols if pd.notna(row[c]))
    if score >= 0.7:
        return ['background-color: #d4edda'] * len(row)
    elif score >= 0.4:
        return ['background-color: #fff3cd'] * len(row)
    else:
        return base


def load_alchemist_latest_summary(cif_path: Path, mode: str):
    """
    After run_alchemy(), inspect alchemist/<basename>/,
    find latest <mode>_timestamp folder, and load summary_<mode>.csv.
    """
    base = cif_path.stem
    root_dir = Path("alchemist") / base
    if not root_dir.exists():
        return None, None

    subdirs = [d for d in root_dir.iterdir()
               if d.is_dir() and d.name.startswith(f"{mode}_")]
    if not subdirs:
        return None, None

    latest_dir = sorted(subdirs)[-1]
    csv_path = latest_dir / f"summary_{mode}.csv"
    if not csv_path.exists():
        return latest_dir, None

    df = pd.read_csv(csv_path)
    return latest_dir, df


# =====================================================
#  HEADER
# =====================================================

st.title("⚛️ Quantum Bricks + Alchemist")
st.caption("Fast physics-based screening and mutation engine for organic semiconductors (no AI, no black boxes).")

tab_qb, tab_alch = st.tabs(["Quantum Bricks", "Alchemist"])


# =====================================================
#  TAB 1: QUANTUM BRICKS
# =====================================================

with tab_qb:
    st.subheader("Quantum Bricks — Band Gap & Mobility Screener")

    st.markdown("""
    **What this does:**

    - Reads CIF files
    - Builds a tight-binding Hamiltonian with Hubbard U
    - Computes:
        - Optical band gap (eV)
        - Max band dispersion (eV) across X/Y/Z
        - Axis of highest mobility
    - Classifies:
        - 🏆 High Mobility (dispersion > 0.4 eV)
        - Low Mobility (dispersion < 0.1 eV)
        - Radicals / metals (gap = 0 eV)
    """)

    uploaded_files_qb = st.file_uploader(
        "Upload one or more CIF files",
        accept_multiple_files=True,
        type=["cif"],
        key="qb_uploader"
    )

    col_l, col_r = st.columns([2, 1])
    with col_l:
        run_qb = st.button("Run Quantum Bricks on uploaded CIFs", type="primary")
    with col_r:
        st.info("You can upload multiple CIFs at once.", icon="📂")

    qb_history_path = RESULTS_ROOT / "quantum_bricks_history.csv"

    if run_qb and uploaded_files_qb:
        results_rows = []
        progress = st.progress(0.0)
        for i, uf in enumerate(uploaded_files_qb):
            try:
                cif_path = save_uploaded_file(uf, "quantum_bricks")
                lattice, atoms = parse_cif(str(cif_path))

                if not atoms:
                    results_rows.append({
                        "RunTimestamp": datetime.now().isoformat(timespec="seconds"),
                        "Filename": uf.name,
                        "Status": "Error (No Atoms)",
                        "Gap (eV)": 0.0,
                        "Dispersion (eV)": 0.0,
                        "Axis": "-",
                        "Atoms": 0,
                    })
                    continue

                solver = TightBindingSolver(lattice, atoms)
                gap, disp, radical, axis = solver.solve()

                status = "Semiconductor"
                if radical:
                    status = "Radical/Metal (0 eV)"
                    gap = 0.0
                elif disp > 0.4:
                    status = "🏆 HIGH MOBILITY"
                elif disp < 0.1:
                    status = "Low Mobility"

                results_rows.append({
                    "RunTimestamp": datetime.now().isoformat(timespec="seconds"),
                    "Filename": uf.name,
                    "Status": status,
                    "Gap (eV)": round(float(gap), 3),
                    "Dispersion (eV)": round(float(disp), 3),
                    "Axis": axis,
                    "Atoms": len(atoms),
                })
            except Exception as e:
                results_rows.append({
                    "RunTimestamp": datetime.now().isoformat(timespec="seconds"),
                    "Filename": uf.name,
                    "Status": f"Error: {e}",
                    "Gap (eV)": 0.0,
                    "Dispersion (eV)": 0.0,
                    "Axis": "-",
                    "Atoms": 0,
                })

            progress.progress((i + 1) / len(uploaded_files_qb))

        df_qb = pd.DataFrame(results_rows)

        # Persist to history
        append_history_csv(qb_history_path, df_qb)

        st.markdown("### Current Run Results")
        st.dataframe(
            df_qb.style.apply(highlight_qb, axis=1),
            use_container_width=True
        )

        csv_bytes = df_qb.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download this run as CSV",
            csv_bytes,
            file_name="quantum_bricks_results_run.csv",
            mime="text/csv"
        )

    # HISTORY + SIMPLE ANALYSIS
    st.markdown("### History & Analysis")

    if qb_history_path.exists():
        df_hist = pd.read_csv(qb_history_path)

        with st.expander("Show full Quantum Bricks history", expanded=False):
            st.dataframe(
                df_hist.style.apply(highlight_qb, axis=1),
                use_container_width=True
            )

            csv_hist = df_hist.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download full Quantum Bricks history",
                csv_hist,
                file_name="quantum_bricks_history.csv",
                mime="text/csv"
            )

        st.markdown("#### Filtered view")
        col1, col2, col3 = st.columns(3)
        with col1:
            min_disp = st.slider("Min dispersion (eV)", 0.0, 1.5, 0.0, 0.05)
        with col2:
            gap_range = st.slider("Band gap range (eV)", 0.0, 5.0, (0.0, 5.0), 0.1)
        with col3:
            only_high_mob = st.checkbox("Only high mobility (disp > 0.4 eV)")

        df_filt = df_hist.copy()
        df_filt = df_filt[
            (df_filt["Dispersion (eV)"] >= min_disp) &
            (df_filt["Gap (eV)"] >= gap_range[0]) &
            (df_filt["Gap (eV)"] <= gap_range[1])
        ]
        if only_high_mob:
            df_filt = df_filt[df_filt["Dispersion (eV)"] > 0.4]

        st.dataframe(
            df_filt.style.apply(highlight_qb, axis=1),
            use_container_width=True
        )

        st.markdown("#### Dispersion vs Gap (history)")
        if not df_hist.empty:
            st.scatter_chart(
                df_hist,
                x="Gap (eV)",
                y="Dispersion (eV)",
            )
    else:
        st.info("No Quantum Bricks history yet. Run a calculation to populate it.", icon="ℹ️")


# =====================================================
#  TAB 2: ALCHEMIST
# =====================================================

with tab_alch:
    st.subheader("Alchemist — Mutation Engine on Top of Quantum Bricks")

    st.markdown("""
    **What this does:**

    - Takes a CIF structure
    - Applies allowed atomic substitutions (N→P/As, O→S/Se, C→Si/Ge, H→F/Cl, etc.)
    - Evaluates each mutant using one of four modes:
        - `fast`  → single screened evaluation (no relaxation or strain)
        - `relax` → geometric relaxation of the unit cell
        - `brute` → brute-force strain scan (compression/expansion)
        - `both`  → relaxation + brute (full exploration)
    - Saves all results under `alchemist/<cif_basename>/mode_timestamp/`
    - Produces `summary_<mode>.csv` and human-readable `.txt` summaries
    """)

    mode = st.selectbox(
        "Select Alchemist mode",
        options=["fast", "relax", "brute", "both"],
        index=3
    )

    uploaded_files_alch = st.file_uploader(
        "Upload one or more CIF files for mutation",
        accept_multiple_files=True,
        type=["cif"],
        key="alch_uploader"
    )

    run_alch = st.button("Run Alchemist on uploaded CIFs", type="primary")

    alch_history_path = RESULTS_ROOT / "alchemist_history.csv"

    if run_alch and uploaded_files_alch:
        progress = st.progress(0.0)
        run_rows_for_history = []

        for i, uf in enumerate(uploaded_files_alch):
            st.markdown(f"#### Running on `{uf.name}`")
            with st.spinner(f"Alchemist ({mode}) running for {uf.name}..."):
                try:
                    cif_path = save_uploaded_file(uf, "alchemist_streamlit")
                    # Run the CLI-like engine function
                    run_alchemy(str(cif_path), mode)

                    run_dir, df_summary = load_alchemist_latest_summary(cif_path, mode)
                    if df_summary is None:
                        st.error("Could not find summary CSV for this run.")
                    else:
                        st.success(f"Completed. Results folder: `{run_dir}`")
                        st.dataframe(
                            df_summary.style.apply(highlight_alchemist, axis=1),
                            use_container_width=True
                        )

                        # Log to global history (one row per mutation)
                        df_to_log = df_summary.copy()
                        df_to_log["CIF"] = uf.name
                        df_to_log["Mode"] = mode
                        df_to_log["RunDir"] = str(run_dir)
                        df_to_log["RunTimestamp"] = datetime.now().isoformat(timespec="seconds")
                        run_rows_for_history.append(df_to_log)
                except Exception as e:
                    st.error(f"Error running Alchemist on {uf.name}: {e}")

            progress.progress((i + 1) / len(uploaded_files_alch))

        if run_rows_for_history:
            df_all_new = pd.concat(run_rows_for_history, ignore_index=True)
            append_history_csv(alch_history_path, df_all_new)

            csv_new = df_all_new.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download this Alchemist run summary as CSV",
                csv_new,
                file_name=f"alchemist_summary_{mode}_run.csv",
                mime="text/csv"
            )

    st.markdown("### Alchemist History & Analysis")

    if alch_history_path.exists():
        df_hist = pd.read_csv(alch_history_path)

        with st.expander("Show full Alchemist history", expanded=False):
            st.dataframe(
                df_hist.style.apply(highlight_alchemist, axis=1),
                use_container_width=True
            )

            csv_hist = df_hist.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download full Alchemist history",
                csv_hist,
                file_name="alchemist_history.csv",
                mime="text/csv"
            )

        st.markdown("#### Filtered view")
        col1, col2, col3 = st.columns(3)
        with col1:
            sel_mode = st.multiselect("Modes", sorted(df_hist["Mode"].unique()), default=list(sorted(df_hist["Mode"].unique())))
        with col2:
            score_min = st.slider("Min score", 0.0, 1.0, 0.0, 0.05)
        with col3:
            cif_filter = st.text_input("Filter by CIF name (contains)", "")

        df_filt = df_hist.copy()
        if sel_mode:
            df_filt = df_filt[df_filt["Mode"].isin(sel_mode)]

        # find all score cols to build a unified "BestScore"
        score_cols = [c for c in df_filt.columns if "Score" in c]
        if score_cols:
            df_filt["BestScore"] = df_filt[score_cols].max(axis=1)
            df_filt = df_filt[df_filt["BestScore"] >= score_min]

        if cif_filter.strip():
            df_filt = df_filt[df_filt["CIF"].str.contains(cif_filter.strip(), case=False)]

        st.dataframe(
            df_filt.style.apply(highlight_alchemist, axis=1),
            use_container_width=True
        )

    else:
        st.info("No Alchemist history yet. Run a mutation to populate it.", icon="ℹ️")
