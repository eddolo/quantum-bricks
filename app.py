import os
import time
import threading
from pathlib import Path
from datetime import datetime

import streamlit as st
import numpy as np
import pandas as pd

# Import your existing engines
from quantum_bricks import parse_cif, TightBindingSolver
from alchemist import run_alchemy

# =====================================================
#  STREAMLIT CONFIG
# =====================================================

st.set_page_config(
    page_title="Quantum Bricks + Alchemist",
    page_icon="⚛️",
    layout="wide",
)

# Root for Streamlit-local persistence
STREAMLIT_ROOT = Path("streamlit_data")
UPLOAD_ROOT = STREAMLIT_ROOT / "uploads"
RESULTS_ROOT = STREAMLIT_ROOT / "results"

UPLOAD_ROOT.mkdir(parents=True, exist_ok=True)
RESULTS_ROOT.mkdir(parents=True, exist_ok=True)


# =====================================================
#  SMALL HELPERS
# =====================================================

def save_uploaded_file(uploaded_file, subdir: str) -> Path:
    """Save a Streamlit UploadedFile to disk and return its path."""
    target_dir = UPLOAD_ROOT / subdir
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / uploaded_file.name
    with open(path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return path


def append_history_csv(path: Path, df_new: pd.DataFrame):
    """Append df_new to CSV, creating the file if it doesn't exist."""
    if path.exists():
        df_old = pd.read_csv(path)
        df_all = pd.concat([df_old, df_new], ignore_index=True)
    else:
        df_all = df_new
    df_all.to_csv(path, index=False)


def highlight_qb(row):
    """Highlight rows in Quantum Bricks results."""
    disp = row.get("Dispersion (eV)", 0.0)
    gap = row.get("Gap (eV)", 0.0)
    base = [''] * len(row)
    if gap == 0.0:
        return ['background-color: #f8d7da'] * len(row)  # red-ish for radicals
    if disp > 0.4:
        return ['background-color: #d4edda'] * len(row)  # green high mobility
    if disp < 0.1:
        return ['background-color: #f0f0f0'] * len(row)  # grey low mobility
    return base


def highlight_alchemist(row):
    """Highlight rows in Alchemist results based on best score."""
    score_cols = [c for c in row.index if "Score" in c]
    base = [''] * len(row)
    if not score_cols:
        return base
    scores = [row[c] for c in score_cols if pd.notna(row[c])]
    if not scores:
        return base
    score = max(scores)
    if score >= 0.7:
        return ['background-color: #d4edda'] * len(row)
    elif score >= 0.4:
        return ['background-color: #fff3cd'] * len(row)
    else:
        return base


def find_latest_alchemist_run_dir(cif_path: Path, mode: str):
    """
    Given a CIF path and mode, find the latest alchemist run directory.

    We support both:
      alchemist/<basename>/<mode>_timestamp/
    and older:
      alchemist/Alchemy_<basename>_timestamp/
    """
    base = cif_path.stem
    root = Path("alchemist")

    if not root.exists():
        return None

    candidates = []

    # New-style: alchemist/<basename>/<mode>_timestamp
    new_root = root / base
    if new_root.exists():
        for d in new_root.iterdir():
            if d.is_dir() and d.name.startswith(f"{mode}_"):
                candidates.append(d)

    # Old-style: alchemist/Alchemy_<basename>_timestamp
    for d in root.iterdir():
        if d.is_dir() and d.name.startswith(f"Alchemy_{base}_"):
            candidates.append(d)

    if not candidates:
        return None

    return sorted(candidates)[-1]


def load_alchemist_summary_and_full(run_dir: Path, mode: str):
    """
    Try to load summary CSV and full log CSV from run_dir.
    Robust to slightly different file naming.
    """
    summary_candidates = [
        run_dir / f"summary_{mode}.csv",
        run_dir / "leaderboard_summary.csv",
        run_dir / "summary.csv",
    ]
    full_candidates = [
        run_dir / f"leaderboard_full_log_{mode}.csv",
        run_dir / "leaderboard_full_log.csv",
        run_dir / "full_log.csv",
    ]

    summary_df = None
    full_df = None

    for p in summary_candidates:
        if p.exists():
            summary_df = pd.read_csv(p)
            break

    for p in full_candidates:
        if p.exists():
            full_df = pd.read_csv(p)
            break

    return summary_df, full_df


# =====================================================
#  UI HEADER
# =====================================================

st.title("⚛️ Quantum Bricks + Alchemist")
st.caption("Fast physics-based screening and mutation engine for organic semiconductors (no AI, no black boxes).")

tab_qb, tab_alch = st.tabs(["Quantum Bricks", "Alchemist"])


# =====================================================
#  TAB 1 — QUANTUM BRICKS
# =====================================================

with tab_qb:
    st.subheader("Quantum Bricks — Band Gap & Mobility Screener")

    st.markdown(
        """
        **What this does:**

        - Reads CIF files
        - Builds a tight-binding Hamiltonian with Hubbard U
        - Computes:
            - Optical band gap (eV)
            - Max band dispersion (eV) along X/Y/Z
            - Axis of highest mobility
        - Classifies:
            - 🏆 High Mobility (dispersion > 0.4 eV)
            - Low Mobility (dispersion < 0.1 eV)
            - ⚠️ Radicals/metals (gap = 0.0 eV)
        """
    )

    uploaded_files_qb = st.file_uploader(
        "Upload one or more CIF files for Quantum Bricks",
        accept_multiple_files=True,
        type=["cif"],
        key="qb_uploader",
    )

    col_l, col_r = st.columns([2, 1])
    with col_l:
        run_qb = st.button("Run Quantum Bricks", type="primary")
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
                    results_rows.append(
                        {
                            "RunTimestamp": datetime.now().isoformat(timespec="seconds"),
                            "Filename": uf.name,
                            "Status": "Error (No Atoms)",
                            "Gap (eV)": 0.0,
                            "Dispersion (eV)": 0.0,
                            "Axis": "-",
                            "Atoms": 0,
                        }
                    )
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

                results_rows.append(
                    {
                        "RunTimestamp": datetime.now().isoformat(timespec="seconds"),
                        "Filename": uf.name,
                        "Status": status,
                        "Gap (eV)": round(float(gap), 3),
                        "Dispersion (eV)": round(float(disp), 3),
                        "Axis": axis,
                        "Atoms": len(atoms),
                    }
                )
            except Exception as e:
                results_rows.append(
                    {
                        "RunTimestamp": datetime.now().isoformat(timespec="seconds"),
                        "Filename": uf.name,
                        "Status": f"Error: {e}",
                        "Gap (eV)": 0.0,
                        "Dispersion (eV)": 0.0,
                        "Axis": "-",
                        "Atoms": 0,
                    }
                )

            progress.progress((i + 1) / len(uploaded_files_qb))

        df_qb = pd.DataFrame(results_rows)

        # Persist run to history
        append_history_csv(qb_history_path, df_qb)

        st.markdown("### Current Run Results")
        st.dataframe(
            df_qb.style.apply(highlight_qb, axis=1),
            use_container_width=True,
        )

        csv_bytes = df_qb.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download this run as CSV",
            csv_bytes,
            file_name="quantum_bricks_results_run.csv",
            mime="text/csv",
        )

    st.markdown("### Quantum Bricks History & Analysis")

    if qb_history_path.exists():
        df_hist = pd.read_csv(qb_history_path)

        with st.expander("Show full Quantum Bricks history", expanded=False):
            st.dataframe(
                df_hist.style.apply(highlight_qb, axis=1),
                use_container_width=True,
            )

            csv_hist = df_hist.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download full Quantum Bricks history",
                csv_hist,
                file_name="quantum_bricks_history.csv",
                mime="text/csv",
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
            (df_filt["Dispersion (eV)"] >= min_disp)
            & (df_filt["Gap (eV)"] >= gap_range[0])
            & (df_filt["Gap (eV)"] <= gap_range[1])
        ]
        if only_high_mob:
            df_filt = df_filt[df_filt["Dispersion (eV)"] > 0.4]

        st.dataframe(
            df_filt.style.apply(highlight_qb, axis=1),
            use_container_width=True,
        )

        st.markdown("#### Dispersion vs Gap (history)")
        if not df_hist.empty:
            st.scatter_chart(df_hist, x="Gap (eV)", y="Dispersion (eV)")
    else:
        st.info("No Quantum Bricks history yet. Run a calculation to populate it.", icon="ℹ️")


# =====================================================
#  TAB 2 — ALCHEMIST (LIVE PER-FILE LEADERBOARD)
# =====================================================

with tab_alch:
    st.subheader("Alchemist — Mutation Engine on Top of Quantum Bricks")

    st.markdown(
        """
        **What this does:**

        - Takes a CIF structure
        - Applies allowed atomic substitutions (N→P/As, O→S/Se, C→Si/Ge, H→F/Cl, etc.)
        - Evaluates each mutant using one of four modes:
            - `fast`  → screened only
            - `relax` → cell relaxation
            - `brute` → strain brute-force
            - `both`  → relax + brute
        - Saves all results under `alchemist/`
        - Produces `summary_<mode>.csv` + human-readable TXT summaries
        """
    )

    mode = st.selectbox(
        "Select Alchemist mode",
        options=["fast", "relax", "brute", "both"],
        index=3,
    )

    uploaded_files_alch = st.file_uploader(
        "Upload one or more CIF files for Alchemist",
        accept_multiple_files=True,
        type=["cif"],
        key="alch_uploader",
    )

    run_alch_btn = st.button("Run Alchemist (sequential, live leaderboards)", type="primary")
    alch_history_path = RESULTS_ROOT / "alchemist_history.csv"

    # ------------- LIVE SEQUENTIAL RUN -------------
    if run_alch_btn and uploaded_files_alch:
        progress = st.progress(0.0)
        all_runs_for_history = []

        # Process CIFs one by one (sequential)
        for i, uf in enumerate(uploaded_files_alch):
            cif_name = uf.name
            cif_key = f"alch_{cif_name}_expanded"

            # Default: expanded while running
            if cif_key not in st.session_state:
                st.session_state[cif_key] = True

            with st.expander(f"{cif_name}", expanded=st.session_state[cif_key]):
                st.markdown(f"**Mode:** `{mode}`")

                # Per-file UI placeholders
                status_placeholder = st.empty()
                summary_placeholder = st.empty()
                full_placeholder = st.empty()

                # Save uploaded CIF for this file
                cif_path = save_uploaded_file(uf, "alchemist_streamlit")

                # Run Alchemist in a background thread
                thread = threading.Thread(target=run_alchemy, args=(str(cif_path), mode))
                thread.start()

                # Live refresh every ~10 seconds until done
                while thread.is_alive():
                    status_placeholder.info("Running Alchemist… refreshing leaderboard every ~10 seconds.")
                    run_dir = find_latest_alchemist_run_dir(cif_path, mode)
                    if run_dir is not None:
                        try:
                            summary_df, full_df = load_alchemist_summary_and_full(run_dir, mode)
                            if summary_df is not None:
                                summary_placeholder.dataframe(
                                    summary_df.style.apply(highlight_alchemist, axis=1),
                                    use_container_width=True,
                                )
                            if full_df is not None:
                                full_placeholder.dataframe(
                                    full_df.style.apply(highlight_alchemist, axis=1),
                                    use_container_width=True,
                                )
                        except Exception:
                            # Ignore partial-write errors and retry next refresh
                            pass

                    time.sleep(10)

                # Thread finished: final load & freeze
                thread.join()
                run_dir = find_latest_alchemist_run_dir(cif_path, mode)

                if run_dir is None:
                    status_placeholder.error("Run finished, but output folder could not be located.")
                else:
                    summary_df, full_df = load_alchemist_summary_and_full(run_dir, mode)

                    if summary_df is not None:
                        status_placeholder.success(f"Completed. Results folder: `{run_dir}`")
                        summary_placeholder.dataframe(
                            summary_df.style.apply(highlight_alchemist, axis=1),
                            use_container_width=True,
                        )

                        # Log to history
                        df_to_log = summary_df.copy()
                        df_to_log["CIF"] = cif_name
                        df_to_log["Mode"] = mode
                        df_to_log["RunDir"] = str(run_dir)
                        df_to_log["RunTimestamp"] = datetime.now().isoformat(timespec="seconds")
                        all_runs_for_history.append(df_to_log)

                    if full_df is not None:
                        full_placeholder.dataframe(
                            full_df.style.apply(highlight_alchemist, axis=1),
                            use_container_width=True,
                        )

                # After completion, auto-collapse this expander (but remember if user re-opens later)
                st.session_state[cif_key] = False

            # Update global progress bar (sequential across CIFs)
            progress.progress((i + 1) / len(uploaded_files_alch))

        # After all CIFs processed: persist history and offer CSV download
        if all_runs_for_history:
            df_all_new = pd.concat(all_runs_for_history, ignore_index=True)
            append_history_csv(alch_history_path, df_all_new)

            csv_new = df_all_new.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download this Alchemist batch summary as CSV",
                csv_new,
                file_name=f"alchemist_summary_{mode}_batch.csv",
                mime="text/csv",
            )

    # ------------- HISTORY VIEW -------------
    st.markdown("### Alchemist History & Analysis")

    if alch_history_path.exists():
        df_hist = pd.read_csv(alch_history_path)

        with st.expander("Show full Alchemist history", expanded=False):
            st.dataframe(
                df_hist.style.apply(highlight_alchemist, axis=1),
                use_container_width=True,
            )

            csv_hist = df_hist.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download full Alchemist history",
                csv_hist,
                file_name="alchemist_history.csv",
                mime="text/csv",
            )

        st.markdown("#### Filtered view")
        col1, col2, col3 = st.columns(3)
        with col1:
            modes_available = sorted(df_hist["Mode"].unique()) if "Mode" in df_hist.columns else []
            sel_mode = st.multiselect("Modes", modes_available, default=modes_available)
        with col2:
            score_min = st.slider("Min score", 0.0, 1.0, 0.0, 0.05)
        with col3:
            cif_filter = st.text_input("Filter by CIF name (contains)", "")

        df_filt = df_hist.copy()
        if "Mode" in df_filt.columns and sel_mode:
            df_filt = df_filt[df_filt["Mode"].isin(sel_mode)]

        score_cols = [c for c in df_filt.columns if "Score" in c]
        if score_cols:
            df_filt["BestScore"] = df_filt[score_cols].max(axis=1)
            df_filt = df_filt[df_filt["BestScore"] >= score_min]

        if cif_filter.strip() and "CIF" in df_filt.columns:
            df_filt = df_filt[df_filt["CIF"].str.contains(cif_filter.strip(), case=False)]

        st.dataframe(
            df_filt.style.apply(highlight_alchemist, axis=1),
            use_container_width=True,
        )
    else:
        st.info("No Alchemist history yet. Run a mutation to populate it.", icon="ℹ️")
