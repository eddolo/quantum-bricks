# 📘 **README for Streamlit App: Quantum Bricks + Alchemist**

## ⚛️ Quantum Bricks + Alchemist

### *Interactive Physics-Based Semiconductor Screening Platform*

---

# 🧱 Overview

This Streamlit app provides an interactive interface to two deterministic physics engines:

---

## **1. Quantum Bricks**

A **fast tight-binding Hamiltonian solver** that computes semiconductor properties directly from CIF crystal structures:

* Optical **Band Gap (eV)**
* **Band Dispersion (eV)** — mobility proxy
* **Transport Axis** (X/Y/Z)
* **Mobility Classification**
* Automatic detection of **radicals/metals**
* No AI. No machine learning. Pure orbital physics.

---

## **2. Alchemist**

A **mutation engine** built on top of Quantum Bricks.
For each CIF, it performs:

* Allowed atomic substitutions (N→P, O→S, C→Si, H→F, etc.)
* Optional **unit-cell relaxation**
* Optional **brute-force strain scanning**
* Combined full exploration (`both` mode)
* Writes all generated mutant CIFs to organized folders
* Produces detailed summary CSVs & TXT reports

Together, these tools allow users to:

* Screen organic semiconductors quickly
* Discover improved variants via atomic mutations
* Explore strain engineering
* Build persistent analysis histories
* Compare materials over time

---

# 🚀 Features

### 🧩 **Two Separate Tabs**

#### **TAB 1 — Quantum Bricks**

* Upload **one or multiple CIF files**
* Instant calculation of:

  * Band gap
  * Dispersion
  * Mobility axis
* Automatic classification:

  * 🏆 High Mobility
  * Medium Mobility
  * Low Mobility
  * ⚠️ Radical / Metal
* Persistent results tracking
* Filtering, scatter plots, and history browser
* CSV export

---

#### **TAB 2 — Alchemist**

* Upload **one or multiple CIFs**
* Choose mutation mode:

  * `fast` → screened only
  * `relax` → geometric relaxation
  * `brute` → strain brute-force
  * `both` → relax + brute-force
* Generates:

  * Mutant CIFs
  * Summary CSVs
  * Human-readable TXT summaries
* Fully persistent history
* Filtering by:

  * Score
  * CIF name
  * Mode
* CSV export

---

# ⚡ Physics Behind Quantum Bricks

Quantum Bricks implements:

* Harrison’s tight-binding rules: **1/d²** interaction scaling
* s and p orbital basis
* On-site energies (**Es**, **Ep**)
* Hubbard U correction (effective screened U)
* Bloch’s theorem k-space scan (X / Y / Z)
* Unit cell reconstruction directly from CIF
* Symmetry operator application
* Distance-based hopping cutoffs

It computes:

### **Optical Band Gap**

`Gap = min_gap_TB + U_eff`

### **Dispersion**

Max valence-band curvature between Γ and X/Y/Z.

### **Mobility Classification**

* **Disp > 0.4 eV** → 🏆 High mobility
* **0.1–0.4 eV** → medium
* **< 0.1 eV** → insulating/molecular

---

# ⚗️ Alchemist — Mutation Engine Details

For every CIF:

### ✔ Allowed substitutions

```
H → F, Cl
C → Si, Ge
N → P, As
O → S, Se
S → O, Se, Te
Se → S, Te
Si → C, Ge
Te → Se, S
```

### ✔ Modes

#### **FAST**

screened only, no relaxation
→ for large batch scans

#### **RELAX**

optimizes cell (a, b, c) to prevent atomic collisions
→ more realistic predictions

#### **BRUTE**

applies ± strain
→ finds max-performance geometry
→ essential for high-mobility organic crystals

#### **BOTH**

relax + brute
→ full exploration
→ best scientific results

---

# 📂 Output Structure

Every run writes results in:

```
uploads/
results/
alchemist/
    <cif_basename>/
        fast_YYYY-MM-DD_HH-MM-SS/
        relax_YYYY-MM-DD_HH-MM-SS/
        brute_YYYY-MM-DD_HH-MM-SS/
        both_YYYY-MM-DD_HH-MM-SS/
```

Inside each mode folder:

```
summary_<mode>.csv
summary_<mode>.txt
<mutation_cif_files>.cif
```

---

# 📊 Summary CSV Format

### **FAST mode**

```
Mutation, Type, Formula, Disp_Screened, Gap_Screened, Score_Screened
```

### **RELAX mode**

```
Mutation, Type, Formula, Disp_Relaxed, Gap_Relaxed, Score_Relaxed
```

### **BRUTE mode**

```
Mutation, Type, Formula, Disp_Brute, Gap_Brute, Score_Brute, Strain_Brute
```

### **BOTH mode**

```
Mutation, Type, Formula,
Disp_Relaxed, Gap_Relaxed, Score_Relaxed,
Disp_Brute, Gap_Brute, Score_Brute, Strain_Brute
```

---

# 🔍 Interpreting Results

### **Band Gap**

* 1–3 eV → ideal semiconductor window
* <1 eV → radical/metal
* > 3.5 eV → insulator

### **Dispersion**

* > 0.4 eV → high mobility
* 0.1–0.4 eV → moderate
* <0.1 eV → poor

### **Score**

Integrates:

* semiconductor window
* mobility
* strain penalty

---

# 💾 Persistent History

The app saves:

```
results/
    quantum_bricks_history.csv
    alchemist_history.csv
```

This allows:

* revisiting past runs
* comparisons
* scatter plots
* global filtering
* consistent notebooks for research workflows


---


# 🧱 Inorganic Physics Engine (QB Inorganics)**

Quantum Bricks now includes a **second physics engine**, dedicated to **inorganic** CIF structures:

### **QB Inorganics — CIF-Based Band Gap & Physics Analyzer**

This engine performs fast, calibrated analysis for **oxides, nitrides, chalcogenides, III–V, II–VI, and perovskites**.

It computes:

* **Crystal family** (e.g., ABO₃_PEROVSKITE, GAN, CDS, MOS2)
* **Calibrated band gap**
* **Base geometric band gap**
* **Band dispersion**
* **Nearest-neighbor bonding**
* **Electronegativity mismatch (Δχ)**
* **Coordination type**
* **Disorder (RMS angle deviation)**
* **Fermi level offset**
* **Dielectric constants (ε∞, ε₀)**
* **Effective masses (mₑ*, mₕ*)**
* **Tolerance factor / octahedral factor (perovskites)**
* **Canonical formula reconstruction**

All results are generated **deterministically**, with no machine learning.

---

## 🧩 **TAB 3 — QB Inorganics (Inorganic Analyzer)**

This tab allows users to:

* Upload **inorganic CIF files**
* View full per-file physics tables
* Download results as CSV
* Browse local history
* Filter by Eg, dispersion, family, or domain confidence
* Scatter charts (Eg vs dispersion)
* Detect out-of-domain structures

---

### ⚠️ **QB Inorganics Limitations (v1.0)**

Full calibrated physics is available only for:

* **ABO₃ Perovskites**
* **III–V** (GaN, AlN, InN, GaAs, InP)
* **II–VI** (CdS, CdSe, ZnO)
* **MX₂ dichalcogenides** (MoS₂, WS₂)

All other structures fall back to:

* `GENERIC`
* `OUT_OF_DOMAIN`

These entries still show geometric descriptors, but **dielectric constants, effective masses, t-factor, and calibration corrections may be missing**.

The analyzer does **not** perform:

* structural relaxation
* phonons
* DFT
* angle or shear strain
* charged-defect physics

It is intended for **fast screening**, not full ab initio simulation.

---

# ⚗️ **NEW: Inorganic Mutation Engine (AC Inorganics)**

A streamlined variant of Alchemist specialized for **inorganic materials**.

### **AC Inorganics performs:**

* Allowed elemental substitutions (e.g., Ti→Zr, O→S/Se/Te, Ga→In/Al)
* Screened evaluation (`fast`)
* Isotropic strain brute-force (`brute`)
* Combined (`both`)
* Ranking via an inorganic-specific scoring model
* Generation of mutant CIFs
* Writing full logs and summary CSVs

---

## 🧩 **TAB 4 — AC Inorganics (Inorganic Mutation Engine)**

This tab offers:

* Multi-file uploads
* Per-file expanders showing live results
* Screened + strain-sweep evaluations
* Automatic summary CSV downloads
* Full ZIP export for each run
* Persistent mutation history with filters
* Color-coded score interpretation

---

### ⚠️ **AC Inorganics Limitations (v1.0)**

To match the physics model:

* Only the families supported by **QB Inorganics** receive calibrated properties.
* Mutants outside these families revert to `GENERIC` results.
* Only **isotropic strain** is applied (a=b=c scaling).
* No geometric relaxation is performed (unreliable for inorganics).
* No angle / shear strain.
* CIF must contain standard labels (Ti1, O2, Ba1…).

Even with these limits, AC Inorganics is extremely fast and well-suited for **high-throughput inorganic semiconductor mutation screening**.

---

# 📂 **Output: AC Inorganics**

```
ac_inorganic/
    <basename>/
        fast_YYYY-MM-DD_HH-MM-SS/
        brute_YYYY-MM-DD_HH-MM-SS/
        both_YYYY-MM-DD_HH-MM-SS/
            full_log_<mode>.csv
            summary_<mode>.csv
            <mutant>.cif
```

Summaries include:

```
file, family, Eg, disp, dielectric, masses,
formula, confidence, mutation, strain_pct, score, variant
```

---

# 📌 Summary

You now have **four major engines** in one Streamlit app:

1. **Quantum Bricks (Organics)**
2. **Alchemist (Organics)**
3. **QB Inorganics (Inorganic Analyzer)**
4. **AC Inorganics (Inorganic Mutation Engine)**

Each is deterministic, fast, reproducible, and designed for **scientific-grade semiconductor discovery**.

---

# 🛠 Installation

### **Install dependencies**

```
pip install -r requirements.txt
```

### **Run the app**

```
streamlit run app.py
```

Streamlit will open in your browser.

---

# 📜 License

This project is distributed under the **MIT License**.