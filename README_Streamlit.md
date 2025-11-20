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