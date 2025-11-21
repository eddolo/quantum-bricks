# 📘 Quantum Bricks & Alchemist v41

### *Fast Physics-Based Screening & Mutation Engine for Organic Semiconductors*

---

## 🧱 What Is This?

This repository contains **two tools** that work together:

1. **Quantum Bricks**
   A deterministic physics engine that computes band gaps, band dispersion, and mobility from CIF crystal structures.

2. **Alchemist**
   A mutation engine that swaps atoms (N→P, O→S, C→Si, etc.), relaxes the structure, applies strain, and finds high-performance semiconductor variants.

These tools provide **DFT-inspired predictions** at a fraction of the computational cost.

---

# 🚀 Quick Start

### **Run bandgap prediction on a CIF file**

```bash
python quantum_bricks.py myfile.cif
```

### **Run Alchemist mutation engine**

```bash
python alchemist.py myfile.cif --both
```

### Supported modes:

```
--fast     → screened (no relaxation)
--relax    → geometric relaxation
--brute    → strain brute-force
--both     → relaxation + brute-force
```

---

# 🧱 Quantum Bricks — Physics Engine

Quantum Bricks is a **white-box** physics model:
it **does not guess** — it **calculates**.

It builds a Harrison-style tight-binding Hamiltonian with:

* s and p orbitals
* distance-dependent hopping (~1/d²)
* Hubbard U correction
* symmetry expansion from the CIF
* Bloch phases over X/Y/Z directions
* element-specific atomic parameters

It computes:

### ✔ **Optical Band Gap (eV)**

Corrected using an effective Hubbard term.

### ✔ **Band Dispersion (eV)**

Max valence-band dispersion → mobility proxy.

### ✔ **Mobility Axis**

X / Y / Z direction giving maximum dispersion.

### ✔ **Classification**

* **High Mobility (Speed King)**: disp > 0.25–0.40 eV
* **Medium Mobility**
* **Low Mobility**

Quantum Bricks is ideal for:

* OFET / OLED material design
* Organic crystal screening
* Fast band structure heuristics
* Early-stage semiconductor discovery

---

# ⚗️ Alchemist — Mutation Engine

Alchemist builds on Quantum Bricks.

For each CIF, it:

1. **Identifies allowed atomic substitutions**
   (e.g., N→P/As, O→S/Se, C→Si/Ge)

2. **Generates mutant CIFs**

3. Evaluates them using the mode you select:

---

## ⚡ FAST mode

```
python alchemist.py myfile.cif --fast
```

Computes only the **screened** (unrelaxed) structure.

* Fastest
* Rough, but good for big batches

---

## 🧘 RELAX mode

```
python alchemist.py myfile.cif --relax
```

Performs geometric relaxation of the cell (a, b, c only) and evaluates the relaxed structure.

* More realistic
* Good for scientific results

---

## 🧩 BRUTE mode

```
python alchemist.py myfile.cif --brute
```

Applies ± strain to the unit cell:

* -10% compression to +15% expansion
* Evaluates bandgap + dispersion for each
* Finds the highest-mobility configuration

Good for high-mobility organic semiconductor design.

---

## 🔮 BOTH mode

```
python alchemist.py myfile.cif --both
```

Computes both:

* Relaxed structure
* Brute-force strain exploration

Provides the most complete insight.

---

# 📂 Output Folder Structure

Each run generates a clean, timestamped folder:

```
alchemist/
    <cif_basename>/
        both_2025-11-20_06-50-22/
            N-P_Mutant_relaxed_both.cif
            N-P_Mutant_Exp0_both.cif
            summary_both.csv
            summary_both.txt
            ...
```

Structure:

* **Top folder** = your CIF name
* **Subfolder** = mode + timestamp
* Inside: all generated CIFs + summary reports

You can run multiple times without overwriting anything.

---

# 📊 Summary Files

Each mode generates a:

### ✔ `summary_<mode>.csv`

Machine-readable results.

### ✔ `summary_<mode>.txt`

Human-readable summary (easy to read/share).

---

## CSV formats

### FAST:

```
Mutation,Type,Formula,Disp_Screened,Gap_Screened,Score_Screened
```

### RELAX:

```
Mutation,Type,Formula,Disp_Relaxed,Gap_Relaxed,Score_Relaxed
```

### BRUTE:

```
Mutation,Type,Formula,Disp_Brute,Gap_Brute,Score_Brute,Strain_Brute
```

### BOTH:

```
Mutation,Type,Formula,
Disp_Relaxed,Gap_Relaxed,Score_Relaxed,
Disp_Brute,Gap_Brute,Score_Brute,Strain_Brute
```

Scores are sorted from best → worst.

---

# 🧪 How to Interpret Results

### ✔ **Dispersion (Mobility Proxy)**

* > 0.4 eV → world-class mobility
* 0.2–0.4 eV → good mobility
* <0.1 eV → poor mobility

### ✔ **Bandgap (Semiconductor Window)**

Ideal: **1–3 eV**
Higher → insulator
Lower → metal or radical state

### ✔ **Score (0 to 1)**

Combines:

* Semiconductor gap window
* Mobility (dispersion)
* Strain penalties

Higher score = better candidate.


---


# 🧱 Quantum Bricks — **Inorganic Engine** (v2.3)

A fast, deterministic physics engine for **inorganic crystals**, built with:

* full symmetry expansion
* periodic 3×3×3 bonding
* ionic / covalent descriptors
* calibrated band-gap corrections
* perovskite t-factor and octahedral factor
* dielectric & effective-mass estimation
* Fermi level alignment

It does **not** rely on machine learning — all values come from **geometric + empirical physics** models.

---

## 🔍 What It Computes

For each CIF:

* **Crystal family** (e.g., ABO₃_PEROVSKITE, GAN, CDS, MOS2)
* **Calibrated band gap**
* **Base geometric band gap** (pre-calibration)
* **Band dispersion** (mobility proxy)
* **Coordination type & disorder metric**
* **Nearest-neighbor bond lengths**
* **Electronegativity mismatch (Δχ)**
* **Goldschmidt tolerance factor (perovskites)**
* **Octahedral factor**
* **Dielectric constants (ε∞, ε₀)**
* **Effective masses (mₑ*, mₕ*)**
* **Fermi level offset**
* **Canonical chemical formula**

This allows **rapid screening** of oxide semiconductors, nitrides, chalcogenides, perovskites, and related materials.

---

## 🏷 Supported Inorganic Families (Calibrated)

The inorganic engine currently provides *full, calibrated descriptors* only for:

* **ABO₃ perovskites**
  (SrTiO₃, BaTiO₃, CaTiO₃, KTaO₃, LaAlO₃…)

* **III–V semiconductors**
  (GaN, AlN, InN, GaAs, InP)

* **II–VI semiconductors**
  (CdS, CdSe, ZnO)

* **MX₂ dichalcogenides**
  (MoS₂, WS₂)

Structures outside these families still produce results, but:

✔ Family = `GENERIC`
✔ Confidence = `OUT_OF_DOMAIN`
✔ Dielectric/masses/t-factor may be blank

This indicates: **no calibration data exists**, but the structural analysis is valid.

---

# ⚗️ Alchemist — **Inorganic Mutation Engine** (AC Inorganics)

A lightweight inorganic counterpart to the organic Alchemist.

AC Inorganics performs:

1. **Element swaps** in the CIF
   (e.g., Ti→Zr, O→S/Se/Te, Ga→Al/In, Ba→Sr/Ca…)

2. **Screened evaluations** (no relaxation)

3. **Strain brute-force** sweeps
   from **−8%** compression to **+18%** expansion
   using isotropic scaling of a/b/c

4. **Scoring & ranking** based on:

   * semiconductor gap window
   * dispersion (mobility)
   * strain penalty
   * family-specific bonuses

5. **Automatic summaries**
   `summary_fast.csv`, `summary_brute.csv`, `summary_both.csv`

---

## 🔧 AC Inorganics — Supported Modes

```
fast   → screened (structure as-is)
brute  → strain sweep only
both   → screened + strain sweep
```

Relaxation mode is intentionally **disabled** for inorganics (no organic-cell relaxation here).

---

## 🔎 Limitations of AC Inorganics (v1.0)

To match expectations clearly:

* Only the families listed above have **full calibrated physics**
* Mutants outside these families fall back to **GENERIC**
* No dielectric/masses/t-factor for out-of-domain structures
* Only **isotropic strain** is applied
* No angle/shear strain
* No relaxation (too unreliable for inorganic CIFs)
* CIF must contain readable labels (Ti1, O2, Ba1…)

Despite this, the engine is extremely fast and excellent for **high-throughput screening**.

---

# 📈 Example AC Inorganics Output Folder

```
ac_inorganic/
    BaTiO3/
        both_2025-11-21_03-56-19/
            Ti_to_Zr_strain-2.cif
            O_to_S_strain4.cif
            Ba_to_Sr_screened.cif
            summary_both.csv
            full_log_both.csv
```

CSV includes:

```
file, status, family, Eg, disp, formula, strain_pct, score, variant
```

---

# 🌐 Streamlit Integration

Both **QB Inorganics** and **AC Inorganics** are included in the Streamlit app:

### **QB Inorganics Tab**

* Upload CIFs
* Perform inorganic band-gap analysis
* View full physics tables
* Download per-run CSV + history CSV

### **AC Inorganics Tab**

* Upload CIFs
* Run mutation engine
* Preview results inside the browser
* Download summary + full logs
* Explore strain diagrams
* Compare mutations per structure

Full real-time mutation monitoring is supported.

---

# 🛠 Requirements

* Python 3.10+
* numpy
* scipy
* pandas

No GPU required.

---

# 📄 License

This project is released under the MIT License.

You are free to:

Use it for personal, academic, commercial, or research purposes

Modify the code

Redistribute it

Incorporate it into other projects

As long as you include the MIT license text and copyright notice.