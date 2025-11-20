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