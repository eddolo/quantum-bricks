Here is a professional description you can use for the \*\*GitHub README\*\*, the \*\*Streamlit sidebar\*\*, or the email to your friend.



It emphasizes that this is a \*\*Physics Engine\*\*, not an AI model.



---



\### \*\*App Name:\*\* Quantum Bricks

\*\*Tagline:\*\* Deterministic High-Throughput Screener for Organic Semiconductors.



\### \*\*Short Description (For the App Sidebar)\*\*

Quantum Bricks is a physics-based engine that calculates the \*\*Electronic Band Structure\*\* of organic crystals instantly. Unlike AI models that "guess" based on training data, this tool solves the \*\*Tight-Binding Hamiltonian\*\* with \*\*Hubbard U correction\*\* for every single crystal structure.



It is designed to screen thousands of CIF files to identify \*\*High-Mobility Semiconductors\*\* (High Dispersion) and filter out insulators or broken structures.



---



\### \*\*Full Description (For GitHub / Documentation)\*\*



\#### \*\*The Problem\*\*

Searching for new organic semiconductors (for OFETs, OLEDs, and Solar Cells) is a bottleneck.

\*   \*\*DFT (Density Functional Theory):\*\* Accurate but too slow (12+ hours per crystal).

\*   \*\*AI/ML Models:\*\* Fast but unreliable for exotic structures (Black Box).



\#### \*\*The Solution: Quantum Bricks\*\*

Quantum Bricks is a \*\*"White Box" Physics Engine\*\*. It condenses the essential physics of electron transport into a lightweight, vector-optimized Python script.



It constructs the quantum mechanical Hamiltonian using \*\*Harrison’s Tight-Binding Rules\*\* ($1/d^2$ scaling) and applies \*\*Bloch’s Theorem\*\* to scan the Brillouin Zone in 3D (X, Y, Z axes). It solves for:

1\.  \*\*Optical Band Gap:\*\* Corrected for Electron-Electron Repulsion (Hubbard U).

2\.  \*\*Band Dispersion (Bandwidth):\*\* The critical indicator of electron mobility.

3\.  \*\*Anisotropy:\*\* Identifies the specific crystal axis (stacking direction) where conduction is fastest.



\#### \*\*Validation\*\*

The engine has been benchmarked against the "Hall of Fame" of organic electronics:

\*   ✅ \*\*DNTT:\*\* Correctly identified as a "Speed King" (Dispersion > 0.7 eV).

\*   ✅ \*\*Rubrene \& Pentacene:\*\* Correctly matched Fundamental Transport Gaps (~1.1 - 3.0 eV).

\*   ✅ \*\*CuPc:\*\* Correctly identified the Y-axis stacking and high mobility.

\*   ✅ \*\*Broken Files:\*\* Automatically detects and flags radical fragments (0.0 eV).



\#### \*\*How to Interpret Results\*\*

Upload your CIF files. The table will highlight promising candidates:



\*   🏆 \*\*High Mobility (Dispersion > 0.4 eV):\*\* These are your world-record candidates (like DNTT). Electrons move ballistically through the crystal stacks.

\*   ⚡ \*\*Semiconductor (Gap 1.0 - 3.0 eV):\*\* Good candidates for Solar Cells or OLEDs.

\*   🧱 \*\*Low Mobility (Dispersion < 0.1 eV):\*\* Likely insulators or isolated molecular crystals.

\*   ⚠️ \*\*Radical/Metal (Gap 0.0 eV):\*\* Indicates either a metallic charge-transfer salt or a broken/incomplete CIF file.

