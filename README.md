# Acoustic Room Geometry Reconstruction & Solvability Simulator

This project is based on and extends the ray-tracing simulator originally created by Luca D'Agostino ([lucad8414/ildars-simulator](https://github.com/lucad8414/ildars-simulator)).

It has been customized, updated, and extended to investigate the **Acoustic Inverse Problem**—specifically proving whether an open/semi-open 2D polygonal room geometry (wall positions, source, and receiver) is **uniquely solvable** using the **Image Source Method (ISM)** and **Euclidean Distance Matrix (EDM)** unmixing.

---

## Key Modifications & Additions

Compared to the base simulator, this repository introduces several improvements and dedicated research tooling:

1. **Ray-Tracing Engine Robustness & Geometry Handling**:
   - Fixed ray propagation in open and semi-open geometric systems (`objects/ray.py` and `engine.py`), preventing infinite loops and missing intersection errors when rays escape into unbounded space.
   - Added support for bounded wall segments (`BoundedLine`) in `custom_case_runner.py`.
   - Improved interactive plotting prompts and file export handling.

2. **Solvability & Uniqueness Proof (`solvability_check.py`)**:
   - Implements the numerical check of the solvability argument for the 3-wall, 1-source, 1-receiver laboratory case.
   - Verifies whether the direct distance $p = \|S - R\|$ is the unique common positive root of geometric quadratics.
   - Reconstructs wall positions from perpendicular bisectors and performs consistency verification on higher-order reflections.
   - Exports verification figures and numerical results to `data/report_figures/`.

3. **Parametric Experimentation & Dataset Generation (`generate_dataset.py`)**:
   - Automated parametric sweeps across reflection orders ($k \in \{1, 2, 3, 4\}$), ray counts ($N \in \{5\,000, 20\,000, 50\,000\}$), and distance rounding precisions ($nd \in \{1, 2\}$).
   - Automatically classifies equal-distance point pairs into true Coxeter symmetries (`Clear`) vs. spurious ghost collisions (`Rest`).
   - Produces structured outputs in `data/solvability_dataset_summary.csv` and `data/solvability_dataset_details.json`.

4. **Report Simulation & Plotting (`report_case_simulation.py` & `generate_plots.py`)**:
   - Reproduces the exact acoustic room scenario presented in the lab report.
   - Generates publication-ready figures in PDF and PNG formats.

---

## Repository Structure

```text
ildars-simulator/
├── data/                               # Generated datasets, plots, and report figures
│   ├── plots/                          # Ray-tracing and image source visualization plots
│   ├── report/                         # Reference report documents
│   ├── report_figures/                 # Publication-ready figures and analysis summaries
│   ├── solvability_dataset_summary.csv # Parametric simulation summary table
│   └── solvability_dataset_details.json# Detailed equal-distance pair classifications
├── objects/                            # Core geometric & ray-tracing objects
│   ├── analytics.py                    # Reflection symmetry classification and algebraic reduction
│   ├── dataset.py                      # Dataset manipulation utilities
│   ├── line.py                         # 2D line representations and intersection geometry
│   ├── point.py                        # 2D point operations
│   ├── ray.py                          # Ray propagation and reflection logic (updated)
│   └── vector.py                       # 2D vector mathematics
├── custom_case_runner.py               # Runner with bounded line support for user geometry
├── engine.py                           # Core simulation engine and ray propagation driver
├── generate_dataset.py                 # Multi-order, multi-ray dataset generation pipeline
├── generate_plots.py                   # Plot generation script for simulation runs
├── report_case_simulation.py           # Simulation runner for lab report figures
├── solvability_check.py                # Analytical and numerical solvability verification
├── GEMINI_CONTEXT.md                   # Full mathematical formulation and empirical findings
├── requirements.txt                    # Python dependencies
└── README.md                           # Project documentation (this file)
```

---

## Installation

Ensure you have Python 3.10+ installed. Install the required dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

### 1. Run the Solvability Check
Numerically verifies the unique solvability proof and generates report figures:
```bash
python solvability_check.py
```
Outputs are saved to `data/report_figures/`.

### 2. Simulate the Lab Report Scenario
Simulates the exact example room configuration and plots received rays and virtual image sources:
```bash
python report_case_simulation.py
```

### 3. Generate the Solvability Dataset
Runs the multi-order ray-tracing sweep across ray counts and precisions:
```bash
python generate_dataset.py
```

### 4. Run the Core Simulator
Run the interactive ray-tracing engine directly:
```bash
python engine.py
```

---

## Attribution & Acknowledgments

- Original simulation architecture and ray-tracing concepts by **Luca D'Agostino**: [lucad8414/ildars-simulator](https://github.com/lucad8414/ildars-simulator).
- Customizations, solvability verification, dataset generators, and bug fixes developed for this acoustic inverse problem research project.