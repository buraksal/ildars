# ILDARS lab: three walls, one first-order reflection

This repository contains the code for my lab report in the course
*Networks and Telematics* at the University of Freiburg
(Chair of Computer Networks and Telematics).

The report studies a two-dimensional ILDARS setup with one source, one
receiver and three reflecting walls. The receiver gets six signals with the
reflection descriptors `()`, `(w1)`, `(w1,w2)`, `(w2,w1)`, `(w1,w3)` and
`(w2,w3)`. The report shows that the source and all three walls can be
reconstructed from the directions of arrival and the time differences of
these signals.

## Origin of the code

The ray-tracing simulator (`engine.py`, `objects/`, `simmer.py`,
`point.ipynb`) was written by **Luca Diringer**
(original repository: https://github.com/lucad8414/ildars-simulator).
I did not write it. I made some changes to it and added the scripts listed below.

## What I added or changed

| File | Purpose |
|------|---------|
| `solvability_check.py` | Main check for the report. Implements the reconstruction method (two quadratic equations for the source distance `p`, walls as perpendicular bisectors), and runs the tests: example room, 5000 random rooms, Jacobian rank, noise study. Only needs `numpy` and `matplotlib`. Independent of the simulator. |
| `report_case_simulation.py` | Runs the ray tracer on the example room of the report and lists the received image points and equal-distance pairs. |
| `custom_case_runner.py`, `generate_plots.py`, `generate_dataset.py` | Earlier scripts for running the simulator on a fixed configuration and producing plots and data sets. |
| `engine.py`, `objects/ray.py` | Small changes to run the simulator on my configurations. |

Note: the data sets from `generate_dataset.py` group image-point pairs by
*rounded* distances from traced rays. These are useful for seeing patterns,
but they are not used as a proof of solvability, because the traced points
are only accurate up to the receiver radius. The solvability argument is in
the report and is checked numerically by `solvability_check.py`.

## How to run

Python 3.10 or newer.

    pip install -r requirements.txt
    python solvability_check.py        # figures and results.json in data/report_figures/
    python report_case_simulation.py   # simulator plot for the example room

## Example room used in the report

Receiver `R = (0, 0)`, source `S = (0, 13)`. The walls are infinite lines through:

- `w1`: (-7, 15) and (-11, -8)
- `w2`: (3, 15) and (9, 7)
- `w3`: (9, 7) and (3, -8)