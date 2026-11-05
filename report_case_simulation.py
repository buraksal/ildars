"""
report_case_simulation.py
=========================

Runs the existing ray-tracing engine (engine.py / objects/) on the example
room that is used in the lab report and saves one plot.

Example room (walls are infinite lines through the given points):
    receiver R = (0, 0), source S = (0, 13)
    w1 through (-7, 15) and (-11, -8)
    w2 through ( 3, 15) and (  9,  7)
    w3 through ( 9,  7) and (  3, -8)

Run:   python report_case_simulation.py
Output: data/report_figures/fig_simulator_order2.(pdf|png)
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import engine
import objects.analytics as analysis
from engine import Engine
from objects.line import Line
from objects.point import Point
from objects.vector import Vector

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "report_figures")
os.makedirs(OUT, exist_ok=True)

ORDER = 2
ROUNDS = 50000


def wall(a, b):
    return Line(Point(list(map(float, a))), Vector(Point(list(map(float, a))), Point(list(map(float, b)))))


def main():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
    plt.figure(figsize=(5.6, 4.6))

    eng = Engine()
    eng.reciever = Point([0.0, 0.0])
    eng.sender = Point([0.0, 13.0])
    eng.walls = [wall((-7, 15), (-11, -8)), wall((3, 15), (9, 7)), wall((9, 7), (3, -8))]
    eng.radius = 0.0125
    engine.ROUNDS = ROUNDS

    # traces all rays, draws the received ones and the image points they produce
    eng.sound_events(max_order=ORDER, show=False)

    labels = sorted(img[1] for img in eng.images)
    print("received image points:", labels)

    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-37, 36)
    plt.ylim(-12, 36)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.title(f"ray tracer, reflection order <= {ORDER}, {ROUNDS} rays")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig_simulator_order2.pdf"))
    plt.savefig(os.path.join(OUT, "fig_simulator_order2.png"), dpi=200)
    plt.close()

    with open(os.path.join(OUT, "simulator_received.txt"), "w") as f:
        f.write("\n".join(labels) + "\n")

    # equal-distance pairs between the traced image points (analysis module of the code base)
    eng.create_edm()
    lines = []
    for nd in (1, 2):
        eng.distance_pairs = {}
        eng.compare_point_pair_distances(nd=nd)
        groups = {k: v for k, v in eng.distance_pairs.items() if len(v) > 1}
        clear, rest, _ = analysis.split_dictionary_with_flags(groups)
        lines.append(f"rounding to {nd} decimal(s):")
        for key in sorted(groups):
            lines.append(f"  {key:8.2f}: explained={clear[key]}  unexplained={rest[key]}")
    print("\n".join(lines))
    with open(os.path.join(OUT, "simulator_equal_distances.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
