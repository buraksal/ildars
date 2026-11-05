import os
import sys
import matplotlib
matplotlib.use('Agg')  # Headless mode for saving plots without popup windows
import matplotlib.pyplot as plt

from objects.point import Point
from objects.vector import Vector
from objects.line import Line
import engine
from engine import Engine


def setup_user_case(eng: Engine):
    eng.reciever = Point([0.0, 0.0])
    eng.sender = Point([-1.0, 11.0])

    w1 = Line(Point([-3.0, 17.0]), Vector(Point([-3.0, 17.0]), Point([-9.0, -9.0])))
    w2 = Line(Point([3.0, 17.0]), Vector(Point([3.0, 17.0]), Point([9.0, -9.0])))
    w3 = Line(Point([6.0, 5.5]), Vector(Point([6.0, 5.5]), Point([4.0, -9.0])))

    eng.walls = [w1, w2, w3]
    eng.radius = 0.0125


def generate_plots():
    print("=" * 80)
    print("GRAFİK ÜRETİCİ BAŞLATILDI (PNG OLUŞTURULUYOR)")
    print("=" * 80)

    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "plots")
    os.makedirs(output_dir, exist_ok=True)

    orders = [1, 2, 3, 4]
    rounds_list = [5000, 20000, 50000]

    for order in orders:
        for r_count in rounds_list:
            print(f"Çiziliyor -> Order: {order}, Rounds: {r_count}...", end="", flush=True)
            
            plt.figure(figsize=(10, 10))

            eng = Engine()
            setup_user_case(eng)
            engine.ROUNDS = r_count

            # Simulate with show=True so eng.sound_events plots all received rays, walls, sender & receiver
            eng.sound_events(max_order=order, show=True)

            # Highlight and label image sources
            for img in eng.images:
                pt, label = img
                plt.scatter(pt.value[0], pt.value[1], color="orange", s=65, zorder=5)
                plt.annotate(
                    label,
                    (pt.value[0], pt.value[1]),
                    textcoords="offset points",
                    xytext=(6, 6),
                    fontsize=8,
                    fontweight="bold",
                    color="darkred"
                )

            plt.autoscale(False)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.title(
                f"Acoustic Simulation (Order: {order}, Rounds: {r_count})\n"
                f"Captured Image Sources: {len(eng.images)} | Receiver Buffer: {eng.radius}",
                fontsize=11,
                fontweight="bold"
            )
            plt.xlabel("X Coordinate")
            plt.ylabel("Y Coordinate")
            plt.grid(True, linestyle="--", alpha=0.5)

            filename = f"sim_order_{order}_rounds_{r_count}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=200, bbox_inches='tight')
            plt.close()

            print(f" -> Kaydedildi: {filename}")

    print("\n" + "=" * 80)
    print(f"Tüm grafikler başarıyla kaydedildi! Klasör: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    generate_plots()
