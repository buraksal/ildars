import os
import sys
import csv
import json
import time

# Workspace imports without touching existing source code
from objects.point import Point
from objects.vector import Vector
from objects.line import Line
import engine
from engine import Engine
import objects.analytics as analysis


def setup_user_case(eng: Engine):
    """
    Configures the user's specific acoustic setup:
    Receiver, Sender and 3 Walls as defined in the user's case.
    """
    eng.reciever = Point([0.0, 0.0])
    eng.sender = Point([-1.0, 11.0])

    w1 = Line(Point([-3.0, 17.0]), Vector(Point([-3.0, 17.0]), Point([-9.0, -9.0])))
    w2 = Line(Point([3.0, 17.0]), Vector(Point([3.0, 17.0]), Point([9.0, -9.0])))
    w3 = Line(Point([6.0, 5.5]), Vector(Point([6.0, 5.5]), Point([4.0, -9.0])))

    eng.walls = [w1, w2, w3]
    eng.radius = 0.0125  # Standard detection buffer


def run_dataset_generation():
    print("=" * 80)
    print("AKUSTİK REKONSTRÜKSİYON / RAPOR DATASET ÜRETİCİSİ")
    print("=" * 80)
    print("Hedef: Yansıma dereceleri (Orders) ve Sinyal sayıları (Rounds) taranarak")
    print("sistemin tekil çözüm (uniqueness / solvability) analiz veri seti üretiliyor...\n")

    orders = [1, 2, 3, 4]
    rounds_list = [5000, 20000, 50000]
    precisions = [1, 2]

    results_summary = []
    detailed_dataset = {}

    total_runs = len(orders) * len(rounds_list)
    current_run = 0

    start_total_time = time.time()

    for order in orders:
        for r_count in rounds_list:
            current_run += 1
            print(f"[{current_run}/{total_runs}] Simüle ediliyor -> Reflection Order: {order}, Rounds: {r_count}...", end="", flush=True)
            
            t0 = time.time()
            # Engine instance with user's exact room setup
            eng = Engine()
            setup_user_case(eng)

            # Override global simulation signal rounds
            engine.ROUNDS = r_count

            # Run ray-tracing without displaying plot popups
            eng.sound_events(max_order=order, show=False)
            eng.create_edm()
            sim_time = time.time() - t0
            print(f" (Süre: {sim_time:.2f} sn, Görüntü Noktası: {len(eng.images)})")

            image_labels = [img[1] for img in eng.images]

            case_key = f"order_{order}_rounds_{r_count}"
            detailed_dataset[case_key] = {
                "reflection_order": order,
                "rounds": r_count,
                "captured_images_count": len(eng.images),
                "image_labels": image_labels,
                "precisions": {}
            }

            for nd in precisions:
                # Reset distance pairs for independent calculation per precision
                eng.distance_pairs = {}
                eng.compare_point_pair_distances(nd=nd)

                # Filter pairs having more than 1 occurrence (equal-distance collisions)
                multi_pairs = {k: v for k, v in eng.distance_pairs.items() if len(v) > 1}
                clear, rest, logs = analysis.split_dictionary_with_flags(multi_pairs)

                clear_count = sum(len(v) for v in clear.values())
                rest_count = sum(len(v) for v in rest.values())

                # Uniqueness / Solvability Criterion:
                # Rest = 0 means all equal distance pairs can be explained by algebraic reflection symmetries.
                is_unique = (rest_count == 0)
                status = "Unique (Solvable)" if is_unique else f"Ambiguous ({rest_count} ghost pairs)"

                results_summary.append({
                    "Order (k)": order,
                    "Rounds (N)": r_count,
                    "Precision (nd)": nd,
                    "Captured Images": len(eng.images),
                    "Equal-Dist Keys": len(multi_pairs),
                    "Clear (Algebraic)": clear_count,
                    "Rest (Ghost Pairs)": rest_count,
                    "Solvability Status": status,
                    "Elapsed Time (s)": round(sim_time, 2)
                })

                detailed_dataset[case_key]["precisions"][f"nd_{nd}"] = {
                    "multi_pairs": {str(k): [list(pair) for pair in v] for k, v in multi_pairs.items()},
                    "clear": {str(k): [list(pair) for pair in v] for k, v in clear.items()},
                    "rest": {str(k): [list(pair) for pair in v] for k, v in rest.items()},
                    "is_unique": is_unique
                }

    total_duration = time.time() - start_total_time

    # Output file paths
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(data_dir, exist_ok=True)

    csv_path = os.path.join(data_dir, "solvability_dataset_summary.csv")
    json_path = os.path.join(data_dir, "solvability_dataset_details.json")

    # 1. Save Summary CSV
    fieldnames = list(results_summary[0].keys())
    with open(csv_path, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results_summary)

    # 2. Save Detailed JSON
    with open(json_path, mode="w", encoding="utf-8") as f:
        json.dump(detailed_dataset, f, indent=2, ensure_ascii=False)

    print("\n" + "=" * 80)
    print("DATASET ÜRETİMİ TAMAMLANDI")
    print(f"Toplam Süre: {total_duration:.2f} saniye")
    print(f"1. CSV Özeti Kaydedildi:  {csv_path}")
    print(f"2. Detaylı JSON Kaydedildi: {json_path}")
    print("=" * 80)

    # Print Table to Terminal
    header_str = f"{'Order':<7} | {'Rounds':<8} | {'nd':<4} | {'Images':<8} | {'Multi-Keys':<12} | {'Clear':<8} | {'Rest':<6} | {'Status':<25}"
    print("\n" + header_str)
    print("-" * len(header_str))
    for row in results_summary:
        print(f"{row['Order (k)']:<7} | {row['Rounds (N)']:<8} | {row['Precision (nd)']:<4} | {row['Captured Images']:<8} | {row['Equal-Dist Keys']:<12} | {row['Clear (Algebraic)']:<8} | {row['Rest (Ghost Pairs)']:<6} | {row['Solvability Status']:<25}")
    print("-" * len(header_str))


if __name__ == "__main__":
    run_dataset_generation()
