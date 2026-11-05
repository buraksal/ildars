import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt

# Add current directory to sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from objects.point import Point
from objects.vector import Vector
from objects.line import Line
from objects.ray import Ray
from engine import Engine

class BoundedLine(Line):
    """
    Line segment bounded strictly between p_start and p_end.
    """
    def __init__(self, p_start: Point, p_end: Point):
        anchor = p_start
        direction = Vector(p_start, p_end)
        super().__init__(anchor, direction)
        self.p_start = p_start
        self.p_end = p_end

    def contains_point(self, p: Point) -> bool:
        min_x = min(self.p_start.value[0], self.p_end.value[0]) - 1e-3
        max_x = max(self.p_start.value[0], self.p_end.value[0]) + 1e-3
        min_y = min(self.p_start.value[1], self.p_end.value[1]) - 1e-3
        max_y = max(self.p_start.value[1], self.p_end.value[1]) + 1e-3
        return (min_x <= p.value[0] <= max_x) and (min_y <= p.value[1] <= max_y)

    def intersection(self, other, show: bool = False):
        inter = super().intersection(other, show)
        if inter is not None and self.contains_point(inter):
            return inter
        return None

    def perp_point(self, p: Point, show: bool = False) -> Point:
        """
        Perpendicular projection onto the infinite mirror line for virtual image source math.
        """
        r = [self.direction.value[0], self.direction.value[1]]
        v1 = Vector(Point([r[1], 2 * r[0]]), Point([2 * r[1], r[0]]))
        v2 = Vector(Point([2 * r[1], r[0]]), Point([r[1], 2 * r[0]]))
        l1 = Line(p, v1)
        l2 = Line(p, v2)
        inter1 = super().intersection(l1)
        inter2 = super().intersection(l2)

        if inter1 is None:
            return inter2
        else:
            return inter1

# Patch Ray.expand to safely handle open-space ray escaping
def open_system_ray_expand(self, walls: list[Line], reciever: tuple[Point, float], show: bool = False) -> list[Vector]:
    if self.expanded:
        return self.values
    self.expanded = True

    prev = None
    block = -1
    self.value = self.value.extend()

    for o in range(0, self.order + 1):
        curr = self.value

        intersections = []
        for j, w in enumerate(walls):
            if j != block:
                intersections.append(w.intersection(curr))
            else:
                intersections.append(None)
        
        intersections.append(curr.intersection(reciever))

        if intersections[-1] is not None:
            self.recieved = True

        vecs = []
        dist = []
        for i, inter in enumerate(intersections):
            if inter is not None:
                vecs.append(Vector(curr.anchor, inter))
                dist.append(abs(vecs[i]))
                if (prev is not None) and abs(vecs[i]) <= abs(Vector(curr.anchor, prev)):
                    dist[i] = math.inf
            else:
                vecs.append(None)
                dist.append(math.inf)
        
        min_d = min(dist)
        if min_d == math.inf:
            # Ray escaped to open space without hitting any wall or receiver
            break

        m = dist.index(min_d)
        block = m
        self.logger.append(m + 1)

        if self.recieved:
            int_point = intersections[-1]
            new_line = None
        else:
            int_point = walls[m].intersection(curr)
            if int_point is None:
                break
            new_line = walls[m].image(int_point, curr.anchor, show, o)

        if prev is None:
            self.values.append(Vector(curr.anchor, int_point))
        else:
            self.values.append(Vector(prev, int_point))

        prev = int_point
        self.value = new_line

        if self.recieved:
            break

    self.expanded = True
    return self.values

# Apply open system ray patch
Ray.expand = open_system_ray_expand

def run_exact_coordinates_simulation():
    print("=" * 60)
    print("RUNNING SIMULATION WITH EXACT USER COORDINATES (2nd Order)")
    print("=" * 60)

    # 1. Receiver at (0, 0)
    receiver = Point([0.0, 0.0])

    # 2. Source at (-1, 11)
    source   = Point([-1.0, 11.0])

    # 3. Left Wall (w1): (-3, 17) to (-9, -9)
    w1 = BoundedLine(Point([-3.0, 17.0]), Point([-9.0, -9.0]))

    # 4. Right Top Wall (w2): (3, 17) to (9, -9)
    w2 = BoundedLine(Point([3.0, 17.0]), Point([9.0, -9.0]))

    # 5. Right Bottom Wall (w3): (6, 5.5) to (4, -9)
    w3 = BoundedLine(Point([6.0, 5.5]), Point([4.0, -9.0]))

    walls = [w1, w2, w3]

    print("Elements configured:")
    print(f"  Receiver: {receiver}")
    print(f"  Source:   {source}")
    print(f"  Wall 1 (Left):       {w1.p_start} -> {w1.p_end}")
    print(f"  Wall 2 (Right Top):  {w2.p_start} -> {w2.p_end}")
    print(f"  Wall 3 (Right Bot):  {w3.p_start} -> {w3.p_end}")

    # Initialize Engine
    eng = Engine()
    eng.walls = walls
    eng.sender = source
    eng.reciever = receiver
    eng.radius = 0.5  # Capture tolerance radius

    # Run ray-tracing simulation up to 2nd-order reflections
    reflection_order = 2
    print("\nSimulating sound rays...")
    eng.sound_events(max_order=reflection_order, show=False)
    eng.create_edm()

    print("\n--- RESULTS SUMMARY ---")
    print(f"Total Image Points Captured: {len(eng.images)}")
    for img in eng.images:
        print(f"  Image Source Label: {img[1]:<10} Position: {img[0]}")

    # Plot simulation
    plt.figure(figsize=(10, 9))
    
    # Plot walls
    colors = ["blue", "green", "magenta"]
    names = ["Wall 1 (Left)", "Wall 2 (Right Top)", "Wall 3 (Right Bot)"]
    for i, w in enumerate(walls):
        plt.plot([w.p_start.value[0], w.p_end.value[0]], [w.p_start.value[1], w.p_end.value[1]], 
                 color=colors[i], label=names[i], linewidth=2.5)
    
    # Plot sender & receiver
    source.plot(color="red", name="Source (-1,11)")
    receiver.plot(color="cyan", name="Receiver (0,0)")

    plt.autoscale(False)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Exact User Coordinates Simulation (2nd Order Reflections)")
    plt.legend(loc="upper left")
    plt.grid(True)
    
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(data_dir, exist_ok=True)
    filename = os.path.join(data_dir, "exact_coords_plot.png")
    plt.savefig(filename, bbox_inches='tight')
    print(f"\nSaved plot image to: {filename}")
    print("=" * 60)

if __name__ == "__main__":
    run_exact_coordinates_simulation()
