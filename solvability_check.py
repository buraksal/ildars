"""
solvability_check.py
====================

Numerical check of the solvability argument for the ILDARS lab case

    3 walls, 1 source, signals  (), (w1), (w1,w2), (w2,w1), (w1,w3), (w2,w3)

Model (same as in the report):
  * receiver R at the origin, walls are infinite lines  n . x = d  (|n| = 1)
  * signal speed c = 1, so time differences are lengths
  * every signal gives a unit direction v_i and a TDOA delta_i (w.r.t. the
    direct signal), and we know which walls it hit (oracle setting)
  * the image point of signal i is  I_i = (p + delta_i) * v_i,  p = |S - R|

Steps that are checked here:
  1. p is a common positive root of two quadratics
        (A)  |I0 - I12|^2 - |I0 - I21|^2 = 0
        (B)  |I13 - I23|^2 - |I0 - I12|^2 = 0
  2. walls from perpendicular bisectors
        w1 = bis(I0, I1),  w2 = bis(I1, I12),  w3 = bis(I1, I13)
  3. remaining signals (w2,w1), (w2,w3) are used as a consistency check

Only numpy and matplotlib are needed.  Run:   python solvability_check.py
Figures and a JSON file with all numbers are written to  data/report_figures/
"""

import json
import os
from itertools import product

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "report_figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 8,
})

LABELS = [(), (1,), (1, 2), (2, 1), (1, 3), (2, 3)]
NAMES = ["()", "(w1)", "(w1,w2)", "(w2,w1)", "(w1,w3)", "(w2,w3)"]
TEX = [r"I_0", r"I_1", r"I_{12}", r"I_{21}", r"I_{13}", r"I_{23}"]
I0, I1, I12, I21, I13, I23 = range(6)
R = np.zeros(2)


# ----------------------------------------------------------------------------
# basic geometry
# ----------------------------------------------------------------------------
def line_from_points(a, b):
    """wall through the points a and b as (n, d) with n . x = d"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    t = (b - a) / np.linalg.norm(b - a)
    n = np.array([t[1], -t[0]])
    return n, float(n @ a)


def reflect(x, wall):
    n, d = wall
    return x - 2.0 * (n @ x - d) * n


def bisector(a, b):
    """perpendicular bisector of the points a and b as (n, d)"""
    n = (b - a) / np.linalg.norm(b - a)
    return n, float(n @ (a + b) / 2.0)


def image_points(S, walls, labels=LABELS):
    """I_(k1,...,km) = s_km( ... s_k1(S))"""
    out = []
    for lab in labels:
        x = np.array(S, float)
        for k in lab:
            x = reflect(x, walls[k - 1])
        out.append(x)
    return out


def measure(S, walls, labels=LABELS):
    """what the receiver sees: unit directions v_i and TDOAs delta_i"""
    imgs = image_points(S, walls, labels)
    dist = np.array([np.linalg.norm(x) for x in imgs])
    v = np.array([x / np.linalg.norm(x) for x in imgs])
    return v, dist - dist[0]


def wall_error(w_est, w_true):
    """angle between the two lines (deg) and difference of their distances to R"""
    (n1, d1), (n2, d2) = w_est, w_true
    if n1 @ n2 < 0:
        n1, d1 = -n1, -d1
    ang = np.degrees(np.arccos(np.clip(n1 @ n2, -1.0, 1.0)))
    return ang, abs(d1 - d2)


# ----------------------------------------------------------------------------
# the method from the report
# ----------------------------------------------------------------------------
def quad_coeffs(v, dl, a, b, c, d):
    """coefficients [x2, x1, x0] of |I_a - I_b|^2 - |I_c - I_d|^2 as a polynomial in p"""
    def sq(i, j):
        u = v[i] - v[j]
        w = dl[i] * v[i] - dl[j] * v[j]
        return np.array([u @ u, 2.0 * (u @ w), w @ w])
    return sq(a, b) - sq(c, d)


def positive_roots(coef, p_min):
    if abs(coef[0]) < 1e-12:                      # degenerates to a linear equation
        r = np.array([-coef[2] / coef[1]])
    else:
        r = np.roots(coef)
    r = r[np.abs(r.imag) < 1e-9].real
    return sorted(float(x) for x in r if x > p_min)


def reconstruct(p, v, dl):
    P = [(p + dl[i]) * v[i] for i in range(len(v))]
    S = P[I0]
    w1 = bisector(P[I0], P[I1])
    w2 = bisector(P[I1], P[I12])
    w3 = bisector(P[I1], P[I13])
    return S, [w1, w2, w3], P


def consistency_residual(p, v, dl):
    """how far the predicted I21 and I23 are from the measured ones"""
    S, (w1, w2, w3), P = reconstruct(p, v, dl)
    I2 = reflect(S, w2)
    r21 = reflect(I2, w1) - P[I21]
    r23 = reflect(I2, w3) - P[I23]
    return float(np.linalg.norm(np.concatenate([r21, r23])))


def solve(v, dl, tol=1e-6):
    p_min = max(0.0, -float(dl.min()))
    A = quad_coeffs(v, dl, I0, I12, I0, I21)
    B = quad_coeffs(v, dl, I13, I23, I0, I12)
    rA, rB = positive_roots(A, p_min), positive_roots(B, p_min)
    common = [a for a in rA if any(abs(a - b) <= tol * max(1.0, a) for b in rB)]
    accepted = [p for p in common if consistency_residual(p, v, dl) <= tol * max(1.0, p)]
    return dict(A=A, B=B, roots_A=rA, roots_B=rB, common=common, accepted=accepted)


# ----------------------------------------------------------------------------
# physical validity of a reflection path (convex room bounded by 3 lines)
# ----------------------------------------------------------------------------
def orient(walls):
    """normals point out of the room, i.e. the room is n . x <= d (R inside)"""
    out = []
    for n, d in walls:
        if n @ R - d > 0:
            n, d = -n, -d
        out.append((n, d))
    return out


def inside(x, walls, tol=1e-7):
    return all(n @ x - d <= tol for n, d in walls)


def path_points(S, walls, lab):
    imgs = [np.array(S, float)]
    for k in lab:
        imgs.append(reflect(imgs[-1], walls[k - 1]))
    pts, target = [R.copy()], R.copy()
    for idx in range(len(lab) - 1, -1, -1):
        n, d = walls[lab[idx] - 1]
        u = imgs[idx + 1] - target
        den = n @ u
        if abs(den) < 1e-12:
            return None
        t = (d - n @ target) / den
        if not 0.0 < t < 1.0:
            return None
        target = target + t * u
        pts.append(target)
    pts.append(np.array(S, float))
    return pts[::-1]


def is_physical(S, walls, lab):
    w = orient(walls)
    pts = path_points(S, w, lab)
    return pts is not None and inside(S, w) and all(inside(q, w) for q in pts[1:-1])


# ----------------------------------------------------------------------------
# measurement model with parameters, used for the Jacobian and the noise test
# theta = [Sx, Sy, phi1, d1, phi2, d2, phi3, d3],  n_k = (cos phi_k, sin phi_k)
# ----------------------------------------------------------------------------
def theta_from(S, walls):
    th = [S[0], S[1]]
    for n, d in walls:
        th += [np.arctan2(n[1], n[0]), d]
    return np.array(th, float)


def walls_from(th):
    return [(np.array([np.cos(th[2 + 2 * k]), np.sin(th[2 + 2 * k])]), th[3 + 2 * k]) for k in range(3)]


def forward(th, labels=LABELS):
    v, dl = measure(th[:2], walls_from(th), labels)
    return np.concatenate([np.arctan2(v[:, 1], v[:, 0]), dl[1:]])


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def num_jacobian(f, th, eps=1e-6):
    cols = []
    for e in np.eye(len(th)):
        df = f(th + eps * e) - f(th - eps * e)
        cols.append(df / (2 * eps))
    return np.array(cols).T


def rank_info(J):
    s = np.linalg.svd(J, compute_uv=False)
    return int(np.sum(s > 1e-8 * s[0])), s


def refine(th0, meas, iters=60):
    """small Levenberg-Marquardt fit of all 8 parameters to all 11 measurements"""
    th, lam = th0.copy(), 1e-3

    def res(t):
        r = forward(t) - meas
        r[:6] = wrap(r[:6])
        return r

    r = res(th)
    for _ in range(iters):
        J = num_jacobian(res, th)
        H = J.T @ J
        step = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-12), -J.T @ r)
        r_new = res(th + step)
        if r_new @ r_new < r @ r:
            th, r, lam = th + step, r_new, lam * 0.3
        else:
            lam *= 10.0
        if np.linalg.norm(step) < 1e-12:
            break
    return th


def solve_noisy(v, dl):
    """closed form under noise: p minimises A(p)^2 + B(p)^2"""
    p_min = max(0.0, -float(dl.min())) + 1e-6
    A = quad_coeffs(v, dl, I0, I12, I0, I21)
    B = quad_coeffs(v, dl, I13, I23, I0, I12)
    F = lambda p: np.polyval(A, p) ** 2 + np.polyval(B, p) ** 2
    grid = np.linspace(p_min, p_min + 200.0, 20001)
    k = int(np.argmin(F(grid)))
    lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
    for _ in range(80):                                   # golden section
        m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if F(m1) < F(m2):
            hi = m2
        else:
            lo = m1
    return 0.5 * (lo + hi)


def random_room(rng):
    while True:
        S = rng.uniform(-10, 10, 2)
        walls = []
        for _ in range(3):
            a = rng.uniform(0, 2 * np.pi)
            walls.append((np.array([np.cos(a), np.sin(a)]), rng.uniform(1, 15)))
        if np.linalg.norm(S) > 0.5 and all(abs(n @ S - d) > 0.5 for n, d in walls):
            return S, walls


# ============================================================================
if __name__ == "__main__":
    results = {}
    rng = np.random.default_rng(2026)

    # ---------------------------------------------------------------- example
    S_true = np.array([0.0, 13.0])
    wall_pts = [((-7, 15), (-11, -8)), ((3, 15), (9, 7)), ((9, 7), (3, -8))]
    walls_true = [line_from_points(a, b) for a, b in wall_pts]

    print("=" * 70)
    print("Example room")
    print("=" * 70)
    all_labels = [()] + [(a,) for a in (1, 2, 3)] + [(a, b) for a, b in product((1, 2, 3), repeat=2) if a != b]
    phys = {str(lab): is_physical(S_true, walls_true, lab) for lab in all_labels}
    print("physically possible paths up to order 2:", [k for k, ok in phys.items() if ok])
    assert all(phys[str(lab)] for lab in LABELS), "example room does not produce all six signals"

    v, dl = measure(S_true, walls_true)
    sol = solve(v, dl)
    p_hat = sol["accepted"][0]
    S_hat, walls_hat, P = reconstruct(p_hat, v, dl)
    errs = [wall_error(we, wt) for we, wt in zip(walls_hat, walls_true)]
    meas_table = []
    for name, vi, di, Pi in zip(NAMES, v, dl, image_points(S_true, walls_true)):
        ang = float(np.degrees(np.arctan2(vi[1], vi[0])))
        meas_table.append(dict(signal=name, angle_deg=ang, delta=float(di), image=Pi.tolist()))
        print(f"  {name:8s}  angle = {ang:9.4f} deg   delta = {di:8.4f}   image = {Pi.round(4)}")
    print("  (A) coefficients:", sol["A"].round(6), " positive roots:", np.round(sol["roots_A"], 6))
    print("  (B) coefficients:", sol["B"].round(6), " positive roots:", np.round(sol["roots_B"], 6))
    print("  accepted p:", p_hat, "  true p:", np.linalg.norm(S_true))
    print("  source error:", np.linalg.norm(S_hat - S_true))
    print("  wall errors (deg, distance):", [(f"{a:.1e}", f"{b:.1e}") for a, b in errs])
    print("  consistency residual at p:", consistency_residual(p_hat, v, dl))
    results["example"] = dict(
        S=S_true.tolist(), wall_points=wall_pts, physical_paths=phys,
        measurements=meas_table, A=sol["A"].tolist(), B=sol["B"].tolist(),
        roots_A=sol["roots_A"], roots_B=sol["roots_B"], p_hat=p_hat,
        source_error=float(np.linalg.norm(S_hat - S_true)),
        wall_errors=[list(map(float, e)) for e in errs],
        consistency_residual=consistency_residual(p_hat, v, dl))

    # -------------------------------------------------------- local rank test
    th_true = theta_from(S_true, walls_true)
    J = num_jacobian(forward, th_true)
    r_all, s_all = rank_info(J)
    r_dir, _ = rank_info(J[:6])
    r_tdoa, _ = rank_info(J[6:])
    print("\nJacobian 11x8: rank", r_all, " singular values", s_all.round(4))
    print("  directions only (6x8): rank", r_dir, "   TDOAs only (5x8): rank", r_tdoa)
    drop = {}
    for k in range(1, 6):                      # the direct signal is needed as reference
        labs = [lab for i, lab in enumerate(LABELS) if i != k]
        Jk = num_jacobian(lambda t: forward(t, labs), th_true)
        drop[NAMES[k]] = dict(measurements=Jk.shape[0], rank=rank_info(Jk)[0])
    print("  without one signal:", drop)
    results["jacobian"] = dict(rank=r_all, singular_values=s_all.tolist(),
                               rank_directions_only=r_dir, rank_tdoa_only=r_tdoa, drop_one=drop)

    # ------------------------------------------------------- random rooms
    N = 5000
    stats = dict(rooms=N, exactly_one_accepted=0, accepted_is_true=0,
                 A_two_roots=0, B_two_roots=0, either_two_roots=0, common_more_than_one=0,
                 max_source_error=0.0, max_wall_angle_error=0.0)
    spurious_example = None
    for _ in range(N):
        S, walls = random_room(rng)
        v_r, dl_r = measure(S, walls)
        s_r = solve(v_r, dl_r)
        p_true = float(np.linalg.norm(S))
        stats["A_two_roots"] += len(s_r["roots_A"]) > 1
        stats["B_two_roots"] += len(s_r["roots_B"]) > 1
        stats["either_two_roots"] += (len(s_r["roots_A"]) > 1) or (len(s_r["roots_B"]) > 1)
        stats["common_more_than_one"] += len(s_r["common"]) > 1
        if len(s_r["accepted"]) == 1:
            stats["exactly_one_accepted"] += 1
            p_r = s_r["accepted"][0]
            if abs(p_r - p_true) < 1e-6 * max(1, p_true):
                stats["accepted_is_true"] += 1
            S_r, W_r, _ = reconstruct(p_r, v_r, dl_r)
            stats["max_source_error"] = max(stats["max_source_error"], float(np.linalg.norm(S_r - S)))
            stats["max_wall_angle_error"] = max(stats["max_wall_angle_error"],
                                                max(wall_error(a, b)[0] for a, b in zip(W_r, walls)))
        if spurious_example is None and 5 < p_true < 15:
            wrong = [r for r in s_r["roots_A"] + s_r["roots_B"] if abs(r - p_true) > 1e-6 * p_true]
            if wrong and all(1.3 * p_true < r < 4.0 * p_true for r in wrong):
                spurious_example = (S, walls, v_r, dl_r, s_r)
    print("\nRandom rooms:", stats)
    results["random_rooms"] = stats

    # ------------------------------------------------------- noise study
    levels = [(0.01, 0.001), (0.1, 0.01), (0.5, 0.05), (1.0, 0.1)]
    trials = 300
    meas_true = forward(th_true)
    noise_rows = []
    for s_ang, s_del in levels:
        e_cf, e_ref, w_cf, w_ref = [], [], [], []
        for _ in range(trials):
            ang = meas_true[:6] + np.radians(rng.normal(0, s_ang, 6))
            dln = np.concatenate([[0.0], meas_true[6:] + rng.normal(0, s_del, 5)])
            vn = np.stack([np.cos(ang), np.sin(ang)], axis=1)
            p_n = solve_noisy(vn, dln)
            S_n, W_n, _ = reconstruct(p_n, vn, dln)
            e_cf.append(np.linalg.norm(S_n - S_true))
            w_cf.append(max(wall_error(a, b)[0] for a, b in zip(W_n, walls_true)))
            th_r = refine(theta_from(S_n, W_n), np.concatenate([ang, dln[1:]]))
            e_ref.append(np.linalg.norm(th_r[:2] - S_true))
            w_ref.append(max(wall_error(a, b)[0] for a, b in zip(walls_from(th_r), walls_true)))
        row = dict(sigma_angle_deg=s_ang, sigma_delta=s_del,
                   src_median_closed=float(np.median(e_cf)), src_p90_closed=float(np.percentile(e_cf, 90)),
                   wall_median_closed=float(np.median(w_cf)),
                   src_median_refined=float(np.median(e_ref)), src_p90_refined=float(np.percentile(e_ref, 90)),
                   wall_median_refined=float(np.median(w_ref)))
        noise_rows.append(row)
        print("noise", row)
    results["noise"] = dict(trials=trials, rows=noise_rows)

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ===================================================================
    # figures
    # ===================================================================
    C_WALL, C_S, C_R, C_IMG = "#1f5fa8", "#c0392b", "#2e8b57", "#e08a00"
    imgs = image_points(S_true, walls_true)

    def draw_walls(ax, walls, xlim, ylim, labels=True):
        for k, (n, d) in enumerate(walls):
            t = np.array([-n[1], n[0]])
            p0 = d * n
            pts = np.array([p0 + s * t for s in (-200, 200)])
            ax.plot(pts[:, 0], pts[:, 1], color=C_WALL, lw=1.8, zorder=2)
        if labels:
            for k, (a, b) in enumerate(wall_pts):
                m = 0.5 * (np.array(a) + np.array(b))
                off = {0: (-1.6, 0.0), 1: (1.2, 0.8), 2: (1.2, -0.6)}[k]
                ax.text(m[0] + off[0], m[1] + off[1], f"$w_{k+1}$", color=C_WALL, fontsize=11,
                        ha="center", va="center")
        ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")


    LABEL_OFF = {I0: (1.0, 1.0), I1: (-3.5, -3.6), I12: (1.0, 1.0), I21: (-3.0, 2.0),
                 I13: (-4.0, -3.8), I23: (1.2, 0.8)}

    def label_walls(ax, walls, ys=(36.0, -8.0, 36.0), dx=(1.2, 1.2, 1.2)):
        for k, (n, d) in enumerate(walls):
            t = np.array([-n[1], n[0]])
            p0 = d * n
            q = p0 + (ys[k] - p0[1]) / t[1] * t
            ax.text(q[0] + dx[k], q[1], f"$w_{k+1}$", color=C_WALL, fontsize=11)

    def room_polygon(walls, box=60):
        """clip a big square by the three half-planes of the room"""
        poly = [np.array(p, float) for p in [(-box, -box), (box, -box), (box, box), (-box, box)]]
        for n, d in orient(walls):
            new = []
            for i in range(len(poly)):
                a, b = poly[i], poly[(i + 1) % len(poly)]
                fa, fb = n @ a - d, n @ b - d
                if fa <= 0: new.append(a)
                if fa * fb < 0:
                    new.append(a + (b - a) * fa / (fa - fb))
            poly = new
        return np.array(poly)

    # ---- Figure: example room with the six paths, and the image points
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.9))
    ax = axes[0]
    poly = room_polygon(walls_true)
    ax.fill(poly[:, 0], poly[:, 1], color="#eef2f7", zorder=0)
    draw_walls(ax, walls_true, (-14, 13), (-10, 18))
    styles = ["-", "-", "--", "--", ":", ":"]
    for i, lab in enumerate(LABELS):
        pts = np.array(path_points(S_true, orient(walls_true), lab))
        ax.plot(pts[:, 0], pts[:, 1], color="0.35", lw=0.9, ls=styles[i], zorder=3)
    ax.scatter(*S_true, s=60, color=C_S, zorder=5); ax.text(0.6, 13.6, "$S$", fontsize=11)
    ax.scatter(*R, s=60, color=C_R, zorder=5); ax.text(0.6, -1.3, "$R$", fontsize=11)
    ax.set_title("(a) room and the six signal paths")
    ax.set_xticks([-10, 0, 10]); ax.set_yticks([-10, 0, 10])

    ax = axes[1]
    draw_walls(ax, walls_true, (-22, 34), (-14, 42), labels=False)
    for (a, b, k) in [(I0, I1, 1), (I1, I12, 2), (I1, I13, 3)]:
        pa, pb = imgs[a], imgs[b]
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="0.25", ls="--", lw=0.9, zorder=3)
        m = 0.5 * (pa + pb)
        ax.scatter(*m, s=14, color=C_WALL, zorder=4)
    for i, Pi in enumerate(imgs):
        ax.plot([0, Pi[0]], [0, Pi[1]], color=C_IMG, lw=0.6, alpha=0.6, zorder=1)
        ax.scatter(*Pi, s=34, color=C_S if i == 0 else C_IMG, zorder=5, edgecolor="k", linewidth=0.4)
        off = LABEL_OFF[i]
        ax.text(Pi[0] + off[0], Pi[1] + off[1], f"${TEX[i]}$", fontsize=10)
    ax.scatter(*R, s=40, color=C_R, zorder=5)
    ax.text(1.2, -3.2, "$R$", fontsize=10)
    label_walls(ax, walls_true)
    ax.set_title("(b) image points seen from $R$")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_example_room.pdf"))
    fig.savefig(os.path.join(OUT, "fig_example_room.png"), dpi=200)
    plt.close(fig)

    # ---- Figure: equal distances used in (A) and (B)
    I2 = reflect(S_true, walls_true[1])
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    draw_walls(ax, walls_true, (-22, 34), (-14, 42), labels=False)
    pairs = [(imgs[I0], imgs[I12], r"$|I_0 I_{12}|$"), (imgs[I0], imgs[I21], r"$|I_0 I_{21}|$"),
             (imgs[I13], imgs[I23], r"$|I_{13} I_{23}|$"), (imgs[I1], I2, r"$|I_1 I_2|$ (not measured)")]
    cols = ["#8e44ad", "#8e44ad", "#16a085", "0.5"]
    lss = ["-", "--", "-", ":"]
    for (a, b, lab), c, ls in zip(pairs, cols, lss):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=c, ls=ls, lw=1.6, label=lab, zorder=3)
    for i, Pi in enumerate(imgs):
        ax.scatter(*Pi, s=30, color=C_S if i == 0 else C_IMG, zorder=5, edgecolor="k", linewidth=0.4)
        off = LABEL_OFF[i]
        ax.text(Pi[0] + off[0], Pi[1] + off[1], f"${TEX[i]}$", fontsize=10)
    ax.scatter(*I2, s=30, facecolor="white", edgecolor="0.3", zorder=5)
    ax.text(I2[0] + 1.0, I2[1] - 3.2, "$I_2$", fontsize=10, color="0.3")
    ax.scatter(*R, s=40, color=C_R, zorder=5)
    ax.text(1.2, -3.2, "$R$", fontsize=10)
    label_walls(ax, walls_true)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_equal_distances.pdf"))
    fig.savefig(os.path.join(OUT, "fig_equal_distances.png"), dpi=200)
    plt.close(fig)

    # ---- Figure: the two quadratics, example room and a room with a spurious root
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.7))
    cases = [("(a) example room", v, dl, sol, np.linalg.norm(S_true))]
    if spurious_example is not None:
        S_sp, W_sp, v_sp, dl_sp, s_sp = spurious_example
        cases.append(("(b) random room", v_sp, dl_sp, s_sp, np.linalg.norm(S_sp)))
        results["spurious_example"] = dict(S=S_sp.tolist(), roots_A=s_sp["roots_A"], roots_B=s_sp["roots_B"],
                                           p_true=float(np.linalg.norm(S_sp)))
    for ax, (title, vv, dd, ss, pt) in zip(axes, cases):
        allr = ss["roots_A"] + ss["roots_B"] + [pt]
        pmax = 1.35 * max(allr)
        p = np.linspace(0, pmax, 600)
        fa, fb = np.polyval(ss["A"], p), np.polyval(ss["B"], p)
        ax.axhline(0, color="0.6", lw=0.8)
        ax.plot(p, fa / np.abs(fa).max(), color="#8e44ad", label="(A)")
        ax.plot(p, fb / np.abs(fb).max(), color="#16a085", ls="--", label="(B)")
        for r in ss["roots_A"]:
            ax.scatter(r, 0, s=26, color="#8e44ad", zorder=4)
        for r in ss["roots_B"]:
            ax.scatter(r, 0, s=60, facecolor="none", edgecolor="#16a085", zorder=4)
        ax.axvline(pt, color=C_S, lw=0.8, ls=":", label="true $p$")
        ax.set_title(title); ax.set_xlabel("$p$")
        ax.set_yticks([0]); ax.set_ylabel("scaled value")
        ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_quadratics.pdf"))
    fig.savefig(os.path.join(OUT, "fig_quadratics.png"), dpi=200)
    plt.close(fig)

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("\nfigures and results written to", OUT)
