#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.22",
#     "scipy>=1.10",
# ]
# ///
"""Synthetic end-to-end test for mesh_munger.py.

Builds a faceted, triangle-soup "island" (a Gaussian hill on a square grid,
rotated 25 degrees, edges dipping below the no-back threshold), runs the
full pipeline, and checks every output:

  * the bounding rectangle recovers the grid's size and rotation
  * welding collapses the soup to one vertex per grid point
  * _no_back keeps only faces at/above the threshold
  * _extended is watertight against _no_back (one open loop: the outer rim),
    has no folded triangles, slopes down outward, and 'ease' ends higher
    than 'linear'
  * skirt UVs continue the model's mapping (outside 0..1)

Run:  python selftest.py   (or: uv run selftest.py)
"""

import contextlib
import io
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mesh_munger as mm  # noqa: E402

SIZE, N, ANGLE = 1000.0, 41, 25.0


def make_island(path):
    g = np.linspace(-SIZE / 2, SIZE / 2, N)
    X, Z = np.meshgrid(g, g)
    Y = 20.0 + 250.0 * np.exp(-(X ** 2 + Z ** 2) / (2 * 250.0 ** 2))
    a = math.radians(ANGLE)
    Xr, Zr = X * math.cos(a) - Z * math.sin(a) + 300, X * math.sin(a) + Z * math.cos(a) - 200
    P = np.stack([Xr, Y, Zr], -1).reshape(-1, 3)
    idx = np.arange(N * N).reshape(N, N)
    q = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel(), idx[1:, :-1].ravel()
    tris = np.concatenate([np.stack([q[0], q[3], q[2]], 1), np.stack([q[0], q[2], q[1]], 1)])
    with open(path, "w") as f:                    # triangle soup, flat normals, colours
        f.write("o island\n")
        k = 0
        for t in tris:
            a_, b_, c_ = P[t]
            n = np.cross(b_ - a_, c_ - a_)
            n /= np.linalg.norm(n)
            assert n[1] > 0
            for p in (a_, b_, c_):
                f.write(f"v {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} 0.5 0.5 0.5\n")
            for _ in range(3):
                f.write(f"vn {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")
            f.write(f"f {k+1}//{k+1} {k+2}//{k+2} {k+3}//{k+3}\n")
            k += 3
    return len(tris)


def run(args):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        mm.main(args)
    return out.getvalue()


def check(cond, msg):
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond:
        check.failed = True


check.failed = False


def main():
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "island.obj")
        n_tris = make_island(src)
        log = run([src, "--extend-distance", "80"])
        base = os.path.join(d, "island")

        print("measurements")
        mesh = mm.read_obj(src)
        rect = mm.min_bounding_rectangle(mm.convex_hull(mm.weld(mesh).v[:, [0, 2]]))
        check(abs(rect.width - SIZE) < 1e-3 and abs(rect.height - SIZE) < 1e-3,
              f"rectangle {rect.width:.3f} x {rect.height:.3f} (expect {SIZE:g})")
        check(abs(rect.angle - ANGLE) < 1e-6 or abs(rect.angle - (90 - ANGLE)) < 1e-6,
              f"rotation {rect.angle:.3f} deg (expect {ANGLE:g} or {90 - ANGLE:g})")

        print("welding")
        welded = mm.weld(mesh)
        check(len(mesh.v) == 3 * n_tris, f"read soup: {len(mesh.v)} vertices")
        check(len(welded.v) == N * N, f"welded to {len(welded.v)} (expect {N * N})")
        check(len(mm.weld(mesh, 0).v) > N * N, "--smooth-angle 0 keeps facets split")

        print("outputs")
        rot = mm.read_obj(base + "_rotated_uvs.obj")
        check(len(rot.faces) == n_tris, "_rotated_uvs keeps every face")
        check(np.allclose(rot.v[:, [0, 2]].min(0), -SIZE / 2, atol=1e-3)
              and np.allclose(rot.v[:, [0, 2]].max(0), SIZE / 2, atol=1e-3),
              "_rotated_uvs is axis-aligned and centered")
        nb = mm.read_obj(base + "_no_back.obj")
        check(nb.v[:, 1].min() >= 50.0 and len(nb.faces) < n_tris,
              f"_no_back keeps {len(nb.faces)} faces, all at Y >= 50")
        hull = mm.read_obj(base + "_hull.obj")
        check(len(hull.v) == 4 and (hull.face_normals()[:, 1] > 0).all(),
              "_hull is a 4-corner, upward-facing quad")

        print("extension")
        ext = mm.read_obj(base + "_extended.obj")
        P, loop, _e, _n, n_loops = mm.outer_border(ext)
        check(n_loops == 1, f"_extended has a single open border ({n_loops})")
        skirt = mm.Mesh(ext.v, ext.faces[len(nb.faces):])
        check((skirt.face_normals()[:, 1] > 0).all(), "no folded skirt triangles")
        rim_r = np.linalg.norm(P[loop][:, [0, 2]], axis=1)
        nbP, nbl = mm.outer_border(nb)[:2]
        nb_r = np.linalg.norm(nbP[nbl][:, [0, 2]], axis=1)
        check(rim_r.min() > nb_r.max() - 5, f"rim is outside the border "
              f"(rim r >= {rim_r.min():.1f}, border r <= {nb_r.max():.1f})")
        check(P[loop][:, 1].max() < nbP[nbl][:, 1].min(), "skirt slopes down outward")
        rim_y_linear = P[loop][:, 1].mean()
        txt = open(base + "_extended.obj").read()
        vts = np.array([l.split()[1:3] for l in txt.splitlines() if l.startswith("vt ")], float)
        check(vts.min() < 0 or vts.max() > 1, "skirt UVs extend past 0..1")
        check("WARNING" not in log, "no fold warning printed")

        run([src, "--extend-distance", "80", "--extend-profile", "ease"])
        e2 = mm.read_obj(base + "_extended.obj")
        P2, l2 = mm.outer_border(e2)[:2]
        check(P2[l2][:, 1].mean() > rim_y_linear,
              f"ease rim Y {P2[l2][:, 1].mean():.1f} > linear rim Y {rim_y_linear:.1f}")

    print("\nFAILED" if check.failed else "\nall checks passed")
    sys.exit(1 if check.failed else 0)


if __name__ == "__main__":
    main()
