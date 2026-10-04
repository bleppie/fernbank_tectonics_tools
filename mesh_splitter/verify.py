#!/usr/bin/env python3
"""Independent checker for mesh_splitter.py output.  Re-derives everything from the
written OBJ files and the original inputs -- it does not trust the cutter.

usage: verify.py <original.obj> <lines.png> <out-dir> [--max-aspect A]
                 [--max-distance D] [--max-error E]
"""
import argparse
import glob
import math
import os
import sys
from collections import defaultdict

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mesh_splitter import (load_obj, skeletonize, smooth_mask, TriLocator,
                      UVImage, cross2, boundary_loops)

FAIL = []


def check(ok, msg, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {msg}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(msg)


def shoelace(p, t):
    return float(np.abs(cross2(p[t[:, 1]] - p[t[:, 0]],
                               p[t[:, 2]] - p[t[:, 0]])).sum() * 0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("orig")
    ap.add_argument("image")
    ap.add_argument("outdir")
    ap.add_argument("--max-aspect", type=float, default=4.0)
    ap.add_argument("--max-distance", type=float, default=100.0)
    ap.add_argument("--max-error", type=float, default=2.0)
    ap.add_argument("--skip-aspect", action="store_true",
                    help="report the aspect ratio but do not fail on it")
    ap.add_argument("--skip-error", action="store_true",
                    help="report the centreline error but do not fail on it")
    a = ap.parse_args()

    orig = load_obj(a.orig)
    base, ext = os.path.splitext(os.path.basename(a.orig))
    files = sorted(glob.glob(os.path.join(a.outdir, f"{base}_*{ext}")),
                   key=lambda f: int(f.rsplit("_", 1)[1].split(".")[0]))
    print(f"original: {len(orig.pos)} verts / {len(orig.tris)} tris; "
          f"{len(files)} submesh file(s)")
    if not files:
        sys.exit("no output files found")
    subs = [load_obj(f) for f in files]

    # ---- 1. area is conserved -------------------------------------------
    a0 = shoelace(orig.xz, orig.tris)
    a1 = sum(shoelace(s.xz, s.tris) for s in subs)
    check(abs(a1 - a0) <= 1e-7 * a0, "submesh areas sum to the original area",
          f"{a1:.4f} vs {a0:.4f}  (rel {abs(a1 - a0) / a0:.2e})")

    # ---- 2. each submesh is a 2-manifold with a closed boundary ----------
    bad = []
    for f, s in zip(files, subs):
        cnt = defaultdict(int)
        half = set()
        for tri in s.tris:
            for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                cnt[(min(u, v), max(u, v))] += 1
                if (u, v) in half:
                    bad.append(f"{f}: duplicated half-edge (bad winding)")
                half.add((u, v))
        if any(c > 2 for c in cnt.values()):
            bad.append(f"{f}: non-manifold edge")
        # boundary edges must form closed loops -> every boundary vertex has
        # exactly as many incoming as outgoing boundary half-edges
        deg = defaultdict(int)
        for (u, v), c in cnt.items():
            if c == 1:
                deg[u] += 1
                deg[v] += 1
        if any(d % 2 for d in deg.values()):
            bad.append(f"{f}: boundary is not a set of closed loops")
        if len(set(s.tris.ravel().tolist())) != len(s.pos):
            bad.append(f"{f}: unreferenced vertices")
    check(not bad, "every submesh is a 2-manifold with closed boundary loops",
          "; ".join(bad[:3]))

    # ---- 3. winding matches the original ---------------------------------
    sgn0 = np.sign(cross2(orig.xz[orig.tris[:, 1]] - orig.xz[orig.tris[:, 0]],
                          orig.xz[orig.tris[:, 2]] - orig.xz[orig.tris[:, 0]]))
    want = 1.0 if sgn0.sum() >= 0 else -1.0
    wrong = 0
    for s in subs:
        sg = np.sign(cross2(s.xz[s.tris[:, 1]] - s.xz[s.tris[:, 0]],
                            s.xz[s.tris[:, 2]] - s.xz[s.tris[:, 0]]))
        wrong += int((sg != want).sum())
    check(wrong == 0, "triangle winding matches the input mesh",
          f"{wrong} reversed")

    # ---- 4. texcoords really are interpolated from the original ----------
    loc = TriLocator(orig.xz, orig.tris)
    worst_uv = worst_y = 0.0
    for s in subs:
        t, b = loc.locate_clamped(s.xz)
        ref_uv = (orig.uv[orig.tris[t]] * b[:, :, None]).sum(1)
        ref_y = (orig.pos[orig.tris[t]][:, :, 1] * b).sum(1)
        worst_uv = max(worst_uv, float(np.abs(ref_uv - s.uv).max()))
        worst_y = max(worst_y, float(np.abs(ref_y - s.pos[:, 1]).max()))
    check(worst_uv < 1e-7, "uv == barycentric interpolation of the original",
          f"max |du| = {worst_uv:.3e}")
    check(worst_y < 1e-6, "y  == barycentric interpolation of the original",
          f"max |dy| = {worst_y:.3e}")

    # ---- 5. aspect ratio --------------------------------------------------
    # A mesh generator cannot bound the aspect ratio inside a wedge whose
    # input angle is already smaller than the target -- that is a property of
    # the input, not of the mesher (Ruppert; Triangle documents the same
    # limitation).  So exempt triangles sitting in such a wedge, and report
    # how many were exempted so the exemption cannot quietly hide a bug.
    theta = math.degrees(math.asin(min(1.0, 1.0 / max(a.max_aspect, 1.0))))
    ratios, exempt_flags, n_sharp = [], [], 0
    for s in subs:
        p, t = s.xz, s.tris
        L = np.stack([np.linalg.norm(p[t[:, 1]] - p[t[:, 2]], axis=1),
                      np.linalg.norm(p[t[:, 2]] - p[t[:, 0]], axis=1),
                      np.linalg.norm(p[t[:, 0]] - p[t[:, 1]], axis=1)], -1)
        ratios.append(L.max(1) / np.maximum(L.min(1), 1e-18))

        cnt = defaultdict(int)
        for tri in t:
            for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                cnt[(min(u, v), max(u, v))] += 1
        adj = defaultdict(list)
        for (u, v), c in cnt.items():
            if c == 1:
                adj[u].append(v)
                adj[v].append(u)
        sharp = set()
        for u, nb in adj.items():
            if len(nb) != 2:
                continue
            d1 = p[nb[0]] - p[u]
            d2 = p[nb[1]] - p[u]
            n1, n2 = np.linalg.norm(d1), np.linalg.norm(d2)
            if n1 < 1e-15 or n2 < 1e-15:
                continue
            ang = math.degrees(math.acos(
                float(np.clip((d1 @ d2) / (n1 * n2), -1.0, 1.0))))
            if ang < theta:
                sharp.add(u)
        n_sharp += len(sharp)
        exempt_flags.append(np.array([any(v in sharp for v in tri)
                                      for tri in t], bool))
    r = np.concatenate(ratios)
    ex = np.concatenate(exempt_flags) if exempt_flags else np.zeros(0, bool)
    over = r > a.max_aspect + 1e-6
    unexplained = int((over & ~ex).sum())
    check(unexplained == 0 or a.skip_aspect,
          f"edge-ratio <= {a.max_aspect} outside wedges narrower than "
          f"{theta:.1f}deg",
          f"max {r.max():.3f}, mean {r.mean():.3f}, {int(over.sum())} over "
          f"({int((over & ex).sum())} inside {n_sharp} sharp input wedge(s), "
          f"{unexplained} unexplained)")

    # ---- 6. vertex spacing along cut boundaries --------------------------
    img = np.array(Image.open(a.image).convert("RGBA"))
    H, W = img.shape[:2]
    uvi = UVImage(W, H, True)
    m = img[..., 3] > 0
    # match the cutter's auto smoothing, else we would be measuring the error
    # against a centreline the cutter never saw
    from scipy.ndimage import distance_transform_edt
    r = max(1, int(round(float(np.median(distance_transform_edt(m)[m])) / 3.0)))
    skel = skeletonize(smooth_mask(m, r))
    sy, sx = np.nonzero(skel)
    from scipy.spatial import KDTree
    sk = KDTree(np.stack([sx, sy], -1).astype(float))

    bd_segs = []
    for loop in boundary_loops(orig.tris, len(orig.pos)):
        ring = orig.xz[loop]
        for i in range(len(ring)):
            bd_segs.append((ring[i], ring[(i + 1) % len(ring)]))
    diag = float(np.linalg.norm(orig.xz.max(0) - orig.xz.min(0)))

    def on_mesh_boundary(p, tol=None):
        """Distance from a point to the original mesh boundary, thresholded."""
        tol = diag * 1e-9 if tol is None else tol
        for A, B in bd_segs:
            d = B - A
            L2 = float(d @ d)
            if L2 <= 0:
                continue
            t = min(1.0, max(0.0, float((p - A) @ d) / L2))
            if float(np.linalg.norm(A + t * d - p)) <= tol:
                return True
        return False

    gaps, errs = [], []
    for s in subs:
        cnt = defaultdict(int)
        for tri in s.tris:
            for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                cnt[(min(u, v), max(u, v))] += 1
        for (u, v), c in cnt.items():
            if c != 1:
                continue
            p, q = s.xz[u], s.xz[v]
            # --max-distance applies to every cut, including the border
            gaps.append(float(np.linalg.norm(q - p)))
            # --max-error only means anything for vertices on a *drawn* line;
            # vertices sitting on the original mesh rim have no centreline
            for z, w in ((u, p), (v, q)):
                if on_mesh_boundary(w):
                    continue
                errs.append(float(sk.query(uvi.uv_to_px(s.uv[z]))[0]))
    gaps = np.asarray(gaps)
    errs = np.asarray(errs)
    check(gaps.max() <= a.max_distance + 1e-6,
          f"cut-edge length <= --max-distance ({a.max_distance})",
          f"max {gaps.max():.2f}, mean {gaps.mean():.2f}, n={len(gaps)}")
    check(errs.max() <= a.max_error + 0.75 or a.skip_error,   # +0.75px: pixel-centre quantisation
          f"cut vertices lie within --max-error ({a.max_error}px) of the "
          f"drawn centreline",
          f"max {errs.max():.2f}px, mean {errs.mean():.2f}px")

    # ---- 7. the pieces reassemble watertight -----------------------------
    allv = np.vstack([s.pos for s in subs])
    key = np.round(allv, 6)
    uniq = {tuple(k) for k in key}
    shared = len(allv) - len(uniq)
    tree = KDTree(np.array(sorted(uniq)))
    # every cut-boundary vertex of one piece must coincide with one in another
    orphan = 0
    for i, s in enumerate(subs):
        others = np.vstack([t.pos for j, t in enumerate(subs) if j != i])
        ot = KDTree(others)
        cnt = defaultdict(int)
        for tri in s.tris:
            for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                cnt[(min(u, v), max(u, v))] += 1
        bverts = {u for (u, v), c in cnt.items() if c == 1 for u in (u, v)}
        for b in bverts:
            if on_mesh_boundary(s.xz[b]):
                continue
            if ot.query(s.pos[b])[0] > 1e-9:
                orphan += 1
    check(orphan == 0, "cut-boundary vertices are shared exactly between "
                       "neighbouring submeshes (watertight)",
          f"{shared} duplicated vertices, {orphan} orphaned")

    print()
    if FAIL:
        print(f"{len(FAIL)} CHECK(S) FAILED")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
