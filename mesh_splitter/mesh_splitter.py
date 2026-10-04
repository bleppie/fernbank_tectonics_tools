#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.22",
#     "scipy>=1.10",
#     "pillow>=9.0",
# ]
# ///
"""
mesh_splitter.py -- cut a 3D OBJ mesh into submeshes along lines drawn on a PNG.

The mesh is treated as a 2D manifold in the XZ plane (Y is carried along and
interpolated, but never used for geometry).  The mesh's texture coordinates are
a bijection onto the image, so every pixel maps to a point on the mesh and vice
versa, piecewise-affinely, one affine map per triangle.

Pipeline
  1. parse OBJ, parse PNG (line = alpha > threshold)
  2. thin the strokes to 1px centrelines (Zhang-Suen), trace them into polylines
  3. prune spurs, simplify with Douglas-Peucker (--max-error, pixels)
  4. snap dangling endpoints onto nearby lines or the mesh rim (--snap)
  5. clip every line to the mesh's UV footprint; add the mesh boundary itself
     as a constraint
  6. split each segment where it crosses a mesh triangle edge, lift to XZ
  7. densify so no two consecutive vertices exceed --max-distance
  8. constrained Delaunay (scipy + encroachment splitting), then Ruppert
     refinement to honour --max-aspect
  9. flood-fill triangles, walls = constraint edges -> N regions
 10. interpolate uv / y / normals from the original mesh, write one OBJ each

Dependencies: numpy, scipy, pillow.  All permissively licensed.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np
from PIL import Image
from scipy.spatial import Delaunay, KDTree

# --------------------------------------------------------------------------
# small utilities
# --------------------------------------------------------------------------

T0 = time.time()
VERBOSE = True


def log(msg):
    if VERBOSE:
        print(f"[{time.time() - T0:6.2f}s] {msg}", flush=True)


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def points_in_loops(pts, loops):
    """Even-odd ray casting.  pts (N,2), loops = list of (K,2) closed rings."""
    inside = np.zeros(len(pts), bool)
    x, y = pts[:, 0], pts[:, 1]
    for ring in loops:
        a = ring
        b = np.roll(ring, -1, axis=0)
        for (ax, ay), (bx, by) in zip(a, b):
            if ay == by:
                continue
            cond = ((ay > y) != (by > y))
            if not cond.any():
                continue
            xint = ax + (y - ay) * (bx - ax) / (by - ay)
            inside ^= cond & (x < xint)
    return inside


# --------------------------------------------------------------------------
# 1. OBJ I/O
# --------------------------------------------------------------------------

class Mesh:
    """Unified mesh: one vertex per distinct (v, vt, vn) triple."""

    def __init__(self, pos, uv, nrm, tris):
        self.pos = pos            # (N,3) float64  x,y,z
        self.uv = uv              # (N,2) float64  u,v in OBJ convention
        self.nrm = nrm            # (N,3) or None
        self.tris = tris          # (M,3) int

    @property
    def xz(self):
        return self.pos[:, [0, 2]]


def load_obj(path):
    V, VT, VN, faces = [], [], [], []
    with open(path, "r") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag == "v":
                V.append([float(x) for x in parts[1:4]])
            elif tag == "vt":
                VT.append([float(x) for x in parts[1:3]])
            elif tag == "vn":
                VN.append([float(x) for x in parts[1:4]])
            elif tag == "f":
                corners = []
                for tok in parts[1:]:
                    bits = (tok.split("/") + ["", ""])[:3]

                    def idx(s, n):
                        if s == "":
                            return -1
                        i = int(s)
                        return i - 1 if i > 0 else n + i
                    corners.append((idx(bits[0], len(V)),
                                    idx(bits[1], len(VT)),
                                    idx(bits[2], len(VN))))
                for k in range(1, len(corners) - 1):      # fan-triangulate
                    faces.append((corners[0], corners[k], corners[k + 1]))
    if not V:
        raise SystemExit(f"{path}: no vertices")
    if not VT:
        raise SystemExit(f"{path}: no texture coordinates (vt) -- required")
    V = np.asarray(V, float)
    VT = np.asarray(VT, float)
    VN = np.asarray(VN, float) if VN else None

    remap, pos, uv, nrm, tris = {}, [], [], [], []
    for tri in faces:
        out = []
        for c in tri:
            key = c
            j = remap.get(key)
            if j is None:
                j = len(pos)
                remap[key] = j
                pos.append(V[c[0]])
                uv.append(VT[c[1]] if c[1] >= 0 else [0.0, 0.0])
                if VN is not None:
                    nrm.append(VN[c[2]] if c[2] >= 0 else [0.0, 1.0, 0.0])
            out.append(j)
        tris.append(out)
    return Mesh(np.asarray(pos, float), np.asarray(uv, float),
                np.asarray(nrm, float) if nrm else None,
                np.asarray(tris, int))


def write_obj(path, pos, uv, nrm, tris, header=""):
    out = []
    if header:
        out.append(f"# {header}")
    out.append(f"# {len(pos)} vertices, {len(tris)} triangles")
    out.extend("v %.6f %.6f %.6f" % tuple(p) for p in pos)
    out.extend("vt %.8f %.8f" % tuple(t) for t in uv)
    if nrm is not None:
        out.extend("vn %.6f %.6f %.6f" % tuple(n) for n in nrm)
        out.extend("f %d/%d/%d %d/%d/%d %d/%d/%d" %
                   (a + 1, a + 1, a + 1, b + 1, b + 1, b + 1, c + 1, c + 1, c + 1)
                   for a, b, c in tris)
    else:
        out.extend("f %d/%d %d/%d %d/%d" %
                   (a + 1, a + 1, b + 1, b + 1, c + 1, c + 1)
                   for a, b, c in tris)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


# --------------------------------------------------------------------------
# 2. image -> skeleton -> polylines
# --------------------------------------------------------------------------

def _nbrs(p):
    """8 neighbours of a padded boolean image, in Zhang-Suen order P2..P9."""
    return (p[0:-2, 1:-1], p[0:-2, 2:], p[1:-1, 2:], p[2:, 2:],
            p[2:, 1:-1], p[2:, 0:-2], p[1:-1, 0:-2], p[0:-2, 0:-2])


def smooth_mask(mask, r):
    """Morphological close-then-open with a disk.

    Thinning is exquisitely sensitive to the roughness of a stroke's edge: a
    lumpy outline turns one stroke into a braid of short junction-to-junction
    fragments rather than a single centreline.  Rounding the outline first
    costs nothing and makes the trace clean.
    """
    if r <= 0:
        return mask
    from scipy.ndimage import binary_dilation, binary_erosion
    g = np.arange(-r, r + 1)
    st = (g[:, None] ** 2 + g[None, :] ** 2) <= r * r
    # border_value=1 on the erosions: a stroke that runs off the edge of the
    # image must stay attached to it, otherwise its endpoint retreats inside
    # the frame and the cut silently stops separating anything.
    m = binary_dilation(mask, st)
    m = binary_erosion(m, st, border_value=1)
    m = binary_erosion(m, st, border_value=1)
    return binary_dilation(m, st)


def skeletonize(mask):
    """Zhang-Suen thinning, fully vectorised.  Operates on the mask bbox."""
    ys, xs = np.nonzero(mask)
    if len(ys) == 0:
        return mask.copy()
    y0, y1 = ys.min(), ys.max() + 1
    x0, x1 = xs.min(), xs.max() + 1
    m = mask[y0:y1, x0:x1].copy()
    it = 0
    while True:
        changed = False
        for step in (0, 1):
            p = np.pad(m, 1)
            N, NE, E, SE, S, SW, W, NW = _nbrs(p)
            seq = (N, NE, E, SE, S, SW, W, NW)
            B = np.zeros(m.shape, np.uint8)
            A = np.zeros(m.shape, np.uint8)
            for i in range(8):
                B += seq[i]
                A += (~seq[i] & seq[(i + 1) % 8]).view(np.uint8)
            if step == 0:
                ok = ~(N & E & S) & ~(E & S & W)
            else:
                ok = ~(N & E & W) & ~(N & S & W)
            rm = m & (A == 1) & (B >= 2) & (B <= 6) & ok
            if rm.any():
                m &= ~rm
                changed = True
        it += 1
        if not changed or it > 500:
            break
    out = np.zeros_like(mask)
    out[y0:y1, x0:x1] = m
    log(f"thinned in {it} iterations: {int(mask.sum())} -> {int(out.sum())} px")
    return out


_OFF = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]


def trace_skeleton(skel):
    """Trace a 1px skeleton into polylines of (x, y) pixel coordinates."""
    H, W = skel.shape
    ys, xs = np.nonzero(skel)
    pix = set(zip(ys.tolist(), xs.tolist()))
    nbr = {}
    for (y, x) in pix:
        n = []
        for dy, dx in _OFF:
            q = (y + dy, x + dx)
            if q not in pix:
                continue
            # drop a diagonal link that merely short-cuts an orthogonal path,
            # otherwise a staircase pixel looks like a junction
            if dy and dx and ((y + dy, x) in pix or (y, x + dx) in pix):
                continue
            n.append(q)
        nbr[(y, x)] = n
    deg = {k: len(v) for k, v in nbr.items()}
    nodes = {k for k, d in deg.items() if d != 2}

    polys, used = [], set()

    def walk(start, nxt):
        path = [start, nxt]
        prev, cur = start, nxt
        while deg[cur] == 2:
            a, b = nbr[cur]
            nx = a if a != prev else b
            prev, cur = cur, nx
            path.append(cur)
            if cur == start:                       # closed loop back to start
                break
        return path

    for s in nodes:
        for n in nbr[s]:
            if (s, n) in used:
                continue
            path = walk(s, n)
            for i in range(len(path) - 1):
                used.add((path[i], path[i + 1]))
                used.add((path[i + 1], path[i]))
            polys.append(path)
    # isolated cycles (every pixel degree 2)
    seen = set()
    for p in polys:
        seen.update(p)
    for s in pix - seen:
        if deg[s] != 2:
            continue
        path = walk(s, nbr[s][0])
        if path[-1] != s:
            path.append(s)
        seen.update(path)
        polys.append(path)

    out = []
    for p in polys:
        arr = np.array([[x, y] for (y, x) in p], float)
        if len(arr) >= 2:
            out.append(arr)
    return out


def prune_spurs(polys, min_len):
    """Drop short branches that dead-end (thinning artefacts)."""
    if min_len <= 0:
        return polys
    changed = True
    while changed:
        changed = False
        ends = defaultdict(int)
        for p in polys:
            ends[tuple(p[0])] += 1
            ends[tuple(p[-1])] += 1
        keep = []
        for p in polys:
            length = np.linalg.norm(np.diff(p, axis=0), axis=1).sum()
            dangling = ends[tuple(p[0])] == 1 or ends[tuple(p[-1])] == 1
            closed = tuple(p[0]) == tuple(p[-1])
            if dangling and not closed and length < min_len:
                changed = True
                continue
            keep.append(p)
        polys = keep
    return polys


def merge_chains(polys):
    """Join polylines end-to-end at nodes where exactly two of them meet.

    Thinning turns one drawn stroke into many fragments (every little nub is
    a junction).  After pruning the nubs, the leftovers are a chain -- merging
    them back lets Douglas-Peucker simplify across the whole stroke instead of
    preserving a vertex at every former junction.
    """
    def k(v):
        return tuple(np.round(v, 6))

    for _ in range(64):
        ends = defaultdict(list)
        for i, p in enumerate(polys):
            if k(p[0]) == k(p[-1]):
                continue                       # already a closed loop
            ends[k(p[0])].append((i, 0))
            ends[k(p[-1])].append((i, 1))
        used, merged = set(), []
        for node, v in ends.items():
            if len(v) != 2:
                continue
            (i, e1), (j, e2) = v
            if i == j or i in used or j in used:
                continue
            a = polys[i] if e1 == 1 else polys[i][::-1]     # a ends at node
            b = polys[j] if e2 == 0 else polys[j][::-1]     # b starts at node
            merged.append(np.vstack([a, b[1:]]))
            used.add(i)
            used.add(j)
        if not merged:
            break
        polys = merged + [p for i, p in enumerate(polys) if i not in used]
    return polys


# --------------------------------------------------------------------------
# 3. simplify / snap
# --------------------------------------------------------------------------

def rdp(pts, eps):
    """Douglas-Peucker, iterative."""
    n = len(pts)
    if n < 3 or eps <= 0:
        return pts
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = pts[i], pts[j]
        ab = b - a
        L = math.hypot(*ab)
        seg = pts[i + 1:j]
        if L < 1e-12:
            d = np.linalg.norm(seg - a, axis=1)
        else:
            d = np.abs(cross2(np.broadcast_to(ab, seg.shape), seg - a)) / L
        k = int(np.argmax(d))
        if d[k] > eps:
            k += i + 1
            keep[k] = True
            stack.append((i, k))
            stack.append((k, j))
    return pts[keep]


def _closest_on_segments(p, A, B):
    """Closest point to p on each segment A[i]->B[i].  Returns (pts, t, dist)."""
    d = B - A
    L2 = (d * d).sum(1)
    t = np.where(L2 > 1e-18, ((p - A) * d).sum(1) / np.maximum(L2, 1e-18), 0.0)
    t = np.clip(t, 0.0, 1.0)
    q = A + t[:, None] * d
    return q, t, np.linalg.norm(q - p, axis=1)


def snap_endpoints(polys, tol, targets=()):
    """Weld dangling endpoints onto a nearby polyline, splitting the target.

    `targets` are extra polylines (the mesh rim) that an endpoint may snap
    *onto* but that are never split or returned -- the rim gets its shared
    vertex inserted later, where it is clipped.
    """
    polys = [p.copy() for p in polys]
    targets = [np.asarray(t, float) for t in targets]

    def touch_count(polys):
        """How many polylines have a vertex at each endpoint location."""
        cnt = defaultdict(int)
        for p in polys:
            for v in np.round(p, 6):
                cnt[tuple(v)] += 1
        return cnt

    for _round in range(8):
        touches = touch_count(polys)
        moved = False
        for pi, p in enumerate(polys):
            for which in (0, -1):
                key = tuple(np.round(p[which], 6))
                if touches[key] > 1 or tuple(p[0]) == tuple(p[-1]):
                    continue                      # already joined / closed
                pt = p[which]
                cands = [(qi, q) for qi, q in enumerate(polys) if qi != pi]
                cands += [(None, q) for q in targets]
                best = None
                for qi, q in cands:
                    _, t, dist = _closest_on_segments(pt, q[:-1], q[1:])
                    k = int(np.argmin(dist))
                    if dist[k] <= tol and (best is None or dist[k] < best[0]):
                        best = (dist[k], qi, k, t[k], q)
                if best is None:
                    continue
                dist0, qi, k, t, q = best
                target = q[k] + t * (q[k + 1] - q[k])
                if dist0 <= 1e-9:
                    continue                      # already exactly on it
                p[which] = target
                if qi is not None and 1e-9 < t < 1 - 1e-9:
                    polys[qi] = np.vstack([q[:k + 1], target[None], q[k + 1:]])
                moved = True
        if not moved:
            break
    return polys


def project_to_loops(pt, loops):
    """Nearest point on the rings: (ring index, edge index, s, distance)."""
    best = (0, 0, 0.0, np.inf)
    for ri, ring in enumerate(loops):
        A = ring
        B = np.roll(ring, -1, axis=0)
        _, t, d = _closest_on_segments(pt, A, B)
        k = int(np.argmin(d))
        if d[k] < best[3]:
            best = (ri, k, float(t[k]), float(d[k]))
    return best


# --------------------------------------------------------------------------
# 4. UV <-> pixels, mesh locator
# --------------------------------------------------------------------------

class UVImage:
    """Bijection between OBJ uv coordinates and pixel coordinates."""

    def __init__(self, w, h, flip_v=True):
        self.w, self.h, self.flip_v = w, h, flip_v

    def uv_to_px(self, uv):
        x = uv[..., 0] * (self.w - 1)
        v = uv[..., 1]
        y = (1.0 - v) * (self.h - 1) if self.flip_v else v * (self.h - 1)
        return np.stack([x, y], -1)

    def px_to_uv(self, px):
        u = px[..., 0] / (self.w - 1)
        y = px[..., 1] / (self.h - 1)
        v = 1.0 - y if self.flip_v else y
        return np.stack([u, v], -1)


class TriLocator:
    """Point location in a 2D triangulation, via a KD-tree over centroids."""

    def __init__(self, pts2, tris):
        self.p = pts2
        self.t = tris
        tp = pts2[tris]
        self.cen = tp.mean(1)
        self.rad = np.linalg.norm(tp - self.cen[:, None], axis=2).max(1)
        self.maxrad = float(self.rad.max())
        self.tree = KDTree(self.cen)
        a, b, c = tp[:, 0], tp[:, 1], tp[:, 2]
        self.a = a
        m = np.stack([b - a, c - a], -1)            # (M,2,2)
        det = m[:, 0, 0] * m[:, 1, 1] - m[:, 0, 1] * m[:, 1, 0]
        det = np.where(np.abs(det) < 1e-18, 1e-18, det)
        inv = np.empty_like(m)
        inv[:, 0, 0] = m[:, 1, 1] / det
        inv[:, 0, 1] = -m[:, 0, 1] / det
        inv[:, 1, 0] = -m[:, 1, 0] / det
        inv[:, 1, 1] = m[:, 0, 0] / det
        self.inv = inv

    def bary(self, tidx, pts):
        d = pts - self.a[tidx]
        inv = self.inv[tidx]
        l1 = inv[:, 0, 0] * d[:, 0] + inv[:, 0, 1] * d[:, 1]
        l2 = inv[:, 1, 0] * d[:, 0] + inv[:, 1, 1] * d[:, 1]
        return np.stack([1 - l1 - l2, l1, l2], -1)

    def locate(self, pts, tol=1e-9):
        """Returns (tri index or -1, barycentric coords) for each point."""
        pts = np.atleast_2d(np.asarray(pts, float))
        n = len(pts)
        out = np.full(n, -1)
        bar = np.zeros((n, 3))
        todo = np.arange(n)
        for k in (8, 32, 128):
            if len(todo) == 0:
                break
            k = min(k, len(self.cen))
            _, cand = self.tree.query(pts[todo], k=k)
            cand = np.atleast_2d(cand)
            for col in range(cand.shape[1]):
                if len(todo) == 0:
                    break
                t = cand[:, col]
                bb = self.bary(t, pts[todo])
                hit = (bb >= -tol).all(1)
                if hit.any():
                    out[todo[hit]] = t[hit]
                    bar[todo[hit]] = bb[hit]
                    keep = ~hit
                    todo, cand = todo[keep], cand[keep]
            if len(todo):
                # guaranteed search: a triangle containing p has its centroid
                # within maxrad of p
                still = []
                for i in todo:
                    cs = self.tree.query_ball_point(pts[i], self.maxrad)
                    if not cs:
                        still.append(i)
                        continue
                    cs = np.asarray(cs)
                    bb = self.bary(cs, np.repeat(pts[i][None], len(cs), 0))
                    j = int(np.argmax(bb.min(1)))
                    if bb[j].min() >= -tol:
                        out[i], bar[i] = cs[j], bb[j]
                    else:
                        still.append(i)
                todo = np.asarray(still, int)
            break
        return out, bar

    def locate_clamped(self, pts):
        """Like locate() but never fails: falls back to the nearest triangle."""
        t, b = self.locate(pts)
        bad = t < 0
        if bad.any():
            idx = np.nonzero(bad)[0]
            for i in idx:
                cs = np.asarray(self.tree.query(pts[i], k=min(64, len(self.cen)))[1]).ravel()
                bb = self.bary(cs, np.repeat(pts[i][None], len(cs), 0))
                j = int(np.argmax(bb.min(1)))
                t[i] = cs[j]
                b[i] = np.clip(bb[j], 0, None)
                b[i] /= b[i].sum()
        return t, b

    def candidates_near_segment(self, a, b):
        mid = (a + b) * 0.5
        r = float(np.linalg.norm(b - a)) * 0.5 + self.maxrad
        return np.asarray(self.tree.query_ball_point(mid, r), int)


def boundary_loops(tris, nv):
    """Ordered boundary rings (vertex indices) of a manifold triangle soup."""
    cnt = defaultdict(int)
    direction = {}
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            cnt[(min(u, v), max(u, v))] += 1
            direction[(min(u, v), max(u, v))] = (u, v)
    nxt = {}
    for e, n in cnt.items():
        if n == 1:
            u, v = direction[e]
            nxt[u] = v
    loops = []
    seen = set()
    for s in list(nxt):
        if s in seen:
            continue
        loop, cur = [], s
        while cur not in seen:
            seen.add(cur)
            loop.append(cur)
            cur = nxt.get(cur)
            if cur is None:
                break
        if len(loop) >= 3:
            loops.append(np.asarray(loop, int))
    return loops


# --------------------------------------------------------------------------
# 5. clipping + lifting
# --------------------------------------------------------------------------

def clip_polyline(poly, loops, rim_eps=1e-6):
    """Clip a polyline (pixel space) to the interior of the given rings.

    Returns (pieces, splits) where splits is {(ring, edge): [t, ...]} marking
    where the rings must be subdivided so the constraints share vertices.
    """
    pieces, splits = [], defaultdict(list)
    cur = []
    for i in range(len(poly) - 1):
        a, b = poly[i], poly[i + 1]
        ts = [0.0, 1.0]
        for ri, ring in enumerate(loops):
            A = ring
            B = np.roll(ring, -1, axis=0)
            d1 = b - a
            d2 = B - A
            den = cross2(np.broadcast_to(d1, d2.shape), d2)
            ok = np.abs(den) > 1e-12
            if not ok.any():
                continue
            w = a - A
            t = np.where(ok, cross2(d2, w) / np.where(ok, den, 1), -1.0)
            s = np.where(ok, cross2(np.broadcast_to(d1, d2.shape), w) / np.where(ok, den, 1), -1.0)
            hit = ok & (t > 1e-12) & (t < 1 - 1e-12) & (s >= -1e-12) & (s <= 1 + 1e-12)
            for ei in np.nonzero(hit)[0]:
                ts.append(float(t[ei]))
                splits[(ri, int(ei))].append(float(np.clip(s[ei], 0, 1)))
        ts = sorted(set(round(t, 12) for t in ts))
        pts = [a + t * (b - a) for t in ts]
        mids = np.array([(pts[k] + pts[k + 1]) * 0.5 for k in range(len(pts) - 1)])
        # strictly inside: a piece lying along the rim would duplicate the
        # mesh-boundary constraint, and overlapping constraints are
        # unrecoverable
        ins = points_in_loops(mids, loops) & (dist_to_loops(mids, loops) > rim_eps)
        for k, keep in enumerate(ins):
            if keep:
                if not cur:
                    cur = [pts[k]]
                cur.append(pts[k + 1])
            else:
                if len(cur) >= 2:
                    pieces.append(np.array(cur))
                cur = []
    if len(cur) >= 2:
        pieces.append(np.array(cur))
    return pieces, splits


def split_at_triangle_edges(poly, loc):
    """Subdivide a polyline so every segment lies inside a single triangle."""
    out = [poly[0]]
    for i in range(len(poly) - 1):
        a, b = poly[i], poly[i + 1]
        cand = loc.candidates_near_segment(a, b)
        ts = []
        if len(cand):
            tp = loc.p[loc.t[cand]]                     # (C,3,2)
            A = tp.reshape(-1, 2)
            B = np.roll(tp, -1, axis=1).reshape(-1, 2)
            d1 = b - a
            d2 = B - A
            den = cross2(np.broadcast_to(d1, d2.shape), d2)
            ok = np.abs(den) > 1e-12
            w = a - A
            t = np.where(ok, cross2(d2, w) / np.where(ok, den, 1), -1.0)
            s = np.where(ok, cross2(np.broadcast_to(d1, d2.shape), w) / np.where(ok, den, 1), -1.0)
            hit = ok & (t > 1e-9) & (t < 1 - 1e-9) & (s >= -1e-9) & (s <= 1 + 1e-9)
            ts = sorted(set(np.round(t[hit], 10).tolist()))
        for tv in ts:
            out.append(a + tv * (b - a))
        out.append(b)
    arr = np.array(out)
    d = np.linalg.norm(np.diff(arr, axis=0), axis=1)
    keep = np.concatenate([[True], d > 1e-9])
    return arr[keep]


def dist_to_loops(pts, loops):
    """Distance from each point to the nearest ring edge."""
    out = np.full(len(pts), np.inf)
    for ring in loops:
        A = ring
        B = np.roll(ring, -1, axis=0)
        d = B - A
        L2 = np.maximum((d * d).sum(1), 1e-18)
        for i in range(len(A)):
            t = np.clip(((pts - A[i]) @ d[i]) / L2[i], 0.0, 1.0)
            out = np.minimum(out, np.linalg.norm(A[i] + t[:, None] * d[i] - pts,
                                                 axis=1))
    return out


def resolve_crossings(chains, max_passes=6):
    """Split chains wherever two of them cross, so constraints only ever meet
    at shared vertices.

    A constrained Delaunay cannot contain two crossing edges, so a crossing
    pair makes recovery impossible -- midpoint splitting then doubles the
    offenders every round instead of converging.  Building the arrangement up
    front removes the failure mode entirely.
    """
    chains = [np.asarray(c, float) for c in chains]
    for _ in range(max_passes):
        segs, owner = [], []
        for ci, c in enumerate(chains):
            for i in range(len(c) - 1):
                segs.append((c[i], c[i + 1]))
                owner.append((ci, i))
        if not segs:
            break
        A = np.array([s[0] for s in segs])
        B = np.array([s[1] for s in segs])
        lo = np.minimum(A, B)
        hi = np.maximum(A, B)
        mn, mx = lo.min(0), hi.max(0)
        span = np.maximum(mx - mn, 1e-12)
        n = int(np.clip(math.sqrt(len(segs)), 1, 128))
        c0 = np.floor((lo - mn) / span * n).astype(int).clip(0, n - 1)
        c1 = np.floor((hi - mn) / span * n).astype(int).clip(0, n - 1)
        buckets = defaultdict(list)
        for k in range(len(segs)):
            for gx in range(c0[k, 0], c1[k, 0] + 1):
                for gy in range(c0[k, 1], c1[k, 1] + 1):
                    buckets[(gx, gy)].append(k)

        splits = defaultdict(list)
        seen = set()
        col_tol = float(np.linalg.norm(mx - mn)) * 1e-9
        e_ov = 1e-9
        for cell in buckets.values():
            for ii in range(len(cell)):
                for jj in range(ii + 1, len(cell)):
                    k, m = cell[ii], cell[jj]
                    if (k, m) in seen:
                        continue
                    seen.add((k, m))
                    if owner[k][0] == owner[m][0] and abs(owner[k][1] - owner[m][1]) <= 1:
                        continue                      # adjacent in one chain
                    d1, d2 = B[k] - A[k], B[m] - A[m]
                    den = float(cross2(d1, d2))
                    if abs(den) < 1e-15:
                        # parallel -- if they are also collinear and overlap,
                        # split both at the overlap so the shared stretch
                        # becomes one identical (de-duplicated) edge
                        L1 = float(np.linalg.norm(d1))
                        if L1 < 1e-15:
                            continue
                        if abs(float(cross2(d1, A[m] - A[k]))) / L1 > col_tol:
                            continue
                        u = d1 / L1
                        for a_pt, b_k in ((A[m], k), (B[m], k)):
                            tt = float((a_pt - A[k]) @ u) / L1
                            if e_ov < tt < 1 - e_ov:
                                splits[b_k].append(tt)
                        L2 = float(np.linalg.norm(d2))
                        if L2 > 1e-15:
                            u2 = d2 / L2
                            for a_pt in (A[k], B[k]):
                                tt = float((a_pt - A[m]) @ u2) / L2
                                if e_ov < tt < 1 - e_ov:
                                    splits[m].append(tt)
                        continue
                    w = A[k] - A[m]
                    t = float(cross2(d2, w)) / den
                    s = float(cross2(d1, w)) / den
                    e = 1e-9
                    if e < t < 1 - e and e < s < 1 - e:
                        splits[k].append(t)
                        splits[m].append(s)
        if not splits:
            break
        log(f"  resolved {len(splits)} crossing constraint segment(s)")
        per_chain = defaultdict(dict)
        for k, ts in splits.items():
            ci, i = owner[k]
            per_chain[ci][i] = sorted(set(np.round(ts, 12)))
        for ci, m in per_chain.items():
            c = chains[ci]
            out = [c[0]]
            for i in range(len(c) - 1):
                for t in m.get(i, []):
                    out.append(c[i] + t * (c[i + 1] - c[i]))
                out.append(c[i + 1])
            chains[ci] = np.array(out)
    return chains


def weld_px(chains, tol):
    """Make coincident pixel-space vertices bit-identical across chains.

    Where two cut lines meet, both chains carry a vertex at the same pixel.
    Lifting them separately can route them through different triangles and
    return world positions that differ in the last few digits -- enough to
    split the junction into a zero-area sliver with a 100000:1 aspect ratio.
    Collapsing them to one representative first makes the junction exact.
    """
    if not chains:
        return chains
    allp = np.vstack(chains)
    tree = KDTree(allp)
    parent = list(range(len(allp)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    n = 0
    for a, b in tree.query_pairs(max(tol, 0.0)):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
            n += 1
    if n:
        log(f"  welded {n} coincident constraint vertex pair(s) "
            f"(within {tol:g} px)")
    rep = np.array(allp)
    for i in range(len(allp)):
        rep[i] = allp[find(i)]
    out, k = [], 0
    for c in chains:
        nc = rep[k:k + len(c)]
        k += len(c)
        keep = np.concatenate([[True],
                               (np.diff(nc, axis=0) != 0).any(1)])
        nc = nc[keep]
        if len(nc) >= 2:
            out.append(nc)
    return out


def densify(pts, max_d):
    """Insert midpoints until no consecutive pair is further apart than max_d."""
    if max_d <= 0:
        return pts
    cur = pts
    for _ in range(24):
        d = np.linalg.norm(np.diff(cur, axis=0), axis=1)
        bad = d > max_d
        if not bad.any():
            break
        out = [cur[0]]
        for i in range(len(cur) - 1):
            if bad[i]:
                out.append((cur[i] + cur[i + 1]) * 0.5)
            out.append(cur[i + 1])
        cur = np.array(out)
    return cur


# --------------------------------------------------------------------------
# 6. constrained Delaunay (pure scipy) + Ruppert refinement
# --------------------------------------------------------------------------

class PointSet:
    def __init__(self, pts, quantum):
        self.q = quantum
        self.pts = []
        self.index = {}
        self.add_many(pts)

    def key(self, p):
        return (int(round(p[0] / self.q)), int(round(p[1] / self.q)))

    def add(self, p):
        k = self.key(p)
        j = self.index.get(k)
        if j is None:
            j = len(self.pts)
            self.index[k] = j
            self.pts.append(np.asarray(p, float))
        return j

    def add_many(self, P):
        return [self.add(p) for p in np.atleast_2d(P)]

    def array(self):
        return np.asarray(self.pts, float)

    def snapshot(self):
        return len(self.pts)

    def restore(self, n):
        """Drop every point added after the snapshot."""
        if n >= len(self.pts):
            return
        del self.pts[n:]
        self.index = {k: v for k, v in self.index.items() if v < n}


def tri_edges(simplices):
    e = np.concatenate([simplices[:, [0, 1]], simplices[:, [1, 2]], simplices[:, [2, 0]]])
    e.sort(axis=1)
    return e


def split_constraint(ps, segs, i, j, encroachers):
    """Split constraint (i,j).  Returns True if the point set actually grew.

    Prefers splitting at an encroaching vertex that lies *on* the segment --
    midpoint bisection can never resolve that case and would loop forever.
    """
    P = ps.pts
    a, b = P[i], P[j]
    d = b - a
    L2 = float(d @ d)
    k = None
    if L2 > 0:
        tol = ps.q * 8.0
        best = tol
        for c in encroachers:
            if c == i or c == j:
                continue
            p = P[c]
            t = float((p - a) @ d) / L2
            if not (1e-9 < t < 1 - 1e-9):
                continue
            perp = abs(float(cross2(d, p - a))) / math.sqrt(L2)
            if perp <= best:
                best, k = perp, c
    if k is None:
        k = ps.add((a + b) * 0.5)
    if k == i or k == j:
        return False
    segs.discard((i, j))
    segs.add((min(i, k), max(i, k)))
    segs.add((min(j, k), max(j, k)))
    return True


def constrained_delaunay(ps, segs, max_rounds=200):
    """Recover `segs` as Delaunay edges by splitting segments.

    Two mechanisms, alternated:
      * encroachment -- a segment is guaranteed to be Delaunay once no other
        vertex lies in its diametral circle, so split encroached segments;
      * explicit repair -- that guarantee only holds for points in general
        position.  Cocircular input (a regular grid is the obvious case) lets
        qhull pick any of several valid triangulations, and it may simply not
        choose ours, so we also split any constraint that is still missing
        from the actual triangulation.
    """
    segs = set(segs)
    tri = Delaunay(ps.array())
    prev_missing = 0
    for rnd in range(max_rounds):
        pts = ps.array()
        tree = KDTree(pts)
        slist = list(segs)
        A = pts[[s[0] for s in slist]]
        B = pts[[s[1] for s in slist]]
        mid = (A + B) * 0.5
        rad = np.linalg.norm(B - A, axis=1) * 0.5
        hits = tree.query_ball_point(mid, rad * (1 - 1e-9))
        nsplit = stuck = 0
        for (i, j), h in zip(slist, hits):
            enc = [k for k in h if k != i and k != j]
            if not enc:
                continue
            if split_constraint(ps, segs, i, j, enc):
                nsplit += 1
            else:
                stuck += 1

        tri = Delaunay(ps.array())
        have = set(map(tuple, tri_edges(tri.simplices)))
        missing = [s for s in segs if s not in have]
        for (i, j) in missing:
            if split_constraint(ps, segs, i, j, []):
                nsplit += 1
            else:
                stuck += 1

        if nsplit == 0:
            if stuck:
                log(f"  WARNING: {stuck} constraint segment(s) could not be "
                    f"split further (duplicate/overlapping input geometry)")
            break
        if len(missing) > max(8, 2 * prev_missing) and prev_missing:
            print(f"ERROR: constraint recovery is diverging "
                  f"({prev_missing} -> {len(missing)} unrecovered edges). "
                  f"This means two cut lines overlap or cross without "
                  f"sharing a vertex; splitting them can never converge.",
                  file=sys.stderr)
            break
        prev_missing = len(missing)
        if rnd % 10 == 0 or len(missing):
            log(f"  CDT round {rnd}: {nsplit} split "
                f"({len(missing)} missing), {len(ps.pts)} pts")
    else:
        log(f"  WARNING: constraint recovery hit the round cap")

    tri = Delaunay(ps.array())
    have = set(map(tuple, tri_edges(tri.simplices)))
    missing = [s for s in segs if s not in have]
    if missing:
        log(f"  WARNING: {len(missing)} constraint edges not recovered -- "
            f"cuts may leak between regions")
    else:
        log(f"  all {len(segs)} constraint edges recovered")
    return tri, segs


def circumcentre(a, b, c):
    ax, ay = a[:, 0], a[:, 1]
    bx, by = b[:, 0], b[:, 1]
    cx, cy = c[:, 0], c[:, 1]
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    d = np.where(np.abs(d) < 1e-18, 1e-18, d)
    a2, b2, c2 = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return np.stack([ux, uy], -1)


def tri_quality(pts, simp):
    """Returns (edge-ratio, radius/shortest-edge, circumcentre) per triangle."""
    a, b, c = pts[simp[:, 0]], pts[simp[:, 1]], pts[simp[:, 2]]
    la = np.linalg.norm(b - c, axis=1)
    lb = np.linalg.norm(c - a, axis=1)
    lc = np.linalg.norm(a - b, axis=1)
    L = np.stack([la, lb, lc], -1)
    lmin = np.maximum(L.min(1), 1e-18)
    lmax = L.max(1)
    area = np.abs(cross2(b - a, c - a)) * 0.5
    R = la * lb * lc / np.maximum(4 * area, 1e-18)
    return lmax / lmin, R / lmin, circumcentre(a, b, c)


def domain_mask(pts, simp, loops_xz):
    """Triangles that are inside the mesh footprint and not degenerate.

    Qhull emits zero-area facets wherever the input has collinear runs (our
    boundary and cut constraints are full of them); they are not real faces.
    """
    area = np.abs(cross2(pts[simp[:, 1]] - pts[simp[:, 0]],
                         pts[simp[:, 2]] - pts[simp[:, 0]])) * 0.5
    cen = pts[simp].mean(1)
    eps = max(area.max() * 1e-10, 1e-18)
    return points_in_loops(cen, loops_xz) & (area > eps), area


def refine(ps, segs, loops_xz, max_aspect, max_area, max_points=500_000,
           min_feature=0.0, max_rounds=60):
    """Ruppert-style refinement to bound triangle aspect ratio."""
    theta = math.degrees(math.asin(min(1.0, 1.0 / max(max_aspect, 1.0))))
    theta = min(theta, 29.0)                     # keep termination sane
    B_target = 1.0 / (2.0 * math.sin(math.radians(theta)))
    log(f"  refining to edge-ratio <= {max_aspect:g} "
        f"(min angle {theta:.1f}deg, radius-edge <= {B_target:.3f})")
    best_bad, worse = None, 0
    for rnd in range(max_rounds):
        pts = ps.array()
        tri = Delaunay(pts)
        simp = tri.simplices
        inside, area = domain_mask(pts, simp, loops_xz)
        er, re, cc = tri_quality(pts, simp)
        lmin = np.minimum.reduce([
            np.linalg.norm(pts[simp[:, 1]] - pts[simp[:, 2]], axis=1),
            np.linalg.norm(pts[simp[:, 2]] - pts[simp[:, 0]], axis=1),
            np.linalg.norm(pts[simp[:, 0]] - pts[simp[:, 1]], axis=1)])
        bad = inside & ((re > B_target) | (er > max_aspect))
        if max_area > 0:
            bad |= inside & (area > max_area)
        nbad = int(bad.sum())
        if nbad == 0:
            log(f"  refinement converged after {rnd} rounds "
                f"({int(inside.sum())} triangles)")
            break
        # Ruppert only terminates for well-conditioned input; small angles
        # between constraints can cascade forever.  Bail out loudly rather
        # than grinding the machine to a halt.
        if best_bad is None or nbad < best_bad:
            best_bad, worse = nbad, 0
        else:
            worse += 1
        if worse >= 3:
            print(f"WARNING: aspect-ratio refinement is not converging "
                  f"({nbad} triangles still violate --max-aspect "
                  f"{max_aspect:g}).  The cut is still correct; relax "
                  f"--max-aspect or use --no-refine.", file=sys.stderr)
            break
        if len(ps.pts) >= max_points:
            print(f"WARNING: hit --max-points ({max_points}); stopping "
                  f"refinement with {nbad} triangle(s) over --max-aspect.",
                  file=sys.stderr)
            break
        order = np.argsort(-np.where(bad, re, -np.inf))[:nbad]
        ok_cc = points_in_loops(cc[order], loops_xz)
        # thin the batch: each insertion must clear the local feature size, so
        # a batched round stays close to one-at-a-time Ruppert
        tree = KDTree(pts)
        keep, acc, acc_tree, acc_n = [], [], None, 0
        for oi, good in zip(order, ok_cc):
            if not good:
                continue
            p = cc[oi]
            # Ruppert inserts the circumcentre unconditionally; the only
            # reason to hold back is (a) the absolute feature floor, which
            # stops a small input angle from cascading into a pile of
            # coincident vertices, and (b) keeping one batch's insertions
            # apart from each other so batching stays close to the
            # one-at-a-time algorithm.
            # Batched insertion needs each new point to clear the local
            # feature size; dropping to the absolute floor alone lets
            # circumcentres land next to existing vertices and *creates*
            # slivers (measured: worst aspect 13 -> 154 on the self-test).
            s = max(lmin[oi] * 0.5, min_feature)
            if tree.query(p)[0] < s:
                continue
            if acc:
                if acc_tree is None or len(acc) - acc_n > 256:
                    acc_tree, acc_n = KDTree(np.asarray(acc)), len(acc)
                if acc_tree.query(p)[0] < s:
                    continue
                tail = acc[acc_n:]
                if tail and np.linalg.norm(np.asarray(tail) - p, axis=1).min() < s:
                    continue
            acc.append(p)
            keep.append(p)
            if len(keep) >= 20000:
                break
        if not keep:
            log(f"  refinement stopped: {int(bad.sum())} triangle(s) cannot be "
                f"improved without violating the local feature size")
            break
        # split any constraint segment encroached by a new point instead
        slist = list(segs)
        P = ps.array()
        A = P[[s[0] for s in slist]]
        Bp = P[[s[1] for s in slist]]
        mid = (A + Bp) * 0.5
        rad = np.linalg.norm(Bp - A, axis=1) * 0.5
        ktree = KDTree(np.asarray(keep))
        enc = ktree.query_ball_point(mid, rad)
        seglen = np.linalg.norm(Bp - A, axis=1)
        blocked = np.zeros(len(keep), bool)
        for (i, j), h, sl in zip(slist, enc, seglen):
            if not h:
                continue
            # splitting below the feature floor produces sub-segments the
            # triangulation can no longer recover -- that would corrupt the
            # cut itself, which matters far more than the aspect bound
            if sl <= 2.0 * min_feature:
                blocked[np.asarray(h, int)] = True
                continue
            split_constraint(ps, segs, i, j, [])
            blocked[np.asarray(h, int)] = True
        for p, bl in zip(keep, blocked):
            if not bl:
                ps.add(p)
        log(f"  refine round {rnd}: {int(bad.sum())} bad tris, "
            f"{len(ps.pts)} points")
    else:
        log("  refinement hit the round cap")
    return segs


# --------------------------------------------------------------------------
# 7. regions
# --------------------------------------------------------------------------

def label_regions(tri, segs, loops_xz):
    pts = tri.points
    simp = tri.simplices
    inside, _ = domain_mask(pts, simp, loops_xz)
    wall = set(segs)
    nb = tri.neighbors
    label = np.full(len(simp), -1)
    cur = 0
    for start in range(len(simp)):
        if not inside[start] or label[start] >= 0:
            continue
        stack = [start]
        label[start] = cur
        while stack:
            t = stack.pop()
            for e in range(3):
                n = nb[t, e]
                if n < 0 or not inside[n] or label[n] >= 0:
                    continue
                # the shared edge is the one opposite vertex e
                u, v = simp[t, (e + 1) % 3], simp[t, (e + 2) % 3]
                if (min(u, v), max(u, v)) in wall:
                    continue
                label[n] = cur
                stack.append(n)
        cur += 1
    return label, inside, cur


def merge_tiny_regions(tri, label, inside, nregions, min_area):
    """Fold negligible regions into their largest neighbour.

    Shallow crossings between cut lines genuinely do enclose slivers of a few
    thousandths of a unit.  They are not useful submeshes, and emitting them
    leaves the neighbouring pieces with unmatched boundary vertices.
    """
    if min_area <= 0 or nregions <= 1:
        return label, nregions
    pts, simp, nb = tri.points, tri.simplices, tri.neighbors
    area = np.abs(cross2(pts[simp[:, 1]] - pts[simp[:, 0]],
                         pts[simp[:, 2]] - pts[simp[:, 0]])) * 0.5
    for _ in range(nregions):
        tot = np.array([area[inside & (label == k)].sum()
                        for k in range(nregions)])
        alive = [k for k in range(nregions) if (inside & (label == k)).any()]
        small = [k for k in alive if tot[k] < min_area]
        if not small:
            break
        k = min(small, key=lambda z: tot[z])
        nbrs = set()
        for t in np.nonzero(inside & (label == k))[0]:
            for e in range(3):
                n = nb[t, e]
                if n >= 0 and inside[n] and label[n] != k:
                    nbrs.add(int(label[n]))
        if not nbrs:
            break
        label[(label == k) & inside] = max(nbrs, key=lambda z: tot[z])
    # compact the labels
    alive = sorted({int(x) for x in label[inside]})
    remap = {k: i for i, k in enumerate(alive)}
    out = label.copy()
    for k, i in remap.items():
        out[(label == k) & inside] = i
    return out, len(alive)


# --------------------------------------------------------------------------
# 8. debug render
# --------------------------------------------------------------------------

PALETTE = [(228, 94, 90), (90, 160, 228), (120, 196, 120), (238, 180, 80),
           (170, 120, 210), (90, 200, 200), (230, 130, 190), (160, 160, 110)]


def render_debug(path, mask, uvi, mesh, pts_xz, simp, label, inside, segs, loc):
    from PIL import ImageDraw
    base = Image.new("RGB", (uvi.w, uvi.h), (255, 255, 255))
    canvas = ImageDraw.Draw(base)
    t, b = loc.locate_clamped(pts_xz)
    uv = (mesh.uv[mesh.tris[t]] * b[:, :, None]).sum(1)
    px = uvi.uv_to_px(uv)
    for i in range(len(simp)):
        if not inside[i]:
            continue
        col = PALETTE[label[i] % len(PALETTE)]
        poly = [tuple(px[j]) for j in simp[i]]
        canvas.polygon(poly, fill=col, outline=(255, 255, 255))
    # the drawn strokes, as a translucent dark wash (their colour is arbitrary)
    over = np.zeros((uvi.h, uvi.w, 4), np.uint8)
    over[mask] = (20, 20, 20, 110)
    base = base.convert("RGBA")
    base.alpha_composite(Image.fromarray(over))
    canvas = ImageDraw.Draw(base)
    for (i, j) in segs:
        canvas.line([tuple(px[i]), tuple(px[j])], fill=(0, 0, 0), width=3)
    base.convert("RGB").resize((min(1400, uvi.w), min(1400, uvi.h)),
                               Image.LANCZOS).save(path)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(
        description="Cut an OBJ mesh into submeshes along lines drawn on a PNG.")
    ap.add_argument("mesh")
    ap.add_argument("image")
    ap.add_argument("-o", "--out-dir", default=".")
    ap.add_argument("--max-error", type=float, default=2.0,
                    help="max deviation (pixels) of the simplified cut from the "
                         "drawn line [2]")
    ap.add_argument("--max-distance", type=float, default=100.0,
                    help="max spacing between vertices along a cut; longer "
                         "spans are split at the midpoint [100]")
    ap.add_argument("--distance-units", choices=("world", "px"), default="world",
                    help="units for --max-distance [world]")
    ap.add_argument("--max-aspect", type=float, default=4.0,
                    help="max longest-edge/shortest-edge ratio of output "
                         "triangles (1 = equilateral) [4]")
    ap.add_argument("--max-area", type=float, default=0.0,
                    help="optional max triangle area, world units [off]")
    ap.add_argument("--weld", type=float, default=0.05,
                    help="collapse constraint vertices closer than this many "
                         "pixels into one [0.05]")
    ap.add_argument("--min-feature", type=float, default=-1.0,
                    help="refinement never places vertices closer together "
                         "than this, world units (-1 = auto, 0 = off)")
    ap.add_argument("--min-region-area", type=float, default=-1.0,
                    help="regions smaller than this are folded into their "
                         "largest neighbour (-1 = auto, 0 = keep all)")
    ap.add_argument("--max-points", type=int, default=500_000,
                    help="safety cap on vertices created by refinement [500k]")
    ap.add_argument("--no-refine", action="store_true",
                    help="skip quality refinement (ignore --max-aspect)")
    ap.add_argument("--snap", type=float, default=-1.0,
                    help="weld dangling line endpoints within this many pixels "
                         "of another line (-1 = auto, ~3x the stroke radius)")
    ap.add_argument("--alpha-threshold", type=int, default=0,
                    help="a pixel is part of a line if alpha > this [0]")
    ap.add_argument("--smooth", type=float, default=-1.0,
                    help="morphological smoothing radius (px) applied to the "
                         "line mask before thinning (-1 = auto, 0 = off)")
    ap.add_argument("--prune", type=float, default=-1.0,
                    help="drop thinning spurs shorter than N px (-1 = auto)")
    ap.add_argument("--flip-v", dest="flip_v", action="store_true", default=True,
                    help="uv (0,0) is bottom-left, the OBJ convention [default]")
    ap.add_argument("--no-flip-v", dest="flip_v", action="store_false",
                    help="uv (0,0) is top-left")
    ap.add_argument("--debug-png", default=None,
                    help="write an overlay of the regions over the input image")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    VERBOSE = not args.quiet

    os.makedirs(args.out_dir, exist_ok=True)
    _b, _e = os.path.splitext(os.path.basename(args.mesh))
    import glob as _glob
    stale = [f for f in _glob.glob(os.path.join(args.out_dir, f"{_b}_*{_e}"))
             if re.fullmatch(r"\d+", os.path.splitext(os.path.basename(f))[0][len(_b) + 1:] or "x")]
    for f in stale:
        os.remove(f)
    if stale:
        log(f"removed {len(stale)} stale output file(s)")

    # ---- load -------------------------------------------------------------
    mesh = load_obj(args.mesh)
    log(f"mesh: {len(mesh.pos)} vertices, {len(mesh.tris)} triangles")
    im = Image.open(args.image)
    arr = np.array(im.convert("RGBA"))
    H, W = arr.shape[:2]
    uvi = UVImage(W, H, args.flip_v)
    mask = arr[..., 3] > args.alpha_threshold
    log(f"image: {W}x{H}, {int(mask.sum())} line pixels")

    xz = mesh.xz
    uv_loc = TriLocator(uvi.uv_to_px(mesh.uv), mesh.tris)   # pixel space
    xz_loc = TriLocator(xz, mesh.tris)                       # world space

    # mesh footprint, in pixels and in world units
    loops_idx = boundary_loops(mesh.tris, len(mesh.pos))
    if not loops_idx:
        raise SystemExit("mesh has no boundary -- expected an open 2D manifold")
    loops_px = [uvi.uv_to_px(mesh.uv[l]) for l in loops_idx]
    log(f"mesh boundary: {len(loops_idx)} loop(s), "
        f"{sum(len(l) for l in loops_idx)} vertices")

    # ---- centrelines ------------------------------------------------------
    if mask.any():
        from scipy.ndimage import distance_transform_edt
        dt = distance_transform_edt(mask)
        stroke_r = float(np.median(dt[mask])) if mask.any() else 1.0
        smooth = args.smooth
        if smooth < 0:
            smooth = max(1, int(round(stroke_r / 3.0)))
            log(f"auto smoothing radius: {smooth} px "
                f"(stroke radius ~{stroke_r:.1f} px)")
        mask = smooth_mask(mask, int(smooth))
        skel = skeletonize(mask)
        if args.snap < 0:
            # an endpoint's position is only certain to within the stroke's
            # own half-width, so scale the weld tolerance to it
            args.snap = max(8.0, 3.0 * stroke_r)
            log(f"auto snap tolerance: {args.snap:.1f} px")
        polys = trace_skeleton(skel)
        log(f"traced {len(polys)} raw polylines")
        prune = args.prune
        if prune < 0:
            dt = distance_transform_edt(mask)
            prune = float(np.median(dt[skel])) * 3.0 if skel.any() else 0.0
            log(f"auto prune length: {prune:.1f} px")
        polys = prune_spurs(polys, prune)
        polys = merge_chains(polys)
        polys = [rdp(p, args.max_error) for p in polys]
        polys = [p for p in polys if len(p) >= 2]
        log(f"{len(polys)} polylines, {sum(len(p) for p in polys)} vertices "
            f"after pruning + simplification")
        rims = [np.vstack([r, r[:1]]) for r in loops_px]
        polys = snap_endpoints(polys, args.snap, targets=rims)
    else:
        polys = []
        if args.snap < 0:
            args.snap = 8.0

    # ---- clip to the mesh footprint, lift to world ------------------------
    constraints_px, ring_splits = [], defaultdict(list)
    for p in polys:
        pieces, sp = clip_polyline(p, loops_px)
        constraints_px.extend(pieces)
        for k, v in sp.items():
            ring_splits[k].extend(v)
    log(f"{len(constraints_px)} cut polylines inside the mesh footprint")

    # a cut that ends exactly on the rim needs the rim to carry that vertex too
    rim_tol = max(args.weld, 1e-6)
    for piece in constraints_px:
        for which in (0, -1):
            ri, ei, sv, dd = project_to_loops(piece[which], loops_px)
            if dd <= rim_tol and 1e-9 < sv < 1 - 1e-9:
                ring_splits[(ri, ei)].append(sv)

    # A cut only separates anything if both ends reach the mesh rim or another
    # cut.  That is only knowable after clipping: a stroke running off the
    # edge of the image looks dangling in pixel space but is perfectly fine.
    touch = defaultdict(int)
    for piece in constraints_px:
        for v in np.round(piece, 6):
            touch[tuple(v)] += 1
    loose = []
    for piece in constraints_px:
        if tuple(np.round(piece[0], 6)) == tuple(np.round(piece[-1], 6)):
            continue                                   # closed loop
        for which in (0, -1):
            pt = piece[which]
            if project_to_loops(pt, loops_px)[3] <= rim_tol:
                continue                               # terminates on the rim
            if touch[tuple(np.round(pt, 6))] > 1:
                continue                               # meets another cut
            loose.append(pt)
    if loose:
        print(f"WARNING: {len(loose)} dangling cut endpoint(s) separate "
              f"nothing (raise --snap to weld them):", file=sys.stderr)
        for q in loose[:12]:
            print(f"    pixel ({q[0]:.1f}, {q[1]:.1f})", file=sys.stderr)

    # the mesh boundary itself is a constraint, split where cuts leave it
    for ri, ring in enumerate(loops_px):
        pts = []
        for ei in range(len(ring)):
            a = ring[ei]
            b = ring[(ei + 1) % len(ring)]
            pts.append(a)
            for s in sorted(set(np.round(ring_splits[(ri, ei)], 12))):
                if 1e-9 < s < 1 - 1e-9:
                    pts.append(a + s * (b - a))
        pts.append(ring[0])
        constraints_px.append(np.array(pts))

    scale = np.linalg.norm(xz.max(0) - xz.min(0)) / max(
        np.linalg.norm(uvi.uv_to_px(mesh.uv).max(0) - uvi.uv_to_px(mesh.uv).min(0)), 1e-9)
    max_d = args.max_distance * (scale if args.distance_units == "px" else 1.0)

    split_px = [split_at_triangle_edges(p, uv_loc) for p in constraints_px]
    split_px = weld_px(split_px, args.weld)

    # lift each distinct pixel point exactly once, so a vertex shared by two
    # chains is bit-identical in world space
    order, seen = [], {}
    for c in split_px:
        for q in c:
            k = (q[0], q[1])
            if k not in seen:
                seen[k] = len(order)
                order.append(k)
    U = np.array(order, float)
    ti, bi = uv_loc.locate_clamped(U)
    WX = (mesh.pos[mesh.tris[ti]] * bi[:, :, None]).sum(1)[:, [0, 2]]

    chains_xz = []
    for c in split_px:
        world = WX[[seen[(q[0], q[1])] for q in c]]
        keep = np.concatenate([[True], (np.diff(world, axis=0) != 0).any(1)])
        world = densify(world[keep], max_d)
        if len(world) >= 2:
            chains_xz.append(world)
    extent = float(np.linalg.norm(xz.max(0) - xz.min(0)))
    min_feature = args.min_feature
    if min_feature < 0:
        min_feature = extent * 1e-4
    chains_xz = resolve_crossings(chains_xz)
    log(f"{sum(len(c) for c in chains_xz)} constraint vertices after densify")

    loops_xz = [xz[l] for l in loops_idx]

    # ---- build the constrained triangulation ------------------------------
    ps = PointSet(xz, quantum=extent * 1e-9)
    orig_ids = ps.add_many(xz)
    segs = set()
    for c in chains_xz:
        ids = ps.add_many(c)
        for i, j in zip(ids[:-1], ids[1:]):
            if i != j:
                segs.add((min(i, j), max(i, j)))
    log(f"CDT input: {len(ps.pts)} points, {len(segs)} constraint segments")
    tri, segs = constrained_delaunay(ps, segs)

    if not args.no_refine:
        mark, segs_good = ps.snapshot(), set(segs)
        segs = refine(ps, segs, loops_xz, args.max_aspect, args.max_area,
                      args.max_points, min_feature)
        tri = Delaunay(ps.array())
        have = set(map(tuple, tri_edges(tri.simplices)))
        missing = [s for s in segs if s not in have]
        if missing:
            log(f"  re-recovering {len(missing)} constraints after refinement")
            tri, segs = constrained_delaunay(ps, segs)
            have = set(map(tuple, tri_edges(tri.simplices)))
            missing = [s for s in segs if s not in have]
        if missing:
            # A correct cut beats a pretty one: throw the refinement away
            # rather than emit regions that leak into each other.
            print(f"WARNING: refinement left {len(missing)} constraint(s) "
                  f"unrecoverable; discarding it and falling back to the "
                  f"unrefined triangulation so the cut stays exact.  "
                  f"Relax --max-aspect or raise --min-feature.",
                  file=sys.stderr)
            ps.restore(mark)
            segs = segs_good
            tri, segs = constrained_delaunay(ps, segs)

    # ---- regions ----------------------------------------------------------
    label, inside, n = label_regions(tri, segs, loops_xz)
    min_region = args.min_region_area
    if min_region < 0:
        total = float(np.abs(cross2(xz[mesh.tris[:, 1]] - xz[mesh.tris[:, 0]],
                                    xz[mesh.tris[:, 2]] - xz[mesh.tris[:, 0]])
                             ).sum() * 0.5)
        min_region = total * 1e-6
    before = n
    label, n = merge_tiny_regions(tri, label, inside, n, min_region)
    if n != before:
        log(f"merged {before - n} negligible region(s) "
            f"(< {min_region:.4g} area) into their neighbours")
    pts2 = tri.points
    simp = tri.simplices
    areas = np.abs(cross2(pts2[simp[:, 1]] - pts2[simp[:, 0]],
                          pts2[simp[:, 2]] - pts2[simp[:, 0]])) * 0.5
    order = sorted(range(n), key=lambda k: -areas[(label == k) & inside].sum())
    rank = {k: i for i, k in enumerate(order)}
    log(f"{n} region(s); total area {areas[inside].sum():.1f}")

    # ---- interpolate attributes, write ------------------------------------
    t, b = xz_loc.locate_clamped(pts2)
    corner = mesh.tris[t]
    pos = (mesh.pos[corner] * b[:, :, None]).sum(1)
    uvs = (mesh.uv[corner] * b[:, :, None]).sum(1)
    nrm = None
    if mesh.nrm is not None:
        nrm = (mesh.nrm[corner] * b[:, :, None]).sum(1)
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12)
    # original vertices keep their exact attributes
    arr_pts = ps.array()
    for src, dst in enumerate(orig_ids):
        pos[dst] = mesh.pos[src]
        uvs[dst] = mesh.uv[src]
        if nrm is not None:
            nrm[dst] = mesh.nrm[src]

    # winding: match the input mesh's orientation in the XZ plane
    ta = xz[mesh.tris]
    sgn = np.sign(cross2(ta[:, 1] - ta[:, 0], ta[:, 2] - ta[:, 0]))
    want = 1.0 if sgn.sum() >= 0 else -1.0

    base, ext = os.path.splitext(os.path.basename(args.mesh))
    written = []
    for k in range(n):
        sel = (label == k) & inside
        f = simp[sel]
        if len(f) == 0:
            continue
        s = np.sign(cross2(pts2[f[:, 1]] - pts2[f[:, 0]], pts2[f[:, 2]] - pts2[f[:, 0]]))
        flip = s != want
        f = f.copy()
        f[flip] = f[flip][:, [0, 2, 1]]
        used, inv = np.unique(f, return_inverse=True)
        tri_local = inv.reshape(f.shape)
        idx = rank[k]
        path = os.path.join(args.out_dir, f"{base}_{idx}{ext or '.obj'}")
        write_obj(path, pos[used], uvs[used],
                  nrm[used] if nrm is not None else None, tri_local,
                  header=f"submesh {idx} of {base}{ext} "
                         f"({len(tri_local)} triangles)")
        ok, why = check_manifold(tri_local, len(used))
        area = areas[sel].sum()
        q = pos[used][:, [0, 2]]
        LL = np.stack([np.linalg.norm(q[tri_local[:, 1]] - q[tri_local[:, 2]], axis=1),
                       np.linalg.norm(q[tri_local[:, 2]] - q[tri_local[:, 0]], axis=1),
                       np.linalg.norm(q[tri_local[:, 0]] - q[tri_local[:, 1]], axis=1)], -1)
        ar = float((LL.max(1) / np.maximum(LL.min(1), 1e-18)).max())
        written.append((idx, path, len(used), len(tri_local), area, ok, why, ar))

    written.sort()
    print()
    print(f"{'file':<44}{'verts':>8}{'tris':>8}{'area':>14}{'aspect':>9}  manifold")
    for idx, path, nv, nt, area, ok, why, ar in written:
        print(f"{os.path.basename(path):<44}{nv:>8}{nt:>8}{area:>14.1f}"
              f"{ar:>9.2f}  {'yes' if ok else 'NO: ' + why}")
    worst = max((w[7] for w in written), default=0.0)
    print(f"{'TOTAL':<44}{'':>8}{sum(w[3] for w in written):>8}"
          f"{sum(w[4] for w in written):>14.1f}{worst:>9.2f}")
    if not args.no_refine and worst > args.max_aspect * 1.001:
        print(f"WARNING: worst triangle aspect ratio is {worst:.2f}, over the "
              f"requested --max-aspect {args.max_aspect:g}.  Some sliver of "
              f"the input geometry is finer than --min-feature "
              f"({min_feature:.4g}); see the refinement notes above.",
              file=sys.stderr)

    if args.debug_png:
        render_debug(args.debug_png, mask, uvi, mesh, pts2, simp, label, inside,
                     segs, xz_loc)
        log(f"wrote {args.debug_png}")
    return 0


def check_manifold(tris, nv):
    cnt = defaultdict(int)
    for a, b, c in tris:
        if a == b or b == c or a == c:
            return False, "degenerate face"
        for u, v in ((a, b), (b, c), (c, a)):
            cnt[(min(u, v), max(u, v))] += 1
    if any(v > 2 for v in cnt.values()):
        return False, "non-manifold edge"
    if len(set(np.asarray(tris).ravel().tolist())) != nv:
        return False, "unused vertex"
    return True, ""


if __name__ == "__main__":
    sys.exit(main())
