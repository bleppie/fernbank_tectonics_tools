#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.22",
#     "scipy>=1.10",
# ]
# ///
"""
Read an OBJ file, project its vertices onto the XZ plane (ignoring Y), compute
the 2D convex hull, print the distance between each pair of successive hull
vertices, and find the minimum-area bounding rectangle around the hull.

Input units are assumed to be millimeters; distances are printed in both
millimeters and feet.

OBJ files are also written alongside the input, all rotated so the bounding
rectangle is axis-aligned and centered on the origin:
  <name>_rotated_uvs.obj  the whole model, with UVs
  <name>_hull.obj         triangulated convex hull, all face normals pointing +Y
  <name>_no_back.obj      just the top surface: sides and back of the closed
                          print model removed (split at the sharp rim; see
                          --no-back-method)
  <name>_extended.obj     _no_back plus a skirt grown outward from its border
                          (sharing the border's vertices, so no seam)
                          by --extend-distance mm, following the terrain's
                          slope there (--extend-profile linear | ease). New
                          vertices share the original UV mapping, so their UVs
                          may fall outside 0..1.

Vertices at the same position are welded by default, so triangle-soup exports
shrink a lot. Corners whose normals differ by less than --smooth-angle (30 deg)
share one averaged normal (smooth shading); sharper creases stay split. Use
--smooth-angle 0 to keep faceted normals, or --no-weld to keep every vertex. Polygons with more
than three sides are fan-triangulated. Input UVs and vertex colours are
ignored: UVs are regenerated from the bounding rectangle.

Usage:
    python mesh_munger.py model.obj
    python mesh_munger.py model.obj --tolerance 5 --precision 1
    python mesh_munger.py model.obj --extend-distance 150 --extend-profile ease
    uv run mesh_munger.py model.obj          # installs numpy/scipy for you
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

MM_PER_FOOT = 304.8
WELD_EPS = 1e-4          # mm; positions closer than this are the same point
NORMAL_EPS = 1e-4        # normals closer than this (per component) are equal


# ---------------------------------------------------------------- mesh type

@dataclass
class Mesh:
    """A triangle mesh with optional per-vertex normals and UVs.

    `markers` are non-geometry OBJ lines (o, g, usemtl, s, mtllib) to emit
    just before face number `pos`, so groups/materials survive filtering.
    """
    v: np.ndarray                         # (n, 3) float positions
    faces: np.ndarray                     # (f, 3) int, 0-based into v
    vn: np.ndarray | None = None          # (n, 3) normals, aligned with v
    uv: np.ndarray | None = None          # (n, 2) UVs, aligned with v
    markers: list = field(default_factory=list)   # [(pos, line), ...]

    def keep_faces(self, keep: np.ndarray) -> "Mesh":
        """Keep faces where `keep` is True and drop now-unused vertices."""
        faces = self.faces[keep]
        used = np.zeros(len(self.v), dtype=bool)
        used[faces.ravel()] = True
        remap = np.cumsum(used) - 1
        before = np.concatenate([[0], np.cumsum(keep)])
        return Mesh(self.v[used], remap[faces],
                    None if self.vn is None else self.vn[used],
                    None if self.uv is None else self.uv[used],
                    [(int(before[p]), line) for p, line in self.markers])

    def append(self, other: "Mesh", marker: str | None = None) -> "Mesh":
        """This mesh followed by `other` (both must carry the same attributes)."""
        n, f = len(self.v), len(self.faces)
        cat = lambda a, b: None if a is None or b is None else np.vstack([a, b])
        markers = list(self.markers)
        if marker:
            markers.append((f, marker))
        markers += [(p + f, line) for p, line in other.markers]
        return Mesh(np.vstack([self.v, other.v]),
                    np.vstack([self.faces, other.faces + n]),
                    cat(self.vn, other.vn), cat(self.uv, other.uv), markers)

    def face_normals(self) -> np.ndarray:
        a, b, c = (self.v[self.faces[:, i]] for i in range(3))
        return np.cross(b - a, c - a)


# ------------------------------------------------------------------ reading

def _parse_block(lines, tag_len, dtype, slow):
    """Parse same-shaped OBJ lines (e.g. all 'v x y z') in one numpy call.
    Returns an (n, k) array, or None if the lines aren't uniform."""
    if not lines:
        return None
    slow = slow or (lambda b: b)
    k = len(slow(lines[0][tag_len:]).split())
    blob = slow(b" ".join(line[tag_len:] for line in lines))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")       # short reads are caught below
            arr = np.fromstring(blob, dtype=dtype, sep=" ")
    except ValueError:
        return None
    if k <= 0 or arr.size != len(lines) * k:
        return None
    return arr.reshape(len(lines), k)


def _fail(path, kind, lines, parse):
    """Find and report the first bad line (only runs when the fast path fails)."""
    for line in lines:
        try:
            parse(line)
        except (ValueError, IndexError):
            sys.exit(f"{path}: can't parse {kind} line: {line.decode(errors='replace')!r}")
    sys.exit(f"{path}: {kind} lines have inconsistent lengths")


def read_obj(path: str) -> Mesh:
    """Read positions, normals and faces into a Mesh (one vertex per unique
    position/normal index pair). Input UVs and vertex colours are ignored."""
    with open(path, "rb") as f:
        data = f.read()

    v_lines, vn_lines, f_lines, markers = [], [], [], []
    negative = False
    for line in data.splitlines():
        head = line[:3]
        if head[:2] == b"v ":
            v_lines.append(line)
        elif head == b"vn ":
            vn_lines.append(line)
        elif head[:2] == b"f ":
            f_lines.append(line)
            if b"-" in line:
                negative = True
        elif line[:1] in (b"o", b"g", b"s", b"u", b"m"):   # o g s usemtl mtllib
            tag = line.split(None, 1)[0] if line.strip() else b""
            if tag in (b"o", b"g", b"s", b"usemtl", b"mtllib"):
                markers.append((len(f_lines), line.decode(errors="replace").rstrip()))
    del data
    if negative:
        sys.exit(f"{path}: relative (negative) face indices aren't supported")
    if not v_lines:
        sys.exit("No vertices found in the file.")
    if not f_lines:
        sys.exit("No faces found in the file.")

    # Positions: first three numbers of each 'v' line (any colour is dropped).
    V = _parse_block(v_lines, 2, float, None)
    if V is not None and V.shape[1] >= 3:
        V = V[:, :3]
    else:
        V = None
    if V is None:                                  # mixed lengths: slow path
        try:
            V = np.array([[float(x) for x in line.split()[1:4]] for line in v_lines])
        except ValueError:
            V = None
    if V is None or V.ndim != 2 or V.shape[1] != 3:
        _fail(path, "vertex", v_lines, lambda l: [float(x) for x in l.split()[1:4]])
    V = V[:, :3]

    VN = _parse_block(vn_lines, 3, float, None)
    if vn_lines and VN is None:
        _fail(path, "normal", vn_lines, lambda l: [float(x) for x in l.split()[1:4]])
    if VN is not None:
        VN = VN[:, :3]

    # Faces: 'a', 'a/b', 'a//c' or 'a/b/c' per corner; '//' becomes an empty 0.
    tok = f_lines[0].split()[1]
    per = tok.count(b"/") + 1
    F = _parse_block(f_lines, 2, np.int64,
                     lambda b: b.replace(b"//", b"/0/").replace(b"/", b" "))
    if F is None or F.shape[1] % per or F.shape[1] < 3 * per:
        sys.exit(f"{path}: faces must all have the same number of corners and "
                 "the same index format (e.g. all 'a//c')")
    corners = F.shape[1] // per
    F = F.reshape(len(F), corners, per)
    pos = F[:, :, 0] - 1
    nrm = F[:, :, 2] - 1 if per == 3 and VN is not None else None
    if nrm is not None and (nrm < 0).any():
        nrm = None                                 # some corners lack normals

    if corners > 3:                                # fan-triangulate
        fan = [[0, i, i + 1] for i in range(1, corners - 1)]
        pos = np.concatenate([pos[:, t] for t in fan])
        nrm = None if nrm is None else np.concatenate([nrm[:, t] for t in fan])
        order = np.argsort(np.tile(np.arange(len(F)), len(fan)), kind="stable")
        pos = pos[order]
        nrm = None if nrm is None else nrm[order]
        markers = [(p * len(fan), line) for p, line in markers]
        print(f"note: fan-triangulated {corners}-sided faces", file=sys.stderr)

    bad = (pos < 0).any(1) | (pos >= len(V)).any(1)
    if nrm is not None:
        bad |= (nrm >= len(VN)).any(1)
    if bad.any():
        sys.exit(f"{path}: {int(bad.sum())} faces reference vertices that don't exist")

    # One mesh vertex per unique (position index, normal index) pair.
    key = pos if nrm is None else pos * len(VN) + nrm
    uniq, inv = np.unique(key.ravel(), return_inverse=True)
    faces = inv.reshape(-1, 3)
    if nrm is None:
        return Mesh(V[uniq], faces, None, None, markers)
    return Mesh(V[uniq // len(VN)], faces, VN[uniq % len(VN)], None, markers)


def weld(mesh: Mesh, smooth_angle: float = 30.0) -> Mesh:
    """Merge vertices at the same position.

    With normals, corners whose normals are within `smooth_angle` degrees of
    the average at that point merge into one vertex with the averaged normal
    (smooth shading); sharper corners, e.g. along a vertical wall, stay
    separate with their own normal. smooth_angle <= 0 only merges vertices
    whose normals already match exactly (keeps a faceted look).
    """
    q = np.round(mesh.v / WELD_EPS).astype(np.int64)
    _, pid = np.unique(q, axis=0, return_inverse=True)
    pid = pid.ravel()
    n_pos = pid.max() + 1

    if mesh.vn is None:
        key = pid[:, None]
        vn = None
    else:
        nrm = mesh.vn / np.maximum(np.linalg.norm(mesh.vn, axis=1, keepdims=True), 1e-12)
        cluster = np.full(len(nrm), -1)
        vn = nrm.copy()
        if smooth_angle > 0:
            # Peel off up to 6 normal clusters per point: seed with a real
            # corner normal, gather corners within the angle, re-average,
            # repeat. Corners left over keep their own normal.
            cos_a = math.cos(math.radians(smooth_angle))
            for c in range(6):
                rest = np.where(cluster < 0)[0]
                if len(rest) == 0:
                    break
                _, seed = np.unique(pid[rest], return_index=True)
                avg = np.zeros((n_pos, 3))
                avg[pid[rest[seed]]] = nrm[rest[seed]]
                for _ in range(3):
                    member = np.einsum("ij,ij->i", nrm[rest], avg[pid[rest]]) >= cos_a
                    new = np.zeros((n_pos, 3))
                    np.add.at(new, pid[rest[member]], nrm[rest[member]])
                    got = np.linalg.norm(new, axis=1) > 1e-12
                    avg[got] = new[got] / np.linalg.norm(new[got], axis=1, keepdims=True)
                member = np.einsum("ij,ij->i", nrm[rest], avg[pid[rest]]) >= cos_a
                cluster[rest[member]] = c
                vn[rest[member]] = avg[pid[rest[member]]]
        qn = np.round(nrm / NORMAL_EPS).astype(np.int64)
        qn[cluster >= 0] = 0
        key = np.column_stack([pid, cluster, qn])

    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    faces = inv[mesh.faces]
    keep = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & \
           (faces[:, 2] != faces[:, 0])           # drop faces that collapsed
    welded = Mesh(mesh.v[first], faces, None if vn is None else vn[first],
                  None if mesh.uv is None else mesh.uv[first], mesh.markers)
    return welded if keep.all() else welded.keep_faces(keep)


# ------------------------------------------------------------------ writing

def _write_rows(f, fmt, arr, chunk=100_000):
    for s in range(0, len(arr), chunk):
        a = arr[s:s + chunk]
        f.write((fmt * len(a)) % tuple(a.ravel().tolist()))


def write_obj(path: str, mesh: Mesh, comment: str):
    """Write v / vt / vn (one of each per vertex) and faces, with markers."""
    have_uv, have_vn = mesh.uv is not None, mesh.vn is not None
    corner = "%d" + ("/%d" if have_uv else ("/" if have_vn else "")) + ("/%d" if have_vn else "")
    reps = 1 + have_uv + have_vn
    face_fmt = "f " + " ".join([corner] * 3) + "\n"
    with open(path, "w", encoding="utf-8") as f:
        for line in comment.splitlines():
            f.write(f"# {line}\n")
        for _p, line in mesh.markers:
            if line.startswith("mtllib"):
                f.write(line + "\n")
        _write_rows(f, "v %.6f %.6f %.6f\n", mesh.v)
        if have_uv:
            _write_rows(f, "vt %.6f %.6f\n", mesh.uv)
        if have_vn:
            _write_rows(f, "vn %.6f %.6f %.6f\n", mesh.vn)
        idx = np.repeat(mesh.faces + 1, reps, axis=1)
        markers = [(p, l) for p, l in mesh.markers if not l.startswith("mtllib")]
        start = 0
        for pos, line in markers + [(len(idx), None)]:
            _write_rows(f, face_fmt, idx[start:pos])
            start = max(start, pos)
            if line is not None:
                f.write(line + "\n")


# ----------------------------------------------------------------- geometry

def convex_hull(xz: np.ndarray) -> np.ndarray:
    """Convex hull of (x, z) points, counter-clockwise, starting from the
    point with the smallest x (then z). Collinear points are excluded."""
    pts = np.unique(xz, axis=0)
    if len(pts) <= 2:
        return pts
    hull = pts[ConvexHull(pts).vertices]
    start = np.lexsort((hull[:, 1], hull[:, 0]))[0]
    return np.roll(hull, -start, axis=0)


def _cross2(o, a, b):
    return (a[..., 0] - o[..., 0]) * (b[..., 1] - o[..., 1]) - \
           (a[..., 1] - o[..., 1]) * (b[..., 0] - o[..., 0])


def simplify_hull(hull: np.ndarray, tolerance: float) -> np.ndarray:
    """Drop hull vertices that bulge less than `tolerance` off the line between
    their neighbours (smallest deviation first), so near-collinear runs become
    one straight edge."""
    pts = hull.copy()
    while tolerance > 0 and len(pts) > 3:
        prev, nxt = np.roll(pts, 1, axis=0), np.roll(pts, -1, axis=0)
        span = np.linalg.norm(nxt - prev, axis=1)
        d = np.where(span > 0, np.abs(_cross2(prev, nxt, pts)) / np.maximum(span, 1e-300),
                     np.linalg.norm(pts - prev, axis=1))
        i = int(np.argmin(d))
        if d[i] > tolerance:
            break
        pts = np.delete(pts, i, axis=0)
    return pts


@dataclass
class Rect:
    """Minimum-area bounding rectangle. (ux, uz) is the reference edge
    direction ('width' axis), (vx, vz) its perpendicular ('height' axis);
    min/max_u/v are the extents along them."""
    ux: float; uz: float; vx: float; vz: float
    min_u: float; max_u: float; min_v: float; max_v: float

    width = property(lambda s: s.max_u - s.min_u)
    height = property(lambda s: s.max_v - s.min_v)
    area = property(lambda s: s.width * s.height)
    angle = property(lambda s: math.degrees(math.atan2(s.uz, s.ux)) % 90.0)
    center = property(lambda s: ((s.min_u + s.max_u) / 2, (s.min_v + s.max_v) / 2))

    @property
    def corners(self):
        to_world = lambda u, v: (u * self.ux + v * self.vx, u * self.uz + v * self.vz)
        return [to_world(self.min_u, self.min_v), to_world(self.max_u, self.min_v),
                to_world(self.max_u, self.max_v), to_world(self.min_u, self.max_v)]

    def to_frame(self, xz: np.ndarray, center=True) -> np.ndarray:
        """Rotate (x, z) so the rectangle is axis-aligned; optionally center it."""
        out = np.stack([xz[:, 0] * self.ux + xz[:, 1] * self.uz,
                        xz[:, 0] * self.vx + xz[:, 1] * self.vz], axis=1)
        return out - np.array(self.center) if center else out

    def uv(self, frame_xz: np.ndarray, uniform: bool) -> np.ndarray:
        """UVs for points in the rotated, centered frame. U runs 0->1 along the
        longer side, V along the shorter (uniform: both scaled by the longer
        side, V centered on 0.5). Points outside the rectangle land outside
        0..1 on the same mapping. U is flipped for Unity's handedness."""
        long_is_x = self.width >= self.height
        long_len, short_len = max(self.width, self.height), min(self.width, self.height)
        long_c = frame_xz[:, 0] if long_is_x else frame_xz[:, 1]
        short_c = frame_xz[:, 1] if long_is_x else frame_xz[:, 0]
        u_long = long_c / long_len + 0.5 if long_len > 0 else np.zeros(len(frame_xz))
        if uniform:
            u_short = short_c / long_len + 0.5 if long_len > 0 else np.full(len(frame_xz), 0.5)
        else:
            u_short = short_c / short_len + 0.5 if short_len > 0 else np.zeros(len(frame_xz))
        u, v = (u_long, u_short) if long_is_x else (u_short, u_long)
        return np.stack([1 - u, v], axis=1)


def min_bounding_rectangle(hull: np.ndarray) -> Rect:
    """Rotating calipers: the minimum-area rectangle around a convex polygon
    has a side flush with one of its edges, so try every edge direction."""
    if len(hull) < 3:
        raise ValueError("Need at least 3 hull points for a bounding rectangle")
    edges = np.roll(hull, -1, axis=0) - hull
    length = np.linalg.norm(edges, axis=1)
    ok = length > 0
    u = edges[ok] / length[ok, None]
    v = np.stack([-u[:, 1], u[:, 0]], axis=1)
    pu, pv = u @ hull.T, v @ hull.T                     # (edges, points)
    area = (pu.max(1) - pu.min(1)) * (pv.max(1) - pv.min(1))
    i = int(np.argmin(area))
    return Rect(u[i, 0], u[i, 1], v[i, 0], v[i, 1],
                pu[i].min(), pu[i].max(), pv[i].min(), pv[i].max())


def rotate_mesh(mesh: Mesh, rect: Rect, uniform_uv: bool) -> Mesh:
    """Rotate so the bounding rectangle is axis-aligned and centered at the
    origin, and give every vertex a UV."""
    xz = rect.to_frame(mesh.v[:, [0, 2]])
    v = np.column_stack([xz[:, 0], mesh.v[:, 1], xz[:, 1]])
    vn = None
    if mesh.vn is not None:
        nxz = rect.to_frame(mesh.vn[:, [0, 2]], center=False)
        vn = np.column_stack([nxz[:, 0], mesh.vn[:, 1], nxz[:, 1]])
    return Mesh(v, mesh.faces, vn, rect.uv(xz, uniform_uv), mesh.markers)


def hull_mesh(frame_hull: np.ndarray, y: float, rect: Rect, uniform_uv: bool) -> Mesh:
    """Fan-triangulated hull at height y, faces up (+Y). The hull is CCW in
    (x, z), which faces down in Y-up right-handed space, so reverse it."""
    ring = frame_hull[::-1]
    n = len(ring)
    v = np.column_stack([ring[:, 0], np.full(n, y), ring[:, 1]])
    faces = np.array([[0, i, i + 1] for i in range(1, n - 1)])
    return Mesh(v, faces, np.tile([0.0, 1.0, 0.0], (n, 1)), rect.uv(ring, uniform_uv),
                [(0, "g hull")])


# ------------------------------------------------------------- top surface

def top_surface(mesh: Mesh, crease_angle: float = 60.0) -> np.ndarray:
    """Faces of the top surface of a closed (3D-print) model.

    Neighbouring upward-facing faces are linked when their normals differ by
    less than `crease_angle`; the linked patch with the largest footprint is
    the top. The sharp rim where the top meets the side walls stops the
    patch, so walls and bottom drop out with no height threshold. Faces that
    touch the bottom plate (a vertex shared with a downward face) are dropped
    too, in case a stray triangle bridges the rim, and so are needle-thin
    triangles, whose normals are numerical noise.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    _, pid = np.unique(np.round(mesh.v / WELD_EPS).astype(np.int64), axis=0,
                       return_inverse=True)
    pid = pid.ravel()
    F = pid[mesh.faces]
    nf, m = len(F), int(pid.max()) + 1
    fn = mesh.face_normals()
    area = np.linalg.norm(fn, axis=1)
    fn = fn / np.maximum(area[:, None], 1e-300)
    # Needle triangles (three nearly collinear points) have meaningless
    # normals; never let them join the top. Quality = area / longest edge^2.
    T = mesh.v[mesh.faces]
    longest = np.max([np.sum((T[:, i] - T[:, (i + 1) % 3]) ** 2, axis=1)
                      for i in range(3)], axis=0)
    needle = area < 0.02 * longest            # 2*area/edge^2 < 0.02 (~1 deg)
    up = (fn[:, 1] > 0) & ~needle

    E = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    key = np.minimum(E[:, 0], E[:, 1]) * m + np.maximum(E[:, 0], E[:, 1])
    fid = np.tile(np.arange(nf), 3)
    o = np.argsort(key, kind="stable")
    key, fid = key[o], fid[o]
    pair = key[1:] == key[:-1]
    a, b = fid[:-1][pair], fid[1:][pair]
    link = up[a] & up[b] & \
        (np.einsum("ij,ij->i", fn[a], fn[b]) >= math.cos(math.radians(crease_angle)))
    graph = coo_matrix((np.ones(link.sum()), (a[link], b[link])), shape=(nf, nf))
    _, label = connected_components(graph, directed=False)
    footprint = np.bincount(label[up], weights=(area * fn[:, 1])[up], minlength=label.max() + 1)
    keep = label == int(np.argmax(footprint))

    bottom = np.zeros(m, dtype=bool)
    bottom[F[(fn[:, 1] < -0.5) & ~needle].ravel()] = True
    keep &= ~bottom[F].any(axis=1)

    # Put back needles that sit fully inside the top (every edge touches kept
    # surface or another such needle), so the surface has no pinholes. Needles
    # hanging off the rim into the walls fail this and stay out.
    cand = needle & ~bottom[F].any(axis=1)
    while True:
        ok = keep | cand
        covered = np.bincount(a[ok[b]], minlength=nf) + np.bincount(b[ok[a]], minlength=nf)
        new = cand & (covered >= 3)
        if (new == cand).all():
            break
        cand = new
    return keep | cand


# --------------------------------------------------------- border extension

def _circular_smooth(arr, sigma):
    """Gaussian-smooth an (n, ...) array along a closed loop (index space).
    Use on evenly spaced samples, so index space == arc length."""
    if sigma <= 0:
        return arr.copy()
    half = min(int(math.ceil(3 * sigma)), (len(arr) - 1) // 2)
    offsets = np.arange(-half, half + 1)
    w = np.exp(-0.5 * (offsets / sigma) ** 2)
    w /= w.sum()
    out = np.zeros_like(arr, dtype=float)
    for o, wi in zip(offsets, w):
        out += wi * np.roll(arr, -o, axis=0)
    return out


def _resample_loop(pts, spacing):
    """Resample a closed polyline at even arc-length spacing. Returns the
    samples and the arc-length fractions (0..1) of samples and originals."""
    seg = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = cum[-1]
    m = max(8, int(round(total / spacing)))
    ts = np.arange(m) * (total / m)
    closed = np.vstack([pts, pts[:1]])
    out = np.column_stack([np.interp(ts, cum, closed[:, j]) for j in range(pts.shape[1])])
    return out, ts / total, cum[:-1] / total


def outer_border(mesh: Mesh):
    """Find the mesh's largest open border loop, by position (ignoring
    normals). Returns (points, loop, avg_edge, ids, n_loops): `loop` indexes
    `points`; `ids` are the matching vertex indices in `mesh`."""
    _, first, pid = np.unique(np.round(mesh.v / WELD_EPS).astype(np.int64), axis=0,
                              return_index=True, return_inverse=True)
    pid = pid.ravel()
    P = mesh.v[first]
    F = pid[mesh.faces]
    E = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    E = E[E[:, 0] != E[:, 1]]
    m = len(P)
    key = np.minimum(E[:, 0], E[:, 1]) * m + np.maximum(E[:, 0], E[:, 1])
    uniq, idx, counts = np.unique(key, return_index=True, return_counts=True)
    avg_edge = float(np.linalg.norm(P[uniq // m] - P[uniq % m], axis=1).mean())

    border = E[idx[counts == 1]]                  # directed, face winding order
    if len(border) < 3:
        raise ValueError("model has no open border to extend")
    # Walk the border edge by edge. Where it touches itself at a vertex
    # (a pinch), cut the walk there so every loop is a simple polygon.
    outgoing = {}
    for a_, b_ in border.tolist():
        outgoing.setdefault(a_, []).append(b_)
    loops = []
    for start in list(outgoing):
        while outgoing.get(start):
            walk, where, cur = [], {}, start
            while True:
                if cur in where:                   # closed a (sub)loop
                    i = where[cur]
                    sub = walk[i:]
                    if len(sub) >= 3:
                        loops.append(np.array(sub))
                    for w in sub:
                        del where[w]
                    del walk[i:]
                if not outgoing.get(cur):
                    break
                where[cur] = len(walk)
                walk.append(cur)
                cur = outgoing[cur].pop()
    if not loops:
        raise ValueError("couldn't trace a closed border loop")

    def xz_area(loop):
        x, z = P[loop, 0], P[loop, 2]
        return abs(np.sum(x * np.roll(z, -1) - np.roll(x, -1) * z))
    loop = max(loops, key=xz_area)
    return P, loop, avg_edge, first[loop], len(loops)


def build_skirt(P, loop, avg_edge, rect, uniform_uv, distance_mm,
                profile="linear", step=None, smooth=40.0, probe=None):
    """Grow rings outward from the border loop, following the terrain's slope.

    Ring 0 is the border itself: the skirt's first row of triangles uses the
    model's own border vertices, so there is no seam. Rings 1..N are resampled at even spacing along the border and pushed
    outward; a zipper strip joins ring 0 to ring 1. The outline is smoothed by
    a Gaussian `smooth` mm wide, fully at the outer ring and not at the border.
    Returns (skirt Mesh, order, stats). The skirt's first `n0` vertices are
    placeholders for the border; `order` gives, for each, its position in
    `loop` (see attach_skirt).
    """
    step = step or avg_edge
    probe = probe or max(30.0, 6 * step)
    n_rings = max(1, int(math.ceil(distance_mm / step)))
    dr = distance_mm / n_rings

    order = np.arange(len(loop))
    x, z = P[loop, 0], P[loop, 2]
    if np.sum(x * np.roll(z, -1) - np.roll(x, -1) * z) < 0:     # work CCW
        order = order[::-1]
    B0 = P[loop[order]]
    n0 = len(B0)

    B, t_new, t_old = _resample_loop(B0, step)
    n = len(B)
    xz = B[:, [0, 2]]
    samples = lambda mm: mm / step

    # Outward directions from a smoothed border (CCW: outward is to the right).
    sm = _circular_smooth(xz, samples(max(smooth, 4 * step)))
    tan = np.roll(sm, -1, axis=0) - np.roll(sm, 1, axis=0)
    tan /= np.maximum(np.linalg.norm(tan, axis=1, keepdims=True), 1e-12)
    out = np.column_stack([tan[:, 1], -tan[:, 0]])

    # Rings grow from a lightly smoothed copy of the border, pushed out just
    # enough to clear every border vertex nearby, so small spikes and notches
    # in the rim can't poke through the first rings.
    base = _circular_smooth(B, samples(max(2 * step, 10.0)))
    btree = cKDTree(B0[:, [0, 2]])
    clear = np.zeros(n)
    for j, near in enumerate(btree.query_ball_point(base[:, [0, 2]], 4 * step)):
        if near:
            d = B0[near][:, [0, 2]] - base[j, [0, 2]]
            along = d @ out[j]
            lateral = np.abs(d @ tan[j])
            ok = lateral <= 2 * step
            if ok.any():
                clear[j] = max(0.0, along[ok].max())
    clear = np.maximum(clear, _circular_smooth(clear, 2.0))
    B = base.copy()
    B[:, 0] += out[:, 0] * clear
    B[:, 2] += out[:, 1] * clear
    xz = B[:, [0, 2]]

    # Slope dy/dr (r outward): sample the surface inland and fit a line
    # through the border point.
    tree = cKDTree(P[:, [0, 2]])
    num, den = np.zeros(n), 0.0
    for d in np.linspace(probe / 6, probe, 6):
        dist, ii = tree.query(xz - out * d, k=8)
        w = 1.0 / np.maximum(dist, 1e-6)
        y = (P[ii, 1] * w).sum(axis=1) / w.sum(axis=1)
        num += -d * (y - B[:, 1])
        den += d * d
    slope = _circular_smooth(num / den, samples(max(smooth, 4 * step)))
    rising = int((slope > 0).sum())
    slope = np.minimum(slope, 0.0)                # never let the skirt climb

    rings = [B0]
    for k in range(1, n_rings + 1):
        r = dr * k
        h = slope * (r - r * r / (2 * distance_mm)) if profile == "ease" else slope * r
        ring = np.column_stack([B[:, 0] + out[:, 0] * r, B[:, 1] + h, B[:, 2] + out[:, 1] * r])
        if smooth > 0:
            t = k / n_rings
            t = t * t * (3 - 2 * t)               # smoothstep: 0 at border, 1 at edge
            ring = (1 - t) * ring + t * _circular_smooth(ring, samples(smooth))
        rings.append(ring)
    V = np.vstack(rings)                          # ring 0: [0, n0); ring k: n0+(k-1)n

    tris = []
    i = j = 0                                     # zipper by arc-length fraction
    while i < n0 or j < n:
        ti = t_old[i + 1] if i + 1 < n0 else 1.0
        tj = t_new[j + 1] if j + 1 < n else 1.0
        if j >= n or (i < n0 and ti <= tj):
            tris.append((i % n0, (i + 1) % n0, n0 + j % n)); i += 1
        else:
            tris.append((i % n0, n0 + (j + 1) % n, n0 + j % n)); j += 1
    q = np.arange(n)
    for k in range(1, n_rings):
        a, b = n0 + (k - 1) * n + q, n0 + (k - 1) * n + (q + 1) % n
        c, d = n0 + k * n + (q + 1) % n, n0 + k * n + q
        tris += np.column_stack([a, b, c]).tolist() + np.column_stack([a, c, d]).tolist()
    skirt = Mesh(V, np.array(tris))

    fn = skirt.face_normals()
    if (fn[:, 1] < 0).sum() > len(fn) / 2:        # match the model: faces up
        skirt.faces = skirt.faces[:, [0, 2, 1]]
        fn = -fn
    flipped = int((fn[:, 1] <= 0).sum())          # folded / self-crossing

    vn = np.zeros_like(V)
    for c in range(3):
        np.add.at(vn, skirt.faces[:, c], fn)
    skirt.vn = vn / np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-12)
    skirt.uv = rect.uv(V[:, [0, 2]], uniform_uv)

    stats = {"rings": n_rings, "border": n0, "per_ring": n, "step": dr,
             "probe": probe, "rising": rising, "flipped": flipped,
             "slope_min": float(slope.min()), "slope_med": float(np.median(slope)),
             "y_min": float(V[:, 1].min())}
    return skirt, order, stats


def attach_skirt(model: Mesh, skirt: Mesh, ring0_ids: np.ndarray) -> Mesh:
    """Append the skirt to the model, pointing its ring-0 corners at the
    model's own border vertices (shared, not duplicated). The model's border
    normals are blended with the skirt's so shading is continuous."""
    n0, n = len(ring0_ids), len(model.v)
    remap = np.concatenate([ring0_ids, n + np.arange(len(skirt.v) - n0)])
    vn = None
    if model.vn is not None:
        vn = model.vn.copy()
        b = vn[ring0_ids] + skirt.vn[:n0]
        vn[ring0_ids] = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-12)
        vn = np.vstack([vn, skirt.vn[n0:]])
    return Mesh(np.vstack([model.v, skirt.v[n0:]]),
                np.vstack([model.faces, remap[skirt.faces]]), vn,
                None if model.uv is None else np.vstack([model.uv, skirt.uv[n0:]]),
                model.markers + [(len(model.faces), "g extension")])


# ------------------------------------------------------------------- output

def suffixed(path, suffix):
    root, ext = os.path.splitext(path)
    return root + suffix + (ext or ".obj")


def feet_inches(mm):
    feet = mm / MM_PER_FOOT
    whole = int(feet)
    return feet, whole, (feet - whole) * 12


def print_bbox(name, v, p):
    print(f"\n{name} bounding box:")
    for axis, col in zip("XYZ", v.T):
        lo, hi = col.min(), col.max()
        print(f"  {axis}: {lo:.{p}f} to {hi:.{p}f} mm   "
              f"(size {hi - lo:.{p}f} mm / {(hi - lo) / MM_PER_FOOT:.{p}f} ft)")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("obj", help="path to the .obj file (units assumed to be mm)")
    parser.add_argument("-p", "--precision", type=int, default=2,
                        help="decimal places to print (default: 2)")
    parser.add_argument("-t", "--tolerance", type=float, default=5.0,
                        help="in mm: hull vertices deviating less than this from "
                             "the line between their neighbours are treated as "
                             "collinear and dropped (default: 5, use 0 to keep "
                             "the exact hull)")
    parser.add_argument("--hull-y", type=float, default=None,
                        help="Y height for the hull mesh (default: the model's "
                             "lowest Y)")
    parser.add_argument("--no-write", action="store_true",
                        help="print measurements only; don't write OBJ files")
    parser.add_argument("--no-weld", dest="weld", action="store_false",
                        help="don't merge vertices at the same position "
                             "(keeps triangle-soup files as they are)")
    parser.add_argument("--smooth-angle", type=float, default=30.0,
                        help="in degrees: when welding, corners whose normals "
                             "are within this angle share one averaged normal "
                             "(smooth shading); sharper creases stay split. "
                             "0 = only merge identical normals (faceted look) "
                             "(default: 30)")
    parser.add_argument("--stretch-uv", dest="uniform_uv", action="store_false",
                        help="stretch U and V independently to 0->1 on each "
                             "axis (distorts non-square footprints). Default "
                             "is uniform: both axes scaled by the bounding "
                             "box's longer axis, so the mapping isn't "
                             "distorted; the shorter axis is centered on 0.5")
    parser.add_argument("--no-back-method", choices=("crease", "height", "normal"),
                        default="crease",
                        help="how _no_back strips a closed print model down "
                             "to its top surface. crease (default): keep the "
                             "smooth upward-facing patch bounded by the sharp "
                             "rim where top meets sides. height: drop faces "
                             "with any vertex below --no-back-threshold. "
                             "normal: drop downward-facing faces")
    parser.add_argument("--crease-angle", type=float, default=60.0,
                        help="in degrees: with --no-back-method crease, faces "
                             "meeting at a sharper angle than this are split "
                             "apart (default: 60)")
    parser.add_argument("--no-back-by-normal", dest="no_back_method",
                        action="store_const", const="normal",
                        help="same as --no-back-method normal")
    parser.add_argument("--no-back-threshold", type=float, default=50.0,
                        help="in mm: with --no-back-method height, vertices "
                             "below this Y value are removed (default: 50)")
    parser.add_argument("--extend-distance", type=float, default=100.0,
                        help="in mm: how far the _extended skirt reaches "
                             "outward from the _no_back border (default: 100)")
    parser.add_argument("--extend-profile", choices=("linear", "ease"), default="linear",
                        help="linear: keep the border's slope all the way out. "
                             "ease: start at the border's slope and flatten "
                             "gradually to level at the outer edge "
                             "(default: linear)")
    parser.add_argument("--extend-step", type=float, default=None,
                        help="in mm: spacing between skirt rings (default: "
                             "the mesh's average edge length)")
    parser.add_argument("--extend-smooth", type=float, default=40.0,
                        help="in mm: how much the skirt's outline is "
                             "smoothed (Gaussian width along the border). "
                             "Full at the outer ring, none at the border "
                             "(default: 40, 0 = follow the border's shape)")
    parser.add_argument("--extend-probe", type=float, default=None,
                        help="in mm: how far inland to sample when measuring "
                             "the slope at the border (default: max(30, 6 x "
                             "ring spacing))")
    parser.add_argument("--no-extend", action="store_true",
                        help="don't write the _extended model")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    p = args.precision
    src = os.path.basename(args.obj)

    mesh = read_obj(args.obj)
    n_read_v, n_faces = len(mesh.v), len(mesh.faces)
    if args.weld:
        mesh = weld(mesh, args.smooth_angle)

    raw_hull = convex_hull(mesh.v[:, [0, 2]])
    if len(raw_hull) < 3:
        sys.exit("Convex hull is degenerate (all vertices lie on a line or point).")
    hull = simplify_hull(raw_hull, args.tolerance)

    print(f"{n_read_v} vertices read, {n_faces} faces, "
          f"{len(raw_hull)} vertices on the convex hull (XZ plane)")
    if args.weld:
        print(f"welded to {len(mesh.v)} vertices"
              + (f" (smoothing angle {args.smooth_angle:g} deg)"
                 if mesh.vn is not None and args.smooth_angle > 0 else ""))
    if len(hull) != len(raw_hull):
        print(f"simplified to {len(hull)} (tolerance {args.tolerance:g} mm)")
    print()

    print("Hull vertices (counter-clockwise, X / Z in mm):")
    for i, (x, z) in enumerate(hull):
        print(f"  {i:>3}: ({x:.{p}f}, {z:.{p}f})")

    print("\nEdge lengths between successive hull vertices:")
    lengths = np.linalg.norm(np.roll(hull, -1, axis=0) - hull, axis=1)
    for i, d in enumerate(lengths):
        feet, whole, inches = feet_inches(d)
        print(f"  {i:>3} -> {(i + 1) % len(hull):<3}  {d:>11.{p}f} mm   "
              f"{feet:>8.{p}f} ft   ({whole}' {inches:.1f}\")")
    total = lengths.sum()
    print(f"\nTotal perimeter: {total:.{p}f} mm   {total / MM_PER_FOOT:.{p}f} ft")

    rect = min_bounding_rectangle(hull)
    long_side, short_side = max(rect.width, rect.height), min(rect.width, rect.height)
    aspect = long_side / short_side if short_side > 0 else float("inf")
    print("\nMinimum bounding rectangle:")
    print(f"  size:     {rect.width:.{p}f} x {rect.height:.{p}f} mm   "
          f"({rect.width / MM_PER_FOOT:.{p}f} x {rect.height / MM_PER_FOOT:.{p}f} ft)")
    print(f"  rotation: {rect.angle:.{p}f} deg from +X axis")
    print(f"  area:     {rect.area:.{p}f} mm^2   ({rect.area / MM_PER_FOOT ** 2:.{p}f} ft^2)")
    print(f"  aspect ratio (long:short): {aspect:.{p}f} : 1")
    print("  corners (X, Z in mm):")
    for i, (x, z) in enumerate(rect.corners):
        print(f"    {i}: ({x:.{p}f}, {z:.{p}f})")

    if args.no_write:
        return

    # Everything written below is in the rotated, axis-aligned, centered frame.
    rot = rotate_mesh(mesh, rect, args.uniform_uv)
    cx, cz = rect.center
    path = suffixed(args.obj, "_rotated_uvs")
    write_obj(path, rot, f"{src} rotated to axis-align its minimum bounding "
                         "rectangle, with UVs added")
    long_axis = "X" if rect.width >= rect.height else "Z"
    uv_note = "uniform (no stretch)" if args.uniform_uv else "stretched to 0->1 on both axes"
    print(f"\nWrote {path}: {len(rot.v)} vertices rotated {rect.angle:.{p}f} deg and "
          f"centered (bounding box was at X={cx:.{p}f}, Z={cz:.{p}f}, now "
          f"axis-aligned and centered at origin), UVs along rotated-{long_axis} "
          f"as U ({uv_note})")

    hull_y = args.hull_y if args.hull_y is not None else float(rot.v[:, 1].min())
    hm = hull_mesh(rect.to_frame(hull), hull_y, rect, args.uniform_uv)
    path = suffixed(args.obj, "_hull")
    write_obj(path, hm, f"convex hull (XZ) of {src}\n"
                        f"units: mm, {len(hm.v)} boundary vertices, all faces +Y")
    print(f"Wrote {path}: {len(hm.v)} vertices, {len(hm.faces)} triangles "
          f"at Y = {hull_y:.{p}f}")

    if args.no_back_method == "crease":
        keep = top_surface(rot, args.crease_angle)
        comment = "top surface only (sides and back removed)"
        method = f"top surface, split at creases sharper than {args.crease_angle:g} deg"
    elif args.no_back_method == "height":
        keep = (rot.v[rot.faces, 1] >= args.no_back_threshold).all(axis=1)
        comment = f"with vertices below Y={args.no_back_threshold:g} mm removed"
        method = f"by-height (threshold Y = {args.no_back_threshold:.{p}f} mm)"
    else:
        keep = rot.face_normals()[:, 1] >= 0.0
        comment = "with downward-facing polygons removed"
        method = "by-normal (downward-facing polygons removed)"
    nb = rot.keep_faces(keep)
    path = suffixed(args.obj, "_no_back")
    write_obj(path, nb, f"{src} {comment}")
    print(f"Wrote {path}: kept {len(nb.faces)} of {len(rot.faces)} faces, "
          f"{len(nb.v)} of {len(rot.v)} vertices (method: {method})")
    if len(nb.v) == 0:
        print(f"\n{os.path.basename(path)} has no surviving vertices; "
              "no bounding box to report.")
        return
    print_bbox(os.path.basename(path), nb.v, p)

    if args.no_extend:
        return
    P, loop, avg_edge, ids, n_loops = outer_border(nb)
    skirt, order, st = build_skirt(P, loop, avg_edge, rect, args.uniform_uv,
                                   args.extend_distance, profile=args.extend_profile,
                                   step=args.extend_step, smooth=args.extend_smooth,
                                   probe=args.extend_probe)
    ext = attach_skirt(nb, skirt, ids[order])
    path = suffixed(args.obj, "_extended")
    write_obj(path, ext, f"{src} {comment}, extended {args.extend_distance:g} mm "
                         f"outward ({args.extend_profile} slope)")
    u0, w0 = skirt.uv.min(axis=0)
    u1, w1 = skirt.uv.max(axis=0)
    n_new = len(skirt.v) - st["border"]
    print(f"\nWrote {path}: _no_back plus a {args.extend_distance:.{p}f} mm "
          f"skirt ({args.extend_profile} slope)")
    print(f"  border: {st['border']} vertices, resampled to {st['per_ring']} per ring"
          + (f" (largest of {n_loops} border loops)" if n_loops > 1 else ""))
    print(f"  rings:  {st['rings']} at {st['step']:.{p}f} mm spacing (mesh average "
          f"edge {avg_edge:.{p}f} mm), {n_new} new vertices (border shared "
          f"with the model), "
          f"{len(skirt.faces)} new triangles")
    print(f"  slope:  median {st['slope_med']:.3f}, steepest {st['slope_min']:.3f} "
          f"(Y per mm outward, sampled {st['probe']:.{p}f} mm inland)")
    print(f"  lowest skirt Y: {st['y_min']:.{p}f} mm;  "
          f"skirt UVs: U {u0:.3f}..{u1:.3f}, V {w0:.3f}..{w1:.3f}")
    if st["rising"]:
        print(f"  note: {st['rising']} border vertices sloped upward outward; "
              "held level there")
    if st["flipped"] and st["flipped"] <= 1e-3 * len(skirt.faces):
        print(f"  note: {st['flipped']} tiny skirt triangles fold over where the "
              "rim has small notches")
    elif st["flipped"]:
        print(f"  WARNING: {st['flipped']} skirt triangles fold over (rings "
              "crossing at a tight concave spot). Try a larger --extend-smooth "
              "or a smaller --extend-distance.")


if __name__ == "__main__":
    main()
