#!/usr/bin/env python3
"""
Read an OBJ file, project its vertices onto the XZ plane (ignoring Y), compute
the 2D convex hull, print the distance between each pair of successive hull
vertices, and find the minimum-area bounding rectangle around the hull.

Input units are assumed to be millimeters; distances are printed in both
millimeters and feet.

Two OBJ files are also written alongside the input:
  <name>_hull.obj     triangulated convex hull, all face normals pointing +Y
  <name>_no_back.obj  the original model with all downward-facing (-Y) polygons
                      and their now-unused vertices removed

Usage:
    python obj_hull_xz.py model.obj
    python obj_hull_xz.py model.obj --tolerance 5 --precision 1
"""

import argparse
import math
import os
import sys

MM_PER_FOOT = 304.8


# ---------------------------------------------------------------- OBJ parsing

def parse_face(tokens):
    """Parse face tokens into a list of (v, vt, vn) 1-based indices; None where absent."""
    verts = []
    for tok in tokens:
        parts = (tok.split("/") + ["", ""])[:3]
        idx = []
        for part in parts:
            idx.append(int(part) if part.strip() else None)
        verts.append(tuple(idx))
    return verts


def resolve(i, count):
    """Turn a 1-based (possibly negative) OBJ index into a 0-based index."""
    if i is None:
        return None
    return i - 1 if i > 0 else count + i


def read_obj(path):
    """Return (lines, v, vt, vn, faces).

    `lines` is the raw file. `faces` is a list of
    (line_number, [(vi, vti, vni), ...]) with 0-based resolved indices.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    v, vt, vn, faces = [], [], [], []
    for lineno, line in enumerate(lines):
        parts = line.split()
        if not parts:
            continue
        tag = parts[0]
        if tag == "v" and len(parts) >= 4:
            try:
                v.append((float(parts[1]), float(parts[2]), float(parts[3])))
            except ValueError:
                print(f"Skipping unparsable vertex on line {lineno + 1}", file=sys.stderr)
                v.append((0.0, 0.0, 0.0))
        elif tag == "vt":
            vt.append(line)
        elif tag == "vn":
            vn.append(line)
        elif tag == "f" and len(parts) >= 4:
            try:
                raw = parse_face(parts[1:])
            except ValueError:
                print(f"Skipping malformed face on line {lineno + 1}", file=sys.stderr)
                continue
            faces.append((lineno, [(resolve(a, len(v)),
                                    resolve(b, len(vt)),
                                    resolve(c, len(vn))) for a, b, c in raw]))
    return lines, v, vt, vn, faces


# ------------------------------------------------------------------ geometry

def cross2(o, a, b):
    """Z-component of the 2D cross product of OA and OB. >0 means counter-clockwise."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull(points):
    """Andrew's monotone chain on (x, z). Returns hull points counter-clockwise."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross2(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross2(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return lower[:-1] + upper[:-1]


def distance(a, b):
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


def perp_distance(p, a, b):
    """Perpendicular distance from p to the infinite line through a and b."""
    ab = distance(a, b)
    if ab == 0.0:
        return distance(p, a)
    return abs(cross2(a, b, p)) / ab


def simplify_hull(hull, tolerance):
    """Drop hull vertices that bulge less than `tolerance` off the line between
    their neighbours, i.e. treat near-collinear runs as a single straight edge.

    Removal is iterative, smallest deviation first, so a long gentle arc
    collapses cleanly instead of depending on list order.
    """
    if tolerance <= 0:
        return hull

    pts = list(hull)
    while len(pts) > 3:
        n = len(pts)
        worst_i, worst_d = None, None
        for i in range(n):
            d = perp_distance(pts[i], pts[(i - 1) % n], pts[(i + 1) % n])
            if worst_d is None or d < worst_d:
                worst_i, worst_d = i, d
        if worst_d is None or worst_d > tolerance:
            break
        pts.pop(worst_i)
    return pts


def min_bounding_rectangle(hull):
    """Minimum-area bounding rectangle of a convex polygon, via rotating calipers.

    For a convex hull, the minimum-area enclosing rectangle always has one
    side flush with one of the hull's edges, so it's enough to test each edge
    direction, project every hull point onto that direction and its
    perpendicular, and keep the smallest-area result.

    Returns a dict with:
      corners  four (x, z) corners, counter-clockwise
      width    extent along the rectangle's long-axis reference edge (mm)
      height   extent along the perpendicular axis (mm)
      angle    rotation of that reference edge from +X, in degrees, in [0, 90)
      area     width * height (mm^2)
      ux, uz   unit vector of the reference edge (the "width" axis)
      vx, vz   unit vector perpendicular to it (the "height" axis)
      min_u, max_u, min_v, max_v   extents along those two axes (mm)
    """
    n = len(hull)
    if n < 3:
        raise ValueError("Need at least 3 hull points for a bounding rectangle")

    best = None
    for i in range(n):
        x1, z1 = hull[i]
        x2, z2 = hull[(i + 1) % n]
        edge_len = distance((x1, z1), (x2, z2))
        if edge_len == 0.0:
            continue
        ux, uz = (x2 - x1) / edge_len, (z2 - z1) / edge_len  # edge direction
        vx, vz = -uz, ux                                     # perpendicular

        min_u = min_v = float("inf")
        max_u = max_v = float("-inf")
        for x, z in hull:
            pu = x * ux + z * uz
            pv = x * vx + z * vz
            min_u, max_u = min(min_u, pu), max(max_u, pu)
            min_v, max_v = min(min_v, pv), max(max_v, pv)

        width, height = max_u - min_u, max_v - min_v
        area = width * height
        if best is None or area < best[0]:
            best = (area, width, height, ux, uz, vx, vz, min_u, max_u, min_v, max_v)

    area, width, height, ux, uz, vx, vz, min_u, max_u, min_v, max_v = best

    def to_world(u, v):
        return (u * ux + v * vx, u * uz + v * vz)

    corners = [to_world(min_u, min_v), to_world(max_u, min_v),
               to_world(max_u, max_v), to_world(min_u, max_v)]

    angle = math.degrees(math.atan2(uz, ux)) % 90.0

    return {"corners": corners, "width": width, "height": height,
            "angle": angle, "area": area,
            "ux": ux, "uz": uz, "vx": vx, "vz": vz,
            "min_u": min_u, "max_u": max_u, "min_v": min_v, "max_v": max_v}


def newell_normal(points):
    """Newell's method: area-weighted normal of a 3D polygon. Robust for
    non-planar and concave polygons."""
    nx = ny = nz = 0.0
    n = len(points)
    for i in range(n):
        ax, ay, az = points[i]
        bx, by, bz = points[(i + 1) % n]
        nx += (ay - by) * (az + bz)
        ny += (az - bz) * (ax + bx)
        nz += (ax - bx) * (ay + by)
    return nx, ny, nz


def rotate_xz(x, z, rect):
    """Rotate a point into the bounding rectangle's frame: its 'u' basis
    vector becomes +X and its 'v' basis vector becomes +Z, so the rectangle
    is axis-aligned afterward."""
    return (x * rect["ux"] + z * rect["uz"],
            x * rect["vx"] + z * rect["vz"])


def parse_vn_line(line):
    parts = line.split()
    return tuple(float(p) for p in parts[1:4])

def rotate_and_uv_point(x, z, rect, uniform_uv=False):
    """Rotate a single (x, z) point into the bounding rectangle's axis-aligned,
    centered frame, and compute its UV.

    By default U runs 0->1 along the rectangle's longer axis and V runs 0->1
    along its shorter axis (stretches whichever axis is shorter). If
    uniform_uv is True, both axes are scaled by the *longer* axis's length
    instead, so mm-per-UV-unit matches on both axes; V is then centered on
    0.5 since it no longer reaches a full 0->1 span.

    Returns (new_x, new_z, u, v).
    """
    width, height = rect["width"], rect["height"]
    min_u, max_u = rect["min_u"], rect["max_u"]
    min_v, max_v = rect["min_v"], rect["max_v"]
    long_axis_is_u = width >= height  # True: rotated-X is the longer side
    long_len = max(width, height)
    center_x = (min_u + max_u) / 2.0
    center_z = (min_v + max_v) / 2.0

    rx, rz = rotate_xz(x, z, rect)
    if long_axis_is_u:
        long_coord, short_coord = rx, rz
        min_long, min_short = min_u, min_v
        short_center = center_z
        short_len = height
    else:
        long_coord, short_coord = rz, rx
        min_long, min_short = min_v, min_u
        short_center = center_x
        short_len = width

    u_long = (long_coord - min_long) / long_len if long_len > 0 else 0.0
    if uniform_uv:
        u_short = 0.5 + (short_coord - short_center) / long_len if long_len > 0 else 0.5
    else:
        u_short = (short_coord - min_short) / short_len if short_len > 0 else 0.0

    u, v = (u_long, u_short) if long_axis_is_u else (u_short, u_long)

    # Flip U because OBJ uses a different coordinate system than Unity and unity flips an axis

    return rx - center_x, rz - center_z, 1 - u, v


def build_rotated_uv_model(v, vn_lines, rect, uniform_uv=False):
    """Rotate every vertex (and vertex normal) so the bounding rectangle is
    axis-aligned, translate so the rectangle's center sits at X=0, Z=0, and
    compute a UV per vertex (see rotate_and_uv_point for the UV convention).

    Returns (new_v, uvs, new_vn, center_x, center_z).
    """
    min_u, max_u = rect["min_u"], rect["max_u"]
    min_v, max_v = rect["min_v"], rect["max_v"]
    center_x = (min_u + max_u) / 2.0
    center_z = (min_v + max_v) / 2.0

    new_v, uvs = [], []
    for x, y, z in v:
        nx, nz, u, uvv = rotate_and_uv_point(x, z, rect, uniform_uv)
        new_v.append((nx, y, nz))
        uvs.append((u, uvv))

    new_vn = []
    for line in vn_lines:
        nx, ny, nz = parse_vn_line(line)
        rnx, rnz = rotate_xz(nx, nz, rect)
        new_vn.append((rnx, ny, rnz))

    return new_v, uvs, new_vn, center_x, center_z


def write_hull_obj(path, hull, hull_uvs, y, source):
    """Write the hull as a fan-triangulated mesh with all normals pointing +Y
    and a UV per boundary vertex."""
    # Hull points are counter-clockwise in (x, z), which winds *downward* in a
    # Y-up right-handed system, so reverse the order to face +Y. Reverse the
    # matching UVs the same way to keep them paired with the right vertex.
    ring = list(reversed(hull))
    ring_uvs = list(reversed(hull_uvs))

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# convex hull (XZ) of {os.path.basename(source)}\n")
        f.write(f"# units: mm, {len(ring)} boundary vertices, all faces +Y\n")
        f.write("vn 0.0 1.0 0.0\n")
        for x, z in ring:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for u, v in ring_uvs:
            f.write(f"vt {u:.6f} {v:.6f}\n")
        f.write("g hull\n")
        for i in range(1, len(ring) - 1):
            f.write(f"f 1/1/1 {i + 1}/{i + 1}/1 {i + 2}/{i + 2}/1\n")
    return len(ring) - 2


def select_faces_by_normal(v, faces):
    """Keep faces whose (Newell) normal does not point downward (-Y)."""
    kept = []
    for lineno, verts in faces:
        pts = [v[a] for a, _b, _c in verts if a is not None and 0 <= a < len(v)]
        if len(pts) < 3:
            continue
        if newell_normal(pts)[1] < 0.0:
            continue
        kept.append((lineno, verts))
    return kept


def select_faces_by_height(v, faces, threshold):
    """Keep faces none of whose vertices fall below the given Y threshold."""
    kept = []
    for lineno, verts in faces:
        idxs = [a for a, _b, _c in verts if a is not None]
        if not idxs:
            continue
        if any(not (0 <= a < len(v)) or v[a][1] < threshold for a in idxs):
            continue
        kept.append((lineno, verts))
    return kept


def write_filtered_obj(path, lines, v, uvs, vn, faces, source, kept_faces, comment):
    """Write only `kept_faces` and the vertices/normals/UVs they still use,
    replaying the original file's structural lines (groups, materials).
    `uvs` is a per-vertex UV list aligned 1:1 with `v`.

    Returns (kept_face_count, total_face_count, kept_vertex_coords), where
    kept_vertex_coords lists the (x, y, z) actually written, in file order.
    """
    used_v, used_vn = set(), set()
    for _lineno, verts in kept_faces:
        for a, _b, c in verts:
            if a is not None:
                used_v.add(a)
            if c is not None:
                used_vn.add(c)

    def remap(used):
        return {old: new for new, old in enumerate(sorted(used))}

    map_v = remap(used_v)
    map_vn = remap(used_vn)
    have_vn = len(used_vn) > 0

    face_by_line = dict(kept_faces)
    face_lines = {lineno for lineno, _ in faces}

    kept_vertex_coords = [v[old] for old in sorted(used_v)]

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# {os.path.basename(source)} {comment}\n")
        for x, y, z in kept_vertex_coords:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for old in sorted(used_v):
            u, uvv = uvs[old]
            f.write(f"vt {u:.6f} {uvv:.6f}\n")
        for old in sorted(used_vn):
            f.write(vn[old] + "\n")

        # Replay the original file for structure (groups, materials), emitting
        # only the faces we kept, with re-mapped indices.
        for lineno, line in enumerate(lines):
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag in ("v", "vt", "vn", "#"):
                continue
            if tag == "f":
                if lineno not in face_lines:
                    continue
                verts = face_by_line.get(lineno)
                if verts is None:
                    continue
                out = []
                for a, _b, c in verts:
                    if a is None:
                        continue
                    sa = str(map_v[a] + 1)
                    sb = sa  # UV is 1:1 with vertex index
                    if have_vn and c is not None:
                        sc = str(map_vn[c] + 1)
                        out.append(f"{sa}/{sb}/{sc}")
                    else:
                        out.append(f"{sa}/{sb}")
                f.write("f " + " ".join(out) + "\n")
            else:
                f.write(line + "\n")

    return len(kept_faces), len(faces), kept_vertex_coords


def write_rotated_uv_obj(path, lines, new_v, uvs, new_vn, faces, source):
    """Write the full original model with vertices (and normals) rotated so
    the bounding rectangle is axis-aligned, and a fresh UV per vertex."""
    face_by_line = dict(faces)

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# {os.path.basename(source)} rotated to axis-align its "
                "minimum bounding rectangle, with UVs added\n")
        for x, y, z in new_v:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for u, vv in uvs:
            f.write(f"vt {u:.6f} {vv:.6f}\n")
        for nx, ny, nz in new_vn:
            f.write(f"vn {nx:.6f} {ny:.6f} {nz:.6f}\n")

        have_vn = len(new_vn) > 0
        for lineno, line in enumerate(lines):
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag in ("v", "vt", "vn", "#"):
                continue
            if tag == "f":
                verts = face_by_line.get(lineno)
                if verts is None:
                    continue
                out = []
                for a, _b, c in verts:
                    if a is None:
                        continue
                    sa = str(a + 1)
                    sb = str(a + 1)  # UV is 1:1 with vertex index
                    if have_vn and c is not None:
                        sc = str(c + 1)
                        out.append(f"{sa}/{sb}/{sc}")
                    else:
                        out.append(f"{sa}/{sb}")
                f.write("f " + " ".join(out) + "\n")
            else:
                f.write(line + "\n")


# ------------------------------------------------------------------- driver

def suffixed(path, suffix):
    root, ext = os.path.splitext(path)
    return root + suffix + (ext or ".obj")


def main():
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
    parser.add_argument("--stretch-uv", dest="uniform_uv", action="store_false",
                        help="stretch U and V independently to 0->1 on each "
                             "axis (distorts non-square footprints). Default "
                             "is uniform: both axes scaled by the bounding "
                             "box's longer axis, so the mapping isn't "
                             "distorted; the shorter axis is centered on 0.5")
    parser.add_argument("--no-back-by-normal", dest="no_back_by_height",
                        action="store_false",
                        help="build the _no_back model by removing "
                             "downward-facing (-Y normal) polygons instead of "
                             "the default: removing vertices (and faces that "
                             "use them) below --no-back-threshold")
    parser.add_argument("--no-back-threshold", type=float, default=50.0,
                        help="in mm: with the default by-height method, "
                             "vertices below this Y value are removed "
                             "(default: 50)")
    parser.set_defaults(uniform_uv=True, no_back_by_height=True)
    args = parser.parse_args()
    p = args.precision

    lines, v, vt, vn, faces = read_obj(args.obj)
    if not v:
        sys.exit("No vertices found in the file.")

    hull = convex_hull([(x, z) for x, _y, z in v])
    if len(hull) < 2:
        sys.exit("Convex hull is degenerate (all vertices project to one point).")

    raw_count = len(hull)
    hull = simplify_hull(hull, args.tolerance)

    print(f"{len(v)} vertices read, {len(faces)} faces, "
          f"{raw_count} vertices on the convex hull (XZ plane)")
    if len(hull) != raw_count:
        print(f"simplified to {len(hull)} (tolerance {args.tolerance:g} mm)")
    print()

    print("Hull vertices (counter-clockwise, X / Z in mm):")
    for i, (x, z) in enumerate(hull):
        print(f"  {i:>3}: ({x:.{p}f}, {z:.{p}f})")

    print("\nEdge lengths between successive hull vertices:")
    total = 0.0
    n = len(hull)
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        d = distance(a, b)
        total += d
        feet = d / MM_PER_FOOT
        whole_ft = int(feet)
        inches = (feet - whole_ft) * 12
        print(f"  {i:>3} -> {(i + 1) % n:<3}  {d:>11.{p}f} mm   "
              f"{feet:>8.{p}f} ft   ({whole_ft}' {inches:.1f}\")")

    print(f"\nTotal perimeter: {total:.{p}f} mm   {total / MM_PER_FOOT:.{p}f} ft")

    rect = min_bounding_rectangle(hull)
    w_ft = rect["width"] / MM_PER_FOOT
    h_ft = rect["height"] / MM_PER_FOOT
    area_ft2 = rect["area"] / (MM_PER_FOOT ** 2)
    long_side, short_side = max(rect["width"], rect["height"]), min(rect["width"], rect["height"])
    aspect_ratio = (long_side / short_side) if short_side > 0 else float("inf")
    print("\nMinimum bounding rectangle:")
    print(f"  size:     {rect['width']:.{p}f} x {rect['height']:.{p}f} mm   "
          f"({w_ft:.{p}f} x {h_ft:.{p}f} ft)")
    print(f"  rotation: {rect['angle']:.{p}f} deg from +X axis")
    print(f"  area:     {rect['area']:.{p}f} mm^2   ({area_ft2:.{p}f} ft^2)")
    print(f"  aspect ratio (long:short): {aspect_ratio:.{p}f} : 1")
    print("  corners (X, Z in mm):")
    for i, (x, z) in enumerate(rect["corners"]):
        print(f"    {i}: ({x:.{p}f}, {z:.{p}f})")

    if args.no_write:
        return

    # Rotate and center the model first; the hull and no-back files are then
    # built from this rotated, axis-aligned, centered geometry.
    new_v, uvs, new_vn, center_x, center_z = build_rotated_uv_model(
        v, vn, rect, uniform_uv=args.uniform_uv)
    ru_path = suffixed(args.obj, "_rotated_uvs")
    write_rotated_uv_obj(ru_path, lines, new_v, uvs, new_vn, faces, args.obj)
    long_axis = "X" if rect["width"] >= rect["height"] else "Z"
    uv_note = "uniform (no stretch)" if args.uniform_uv else "stretched to 0->1 on both axes"
    print(f"\nWrote {ru_path}: {len(new_v)} vertices rotated {rect['angle']:.{p}f} "
          f"deg and centered (bounding box was at X={center_x:.{p}f}, "
          f"Z={center_z:.{p}f}, now axis-aligned and centered at origin), "
          f"UVs along rotated-{long_axis} as U ({uv_note})")

    hull_data = [rotate_and_uv_point(x, z, rect, args.uniform_uv) for x, z in hull]
    rotated_hull = [(nx, nz) for nx, nz, _u, _v in hull_data]
    hull_uvs = [(u, v) for _nx, _nz, u, v in hull_data]

    hull_y = args.hull_y if args.hull_y is not None else min(y for _x, y, _z in new_v)
    hull_path = suffixed(args.obj, "_hull")
    tris = write_hull_obj(hull_path, rotated_hull, hull_uvs, hull_y, args.obj)
    print(f"Wrote {hull_path}: {len(rotated_hull)} vertices, {tris} triangles "
          f"at Y = {hull_y:.{p}f} (from the rotated, centered model)")

    rotated_vn_lines = [f"vn {nx:.6f} {ny:.6f} {nz:.6f}" for nx, ny, nz in new_vn]
    nb_path = suffixed(args.obj, "_no_back")

    if args.no_back_by_height:
        kept_faces = select_faces_by_height(new_v, faces, args.no_back_threshold)
        comment = (f"with vertices below Y={args.no_back_threshold:g} mm "
                   "removed")
        method_note = f"by-height (threshold Y = {args.no_back_threshold:.{p}f} mm)"
    else:
        kept_faces = select_faces_by_normal(new_v, faces)
        comment = "with downward-facing polygons removed"
        method_note = "by-normal (downward-facing polygons removed)"

    kept_f, all_f, kept_vertex_coords = write_filtered_obj(
        nb_path, lines, new_v, uvs, rotated_vn_lines, faces, args.obj,
        kept_faces, comment)
    kept_v, all_v = len(kept_vertex_coords), len(new_v)
    print(f"Wrote {nb_path}: kept {kept_f} of {all_f} faces, "
          f"{kept_v} of {all_v} vertices (from the rotated, centered model, "
          f"method: {method_note})")

    if kept_vertex_coords:
        xs = [pt[0] for pt in kept_vertex_coords]
        ys = [pt[1] for pt in kept_vertex_coords]
        zs = [pt[2] for pt in kept_vertex_coords]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        min_z, max_z = min(zs), max(zs)
        size_x, size_y, size_z = max_x - min_x, max_y - min_y, max_z - min_z
        print(f"\n{os.path.basename(nb_path)} bounding box:")
        print(f"  X: {min_x:.{p}f} to {max_x:.{p}f} mm   "
              f"(size {size_x:.{p}f} mm / {size_x / MM_PER_FOOT:.{p}f} ft)")
        print(f"  Y: {min_y:.{p}f} to {max_y:.{p}f} mm   "
              f"(size {size_y:.{p}f} mm / {size_y / MM_PER_FOOT:.{p}f} ft)")
        print(f"  Z: {min_z:.{p}f} to {max_z:.{p}f} mm   "
              f"(size {size_z:.{p}f} mm / {size_z / MM_PER_FOOT:.{p}f} ft)")
    else:
        print(f"\n{os.path.basename(nb_path)} has no surviving vertices; "
              "no bounding box to report.")


if __name__ == "__main__":
    main()
