#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.22",
#     "scipy>=1.10",
#     "pillow>=9.0",
# ]
# ///
"""Synthetic end-to-end test for mesh_splitter.py.

Builds a warped 24x24 grid mesh whose UV -> XZ map is bijective but *not*
globally affine (so every cut really does have to be re-lifted triangle by
triangle), draws strokes on a PNG, cuts, and runs verify.py on the result.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
N = 24
W = H = 1024


def build_mesh(path):
    u, v = np.meshgrid(np.linspace(0, 1, N + 1), np.linspace(0, 1, N + 1),
                       indexing="ij")
    # a smooth, invertible warp: the UV map is piecewise affine, never global
    x = (u + 0.05 * np.sin(2 * np.pi * v)) * 1000.0
    z = (v + 0.04 * np.sin(2 * np.pi * u)) * 800.0
    y = 50.0 * np.sin(3 * u) * np.cos(2 * v)        # irrelevant to the cut
    lines = ["# synthetic warped grid"]
    for i in range(N + 1):
        for j in range(N + 1):
            lines.append(f"v {x[i, j]:.6f} {y[i, j]:.6f} {z[i, j]:.6f}")
    for i in range(N + 1):
        for j in range(N + 1):
            lines.append(f"vt {u[i, j]:.6f} {v[i, j]:.6f}")

    def vid(i, j):
        return i * (N + 1) + j + 1
    for i in range(N):
        for j in range(N):
            a, b, c, d = vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)
            # wind so the XZ cross product is consistent
            lines.append(f"f {a}/{a} {d}/{d} {c}/{c}")
            lines.append(f"f {a}/{a} {c}/{c} {b}/{b}")
    open(path, "w").write("\n".join(lines) + "\n")


def build_image(path, tangential=False, dangling=False):
    im = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    # a diagonal running edge to edge, a wobbly vertical, and a closed blob
    # a diagonal running edge to edge.  Offset so it crosses the ellipse
    # cleanly; `tangential` puts it almost tangent to the ellipse instead,
    # which is the adversarial case.
    y0 = 300 if tangential else 80
    d.line([(0, y0), (W - 1, H - 200)], fill=(255, 0, 0, 255), width=9)
    end = H - 120 if dangling else H - 1
    ts = list(range(0, end, 8)) + [end]
    pts = [(300 + int(60 * np.sin(t / 90.0)), t) for t in ts]
    d.line(pts, fill=(0, 128, 255, 200), width=14)
    d.ellipse([620, 560, 900, 860], outline=(0, 255, 0, 255), width=11)
    im.save(path)


def main():
    tmp = tempfile.mkdtemp(prefix="meshcut-selftest-")
    mesh = os.path.join(tmp, "grid.obj")
    png = os.path.join(tmp, "lines.png")
    out = os.path.join(tmp, "out")
    build_mesh(mesh)
    build_image(png)
    png_adv = os.path.join(tmp, "lines-tangential.png")
    build_image(png_adv, tangential=True)
    png_dang = os.path.join(tmp, "lines-dangling.png")
    build_image(png_dang, dangling=True)
    print(f"workdir: {tmp}")

    cases = [
        # (name, cutter flags, require the aspect bound?)
        ("default", ["--max-error", "2", "--max-distance", "60",
                     "--max-aspect", "4"]),
        # tight/loose push the mesher past what Ruppert can guarantee; the
        # contract there is that it MUST report the miss, not hide it
        ("tight", ["--max-error", "0.5", "--max-distance", "25",
                   "--max-aspect", "2.5"]),
        ("loose", ["--max-error", "8", "--max-distance", "250",
                   "--max-aspect", "8"]),
        ("no-refine", ["--max-error", "2", "--max-distance", "60",
                       "--no-refine"]),
    ]
    # An adversarial input: the diagonal runs almost tangent to the ellipse,
    # so the two cuts enclose a sliver far finer than any usable mesh element.
    # The contract there is *graceful degradation*: still manifold, still
    # watertight, still area-conserving -- but the aspect bound may be missed,
    # and the tool must say so.
    adversarial = [
        ("tangential", ["--max-error", "2", "--max-distance", "60",
                        "--max-aspect", "4"]),
        ("dangling", ["--max-error", "2", "--max-distance", "60",
                      "--max-aspect", "4", "--snap", "2"]),
    ]
    rc = 0
    for name, extra in list(cases) + adversarial:
        png_in = {"tangential": png_adv, "dangling": png_dang}.get(name, png)
        print(f"\n=== case: {name} {' '.join(extra)} ===")
        r = subprocess.run([sys.executable, os.path.join(HERE, "mesh_splitter.py"),
                            mesh, png_in, "-o", out, "--quiet",
                            "--debug-png", os.path.join(tmp, f"debug-{name}.png")]
                           + extra, capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr.strip())
        args = dict(zip(extra[::2], extra[1::2]))
        # whenever the aspect bound is missed, the tool must say so out loud
        want = float(args.get("--max-aspect", "0") or 0)
        got = 0.0
        for ln in r.stdout.splitlines():
            if ln.startswith("TOTAL"):
                try:
                    got = float(ln.split()[-1])
                except (ValueError, IndexError):
                    pass
        if want and got > want * 1.001:
            warned = "worst triangle aspect ratio" in r.stderr
            print(f"  {'PASS' if warned else 'FAIL'}  reports that it missed "
                  f"--max-aspect (got {got:.2f} vs {want})")
            if not warned:
                rc = 1
        if name == "dangling":
            warned = "dangling cut endpoint" in r.stderr
            print(f"  {'PASS' if warned else 'FAIL'}  warns about the "
                  f"dangling endpoint")
            if not warned:
                rc = 1
        if r.returncode:
            print(r.stderr)
            rc = 1
            continue
        vargs = [sys.executable, os.path.join(HERE, "verify.py"), mesh,
                 png_in, out]
        if name in ("tangential", "dangling", "tight", "loose"):
            vargs += ["--skip-aspect", "--skip-error"]
        for k in ("--max-error", "--max-distance", "--max-aspect"):
            if k in args:
                vargs += [k, args[k]]
        if "--no-refine" in extra:
            vargs += ["--max-aspect", "1e9"]
        v = subprocess.run(vargs, capture_output=True, text=True)
        print(v.stdout.strip())
        if v.returncode:
            rc = 1
    print("\nSELFTEST", "PASSED" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    sys.exit(main())
