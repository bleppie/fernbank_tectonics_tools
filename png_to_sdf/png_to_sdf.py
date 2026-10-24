#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.22",
#     "scipy>=1.10",
#     "opencv-python-headless>=4.5",
#     "OpenEXR>=3.2",
# ]
# ///
"""
png_to_sdf.py -- signed distance field of a PNG's alpha mask, written as a
single-channel 32-bit float EXR.

Pixels with alpha > 0.5 are "inside". Each output pixel is its distance (in
pixels) to the mask's edge, negative inside and positive outside, divided by
max(width, height). The edge sits halfway between an inside and an outside
pixel, so it's where the field crosses 0.

Usage:
    python png_to_sdf.py shape.png                 # writes shape_sdf.exr
    python png_to_sdf.py shape.png -o out.exr
    python png_to_sdf.py shape.png --positive-inside
    uv run png_to_sdf.py shape.png                 # installs the libraries for you
"""

import argparse
import os
import sys

import cv2
import numpy as np
import OpenEXR
from scipy.ndimage import distance_transform_edt


def read_alpha(path: str) -> np.ndarray:
    """Alpha channel as floats in 0..1 (8- or 16-bit PNGs)."""
    img = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if img is None:
        sys.exit(f"Couldn't read {path}")
    if img.ndim != 3 or img.shape[2] not in (2, 4):
        sys.exit(f"{path} has no alpha channel")
    alpha = img[:, :, -1]
    return alpha.astype(np.float64) / np.iinfo(alpha.dtype).max


def signed_distance(mask: np.ndarray) -> np.ndarray:
    """Distance in pixels to the mask edge: negative inside, positive outside.

    distance_transform_edt gives each nonzero pixel's distance to the nearest
    zero pixel (centre to centre). Taking 0.5 off puts the edge halfway
    between the last inside and first outside pixel.
    """
    if mask.all() or not mask.any():
        sys.exit("The alpha mask is all inside or all outside; there's no edge to measure from.")
    outside = distance_transform_edt(~mask)   # outside pixels: distance to nearest inside pixel
    inside = distance_transform_edt(mask)     # inside pixels: distance to nearest outside pixel
    return np.where(mask, -(inside - 0.5), outside - 0.5)


def write_exr(path: str, data: np.ndarray, channel: str):
    header = {"compression": OpenEXR.ZIP_COMPRESSION, "type": OpenEXR.scanlineimage}
    with OpenEXR.File(header, {channel: data.astype(np.float32)}) as f:
        f.write(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("png", help="input PNG with an alpha channel")
    parser.add_argument("-o", "--output", help="output EXR (default: <input>_sdf.exr)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="alpha above this is inside (default: 0.5)")
    parser.add_argument("--positive-inside", action="store_true",
                        help="make inside positive and outside negative")
    parser.add_argument("--channel", default="Y",
                        help="EXR channel name (default: Y, read as greyscale by "
                             "most tools; use R if your tool expects red)")
    args = parser.parse_args(argv)

    alpha = read_alpha(args.png)
    mask = alpha > args.threshold
    h, w = mask.shape
    sdf = signed_distance(mask) / max(h, w)
    if args.positive_inside:
        sdf = -sdf

    out = args.output or os.path.splitext(args.png)[0] + "_sdf.exr"
    write_exr(out, sdf, args.channel)
    print(f"{args.png}: {w} x {h}, {mask.mean() * 100:.1f}% inside")
    print(f"wrote {out}: channel {args.channel}, float32, range "
          f"{sdf.min():.5f} .. {sdf.max():.5f} (pixels / {max(h, w)})")


if __name__ == "__main__":
    main()
