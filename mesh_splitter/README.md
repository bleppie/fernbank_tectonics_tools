# mesh_splitter

Cut a 3D OBJ mesh into submeshes along lines drawn on a PNG.

```
uv run mesh_splitter.py mesh.obj lines.png -o out/ --debug-png debug.png
```

Dependencies: `numpy`, `scipy`, `pillow`. All permissively licensed (BSD / HPND)
— no Triangle, no GPL, fine for commercial use.

## Running it

Each script carries a [PEP 723](https://peps.python.org/pep-0723/) inline
metadata block, so with [uv](https://docs.astral.sh/uv/) there is no setup step
at all — no virtualenv to create, nothing to install:

```
uv run mesh_splitter.py mesh.obj lines.png -o out/
uv run verify.py mesh.obj lines.png out/
uv run selftest.py
```

uv reads the dependency list out of the file itself and builds a cached,
throwaway environment. That path resolves fresh each time (fast, but not
pinned). For a **pinned, reproducible** environment use the lock file instead:

```
uv sync                       # creates .venv from uv.lock
uv run python mesh_splitter.py mesh.obj lines.png -o out/
```

The distinction matters: `uv run script.py` honours the script's own inline
block and deliberately ignores `pyproject.toml`; `uv run python script.py` uses
the locked project environment.

Without uv, a plain virtualenv works the same as it always has:

```
python3 -m venv .venv && . .venv/bin/activate
pip install numpy scipy pillow
python mesh_splitter.py mesh.obj lines.png -o out/
```

## What it assumes

* The mesh is a simple 2D manifold; geometry is computed from **x and z** only
  (y is carried along and interpolated, never used for a geometric decision).
* The texture coordinates are a **bijection** onto the image, so each triangle
  gives one affine map between pixel space and world space.
* A pixel belongs to a line iff `alpha > --alpha-threshold`. Colour is ignored.
* The **mesh boundary** closes every region. The image border is *not* treated
  as a line: a cut running off the edge of the image is simply clipped to the
  mesh's UV footprint, which is the edge that actually bounds the surface.

## Parameters

| flag | meaning | default |
|---|---|---|
| `--max-error` | max deviation, **in pixels**, of the simplified cut from the drawn line | 2 |
| `--max-distance` | max spacing between vertices along a cut; longer spans are split at the midpoint, recursively | 100 |
| `--distance-units` | whether `--max-distance` is `world` or `px` | world |
| `--max-aspect` | max longest-edge / shortest-edge ratio of output triangles (1 = equilateral) | 4 |
| `--max-area` | optional max triangle area, world units | off |
| `--no-refine` | skip quality refinement entirely | off |
| `--snap` | weld a dangling cut endpoint onto another cut or the mesh rim within this many pixels; `-1` = auto (~3x the stroke radius) | auto |
| `--weld` | collapse constraint vertices closer than this many pixels into one | 0.05 |
| `--min-feature` | refinement never places vertices closer together than this (world units); `-1` = auto | auto |
| `--min-region-area` | fold regions smaller than this into their largest neighbour; `-1` = auto (1e-6 of the mesh area) | auto |
| `--smooth` | morphological smoothing radius (px) applied to the mask before thinning; `-1` = auto | auto |
| `--prune` | drop thinning spurs shorter than this (px); `-1` = auto | auto |
| `--alpha-threshold` | a pixel is a line if alpha exceeds this | 0 |
| `--no-flip-v` | uv (0,0) is top-left instead of the OBJ default bottom-left | off |
| `--max-points` | safety cap on vertices created by refinement | 500000 |
| `--debug-png` | write an overlay: grey = drawn strokes, colours = regions, black = computed cuts | off |

`--max-aspect` is converted to a minimum angle, `asin(1/A)`: `A=2 → 30°`,
`A=3 → 19.5°`, `A=4 → 14.5°`. Capped at 29° because Ruppert refinement is only
guaranteed to terminate below ~20.7°.

## Output

`out/<name>_0.obj` … `<name>_{n-1}.obj`, ordered by descending area. Each is a
2-manifold with `v`, `vt` (and `vn` if the input had them), with attributes
barycentrically interpolated from the original mesh, and triangle winding
matched to the input. Neighbouring pieces share identical vertices along a cut,
so the set reassembles watertight.

## How it works

1. `alpha > threshold` → binary mask → morphological smoothing.
2. **Zhang–Suen thinning** to 1px centrelines. Contour tracing would return
   both sides of every stroke; thinning returns the line itself.
3. Trace the skeleton into a graph. Diagonal links that merely short-cut an
   orthogonal path are dropped, otherwise every staircase pixel looks like a
   junction. Prune spurs, merge the leftover chains, Douglas–Peucker simplify.
4. Snap dangling endpoints onto another cut or onto the mesh rim.
5. Clip to the mesh's UV footprint; add the mesh boundary as a constraint, with
   a shared vertex wherever a cut meets or ends on it. Pieces lying *on* the
   rim are dropped — they would duplicate it. Only now is it knowable which
   endpoints are really dangling, so that is where the warning comes from.
6. Split each segment where it crosses a mesh triangle edge, then lift to world
   space through that triangle's affine map, so a straight pixel line becomes
   the correct piecewise-linear curve in xz.
7. Densify to `--max-distance`; resolve any crossing or overlapping
   constraints into a proper arrangement.
8. Constrained Delaunay via scipy, recovering constraints by splitting
   encroached segments, plus explicit repair for cocircular input (a regular
   grid lets qhull pick a different valid triangulation than ours).
9. Ruppert refinement for `--max-aspect`.
10. Flood-fill triangles with constraint edges as walls → the regions.
11. Interpolate uv / y / normals from the original mesh and write.

## Results on the sample data

`Topo-Mountain-range_hull.obj` + `Sample_Plates.png` (2049x2049, four strokes
each running off the edge of the image), at defaults: **5 submeshes**, 1554
triangles, ~30 s, worst aspect ratio 3.97 against a requested 4. All nine
independent checks pass. `verify.py` measures the cut vertices at a maximum of
1.53 px from the drawn centreline (mean 0.21 px) against a requested 2 px.

## Failure modes it tells you about

* **Dangling cut endpoints.** A stroke that stops short of another cut or of
  the mesh rim does not separate anything — the flood fill walks around the tip
  and you get one region where you expected two. The tool prints the pixel
  coordinates of every loose end; raise `--snap` to weld them. This is checked
  after clipping, so a stroke that runs off the edge of the image is correctly
  *not* reported: it reaches the rim, which is all that matters.
* **Crossing / overlapping cuts.** A constrained triangulation cannot contain
  two crossing edges. These are resolved into an arrangement up front; if
  something still cannot be recovered, the tool says so rather than splitting
  forever.
* **Non-converging refinement.** Ruppert only terminates for well-conditioned
  input; small angles between constraints can cascade. The tool detects this,
  warns, and stops — the cut is still correct, only the aspect bound is missed.
  Relax `--max-aspect` or pass `--no-refine`. If refinement ever leaves a
  constraint unrecoverable, the whole refinement is discarded and the
  unrefined triangulation is used instead: a correct cut beats a pretty one.
* **The aspect bound is a target, not a guarantee.** Under aggressive
  parameter combinations a handful of triangles can miss it (on the self-test:
  3.27 against a requested 2.5). Whenever that happens the tool prints the
  worst achieved ratio and says so -- it never reports success quietly. The
  summary table shows the achieved aspect ratio per output file.

## A caveat worth stating plainly

Passing every check does not prove the cut is the one you *meant*. Nothing in
the checks knows how many regions you expected. While building this, a change
to the mask smoothing pulled stroke endpoints ~9 px back from the image
border, so one cut stopped reaching it and two regions silently merged -- 5
submeshes became 4, and all nine checks still passed. Look at the
`--debug-png` and count the regions.

## Checking the output

`verify.py` re-derives everything from the written OBJ files and checks it
against the original inputs — it does not trust the cutter:

```
python3 verify.py mesh.obj lines.png out/ --max-aspect 4 --max-distance 100 --max-error 2
```

`selftest.py` generates its own fixtures (no sample data is committed) — it
builds a warped grid mesh whose UV→XZ map is bijective but not
globally affine, draws strokes on it, and runs the cutter and the verifier
across six parameter sets, including two adversarial inputs (two cuts
crossing almost tangentially, and a cut left dangling) where the contract is
graceful degradation plus an honest warning rather than a clean result.

The nine checks: area conservation, 2-manifold with closed boundary loops,
winding matched to the input, uv exactly barycentric, y exactly barycentric,
aspect ratio, cut-edge spacing, centreline error, and watertight sharing of
vertices between neighbouring pieces.
