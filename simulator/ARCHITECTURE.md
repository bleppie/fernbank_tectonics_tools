# Architecture

Why this code is shaped the way it is.

Most of the decisions here look arbitrary — or worse, look like premature optimisation —
until you know what they are for. So this starts with the exhibit.

---

## 1. The exhibit

A **12-foot-wide top-down projection** of a region of a planet, in a museum gallery.
Visitors push tectonic plates around with their hands and watch the consequences: trenches
open, volcanic arcs light up, mountain ranges grow where continents collide.

The constraints that follow from that, and that shape everything downstream:

| The exhibit is… | So the software must… |
|---|---|
| A **3840 × 2160 projector**, 12 ft wide, 16:9 | Show hexagons small enough not to read as a mosaic — but *not* by brute-force resolution |
| **Driven by visitors pushing plates** | Respond to forces applied at arbitrary points, and make a push produce *geology* rather than translation |
| Arranged as **small mobile plates around an immovable centre** | Support anchored plates that interact fully but never move |
| Something visitors will **push and then pull back** | Make actions reversible — the same crust must return, not fresh sea floor |
| A **region of a larger world**, not a whole planet | Have edges, and a physically honest story for what happens at them |
| Running **unattended, ~8 hours a day**, rebooted each morning | Never diverge, never stutter, never leak; but not need to survive for weeks |
| Expected to **reset itself** when the room is empty | Have a stored default state and a way to detect "nothing is happening any more" |

Two of these deserve emphasis because they are counter-intuitive.

**The anchored centre is what makes a push interesting.** An immovable plate is an
infinitely heavy collision partner, so a plate a visitor shoves into it decelerates hard and
spends that energy building relief. Without it, a push mostly just slides the map around.

**The daily reboot relaxes one requirement and not another.** It removes the *accumulation*
risks — heap growth, Large Object Heap fragmentation, slow drift. It does nothing about
*rate*-driven problems: a garbage collection every ten seconds is just as visible to a
visitor in minute three as in hour seven. That distinction is why allocation is still taken
seriously here despite the short uptime.

---

## 2. Origin, and what was thrown away

This is a port of Concord Consortium's
[Tectonic Explorer](https://github.com/concord-consortium/tectonic-explorer) — TypeScript,
three.js, and a geodesic sphere — to C# on a bounded 2D plane. It is **not a faithful port**.

**Carried over:** the rock taxonomy and its numeric enum values, the crust-column and
isostasy model, the subduction-vs-orogeny decision rules, age-as-distance-travelled, and the
broad shape of a simulation step. That is the part that is genuinely hard-won geology.

**Deleted, and why:**

| Deleted | LOC | Why |
|---|---|---|
| Geodesic sphere (`peels/`), 12 pentagons, peel-seam adjacency, kd-tree, 205 KB Voronoi raster | ~1,650 | All of it existed to answer *"which cell is at this point"* on a sphere. On a flat regular lattice that is a division and a rounding. |
| Physics integrators (Verlet / RK4 / Euler over quaternions) | ~200 | `Rigidbody2D` does this |
| MobX store layer, Web Worker, transferable buffers, output throttling | ~1,900 | Artifacts of running the model on another thread in a browser. With shared memory the renderer reads the model directly. |
| Gauss–Seidel update semantics | — | Replaced with Jacobi. See §5.3 — this is the one replacement that changes *results*, not just structure. |

---

## 3. The shape of the code

```
core/       the simulation. NEVER references UnityEngine. Plain C#, Unity.Mathematics only.
headless/   a Unity-free host: runs the sim, renders PNGs, asserts invariants.
unity/      the bridge: Rigidbody2D, walls, mesh renderer, settings asset.
tests/      28 tests, ~6 s.
```

The `core` / `unity` split is the most important structural decision, and it is not
tidiness. The geology is where the subtle bugs live, and its failures take *thousands of
steps* to appear — a conservation leak at a one-cell rim, a front propagating at the wrong
speed. Those are close to impossible to debug through a running exhibit and trivial to
debug in a headless loop that renders a PNG and checks an invariant.

Enforced three ways: `core/Tectonic.Core.asmdef` sets `noEngineReferences: true`, the
headless project builds the same sources without Unity present, and CI runs both.

---

## 4. A simulation step

`World.Step()` is deliberately a twelve-line table of contents. The phases, in order:

| # | Phase | Writes |
|---|---|---|
| 1 | `RefreshMassProperties` | plate mass, inertia, centre of mass, world AABB |
| 2 | `UpdateZones` | which cells are outside the visible region (frozen) |
| 3 | `UpdateBoundaryFlags` | which cells sit on a plate edge |
| 4 | `DetectCollisions` | transient contact state, **into `Read`** |
| 5 | `UpdateContinentBuffers` | the "approaching continent" halo |
| 6 | `AccumulateForces` | `Plate.Force` / `Plate.Torque`, for the host to apply |
| 7 | `RunGeology` | the Jacobi gather: `Read` → `Write` |
| 8 | `GenerateNewCrust` | divergent boundaries and the wall spreading ridges |

Then the buffers swap.

Note that **contact detection runs before the geology, and writes into `Read`.** Phases
before the gather cannot write `Write` — the gather rebuilds every cell from `Read` and
would silently discard it. That bug shipped once: the "ocean squeezed between two closing
continents is consumed" rule never fired for weeks because it wrote `Write` from phase 4.
Phases that need to affect the next state now flag the cell (`CellFlags.Doomed`) and the
gather honours it.

---

## 5. The five decisions

### 5.1 Everything is a sparse overlay on a hex lattice

A plate is **a set of cells on its own local hex lattice, plus a rigid transform**. Cells
never move relative to their plate; all motion is the transform. Two plates' lattices do not
align in world space — exactly as in the original, where each plate held its own quaternion
over a shared sphere.

Overlap is tested by inverse-transforming a world point into the other plate's frame and
rounding to an axial coordinate: `Plate.TryCellAtWorld`. O(1), allocation-free.

Hex rather than square is deliberate. Most of the geology is diffusion-shaped — sediment
transport, erosion, bending propagation, metamorphic halos. On a 4-neighbour square grid
those grow visible diamond artifacts; on 8-neighbour, square ones. Hex is isotropic enough
that they don't read.

### 5.2 Unity owns integration and walls. The core owns plate-vs-plate.

| | Owner | Why |
|---|---|---|
| Plate motion | `Rigidbody2D` | Unity integrates; the core only emits force and torque |
| Plate ↔ wall | Unity's solver | Non-penetration against a static wall is what a solver is *for* |
| Plate ↔ plate | The core | Plates **must interpenetrate** — subduction is one plate sliding under another, and a contact solver exists to prevent exactly that |

The layer collision matrix enforces the split: `Plates ↔ Walls` on, `Plates ↔ Plates` off.
Each plate's `PolygonCollider2D` is a **convex hull**, not its true outline, rebuilt every
~30 steps — it only ever touches walls, off-screen, where contact geometry is invisible.

`Rigidbody2D.linearDamping` stays at **zero**. Drag is computed per cell in the core because
it depends on each cell's own velocity and area. Doing it in both places damps the plates to
a standstill inside 130 steps — observed, not theorised.

### 5.3 The geology sweep is Jacobi, not Gauss–Seidel

Every operation is a **gather**: a cell computes its new state from its own and its
neighbours' *old* state, and writes only itself. **No cell ever writes into a neighbour.**

The original does the opposite — cells push sediment and metamorphism into their neighbours,
so cells processed later in a step see earlier cells' writes. Its own source comments admit
this makes results depend on `Map` insertion order and diverge visibly over long runs.

Reformulating scatter as gather buys three things:

- **Order independence**, so results are reproducible (tested: `DeterminismTests`)
- **Free parallelism** — the loop can go to `Parallel.For`, or later to Burst, unchanged
- It is simply **easier to reason about** — no question of what a neighbour's state means

Mass is conserved by making every transfer symmetric: a cell's outflow and its neighbour's
matching inflow are computed from the same old state, so they agree exactly.

Randomness is derived from `hash(seed, step, plate, cell)` rather than drawn from a shared
stream — a shared stream would reintroduce order dependence through the back door.

This is the decision that **could not be deferred**. Switching later changes visible
behaviour (fronts propagate at different rates), which would mean re-tuning geology
constants late, against an exhibit that had already been art-directed.

### 5.4 The margin is frozen, not destroyed

The simulated region is larger than the visible one. Crust outside the visible rectangle is
**retained and transported, but its geology is skipped**.

This is the reversibility guarantee. A visitor can shove a continent off-screen and pull it
back, and it returns *bit-identical* rather than as fresh sea floor. Destroying off-screen
crust would be cheaper and would break that.

Memory is bounded by the walls — roughly 1.7× the visible cell count, and it cannot grow
past that. Frozen cells also skip the dominant per-cell cost, so most of the extra area is
free.

Genuinely empty space — a plate retreating from a wall — is filled with new ocean through
the same divergent-boundary path used everywhere else. **The walls are off-screen spreading
ridges.** They also supply a gentle **ridge push**, which is the ambient drive that keeps the
map evolving between visitors, with a physical justification rather than as a fudge.

Without that drive the exhibit has a heat death: walls and drag drain momentum with no
return path, so the map settles. `World.IsSettled` reports when it has — the cue to begin
the slow reset, and a more reliable signal than a presence sensor alone, since a visitor can
stand and watch while everything grinds to a halt.

### 5.5 Spatial constants are authored in world units

`CellsAcrossVisibleWidth` is the resolution knob, set at startup. Everything spatial is
authored in **world units**, and per-cell values are derived in `TectonicConfig.Rebuild()`.

This is not cosmetic. A constant like *"the bending front decays by 0.25 per neighbour"*
silently means a different physical distance at every resolution — double the grid and the
trench comes out half as wide. Same for subduction width, continent buffers, metamorphic
halos, surface transport and per-cell event probabilities.

Two consequences worth knowing:

- Surface transport is a **diffusion**, so its per-cell coefficient goes as `1/spacing²` and
  hits an explicit-stability limit on fine grids. Above ~240 cells across it is **clamped**
  and `ErosionClampedForStability` is raised. The correct fix is a smaller geology timestep,
  which is affordable because the geology tick is decoupled from the physics tick. Clamping
  is the safe default: an unattended exhibit must never diverge.
- Event probabilities scale with **cell area**, so the number of earthquakes visible on
  screen stays constant instead of quadrupling when the grid gets finer.

---

## 6. Choosing a resolution

The instinct at 4K is to raise resolution until hexes stop being visible. **Don't** — that
buys smoothness with CPU, and it is the expensive way. At 320 cells across, one hexagon is
~12 px on the projector, which *flat-shaded* reads as a mosaic. But the mesh does not have to
be flat-shaded.

Corner vertices carry the **average of the three cells meeting there** (`CellPalette.CornerColor`),
so the fragment shader interpolates across cell boundaries and the grid dissolves. Same cell
count, same step cost, no visible hexes. **Visual resolution is decoupled from simulation
resolution**, and the smoothness is paid for by the GPU.

So pick the resolution for the **geology** — how fine a trench, how narrow an arc, how small
an island you need to represent — not for how smooth it looks.

The same palette and interpolation are used by the headless PNG renderer, so what gets
verified offline is what the projector shows. That is why `CellPalette` lives in `core/`
rather than in either renderer: it used to be duplicated, and a duplicated constant elsewhere
(`FoldRate`) had already silently drifted between the two.

---

## 7. Determinism, and where it stops

The core is fully deterministic: same seed, same result, independent of thread count or
iteration order. All three are tested.

**The system is not.** Once `Rigidbody2D` owns integration, Unity's solver varies across
platforms and frame rates.

That has a concrete consequence for the reset behaviour: **the reset must be a state
restore, not a re-simulation.** Replaying inputs from the default scenario will not reproduce
a given map. This is why save/load is the most urgent piece of unbuilt work — see §9.

---

## 8. What has been measured

On four shared container cores, so treat the *shape* as real and the absolute milliseconds as
pessimistic.

**Resolution** (2 plates):

| across | hex @ 12 ft | cells | ms/step | load @ 10 Hz |
|---|---|---|---|---|
| 120 | 1.20" | 9.5 k | 5.2 | 5% |
| 240 | 0.60" | 38 k | 11.1 | 11% |
| 320 | 0.45" | 67 k | 21.8 | 22% |
| 560 | 0.26" | 204 k | 78.7 | 79% |

**Plate count** (240 across) is essentially flat from 3 to 17 plates, once plates cache a
world-space AABB for the collision broad phase. Before that it grew linearly.

**Allocation:** ~600 B/step steady-state. Not zero — the residual is plate-rect resizing —
but the megabyte-scale Large Object Heap churn is gone. `Plate.TryReorigin` slides a plate's
rect within its existing arrays instead of reallocating, because a plate consumed at one
margin and growing at the other keeps its *size* while its *origin* drifts. That took
22.9 KB/step to 612 B/step.

**Threading:** `ParallelGeology` defaults **off**. The isolated geology microbenchmark says
1.02–1.56×, but a *full step* below ~320 across is ~1.0×, and `Parallel.For` allocates
~14 KB/step against ~600 B serial. Worth enabling above ~320, where it reaches ~1.3×.

The loop is **memory-bandwidth-bound, not compute-bound**: `Cell` is ~104 bytes and each
cell touches six neighbours, so a step streams ~700 bytes per cell through cache. More
threads do not buy bandwidth. If performance ever becomes a problem, the lever is
**shrinking the working set** — splitting `Cell` into hot/cold arrays, narrowing fields that
don't need 32 bits — not more parallelism. That is also most of what a Burst conversion
would actually buy.

### On DOTS

Not now, and probably never for the whole project. **There is no production 2D physics in
DOTS** — `com.unity.2d.entities.physics` never left preview and was tied to the discontinued
Project Tiny; `com.unity.physics` is 3D. Going full ECS would forfeit `Rigidbody2D`, which is
the thing Unity is being asked to provide.

The endpoint, if ever, is hybrid: plates stay GameObjects (there are ~10 and their physics
cost is nil), and only the per-cell geology loop becomes jobs. That loop is already shaped
for it — blittable structs, indexed by cell, gather-only. `Cell[]` becomes
`NativeArray<Cell>` and `GatherCell` becomes an `IJobParallelFor` with no change to the
physics.

---

## 9. Known limitations

- **No volcanic island arcs.** Oceanic columns have a blanket thickness ceiling
  (`MaxOceanicThickness`) so that ocean crust in a multi-plate jam doesn't accumulate arc
  volcanics until it breaches sea level and whole plates render as land. That ceiling also
  prevents legitimate islands. They would need concentrated uplift at a few cells rather
  than a blanket cap.
- **The Unity bridge is less verified than the core.** It cannot be exercised by the headless
  harness, so it is reviewed rather than tested.
- **Per-cell randomness is seeded from the array index**, which changes when a plate's rect
  grows or slides. Determinism across identical runs holds, but the same scenario with
  different initial rect padding produces different geology — and any future save/restore
  must preserve the rect exactly.
- **Erosion is clamped above ~240 cells across** (§5.5).

## 10. Deliberately not built

- **Save / restore.** The most urgent, because of §7 — the reset needs it. `Cell` is
  blittable, so this is a straight binary dump.
- **Cross-sections.** In 2D the data path is a line walk reading rock columns — roughly 50
  lines, versus 481 in the original where great-circle sampling made it hard. The expensive
  part is drawing stratigraphy.
- **Reset animation**, presence detection, visitor input hardware.
- **Plate splitting and merging.** The original divides plates by crust age and welds
  co-moving ones. Neither is needed yet.
- **`SurfaceTransport` and `ContactRules` extraction.** Both are pure functions currently
  welded into long methods. Extracting them would let conservation and the subduct/orogeny
  truth table be tested directly instead of by running 400 steps and inspecting a PNG. Worth
  doing the next time the geology is touched.

---

## 11. Before you change anything

Read [`CONTRIBUTING.md`](CONTRIBUTING.md). It is ten lines of rules and a list of things that
look like obvious cleanups but are load-bearing measured decisions — a `List<RockLayer>`
instead of a fixed array, a seed copy between the double buffers, collider-based plate
contact. Each has a measurement attached and a comment explaining itself.
