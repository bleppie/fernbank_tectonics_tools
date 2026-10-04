# Contributing

Three rules. They are short because they are the ones that, if broken, produce bugs that
look like geology rather than like errors.

## 1. `core/` must never reference `UnityEngine`

The simulation is a plain C# library so it can be built, run and tested without Unity. That
is not tidiness — it is how the geology gets verified at all, since the interesting failures
take thousands of steps to appear and are impossible to debug through a running exhibit.

`core/Tectonic.Core.asmdef` sets `noEngineReferences: true`, and the headless build fails
loudly if this is violated. `Unity.Mathematics` is fine (it is pure C#); `UnityEngine` is not.

## 2. Geology is gather-only — never write a neighbour

Inside `World.Geology.cs`, a cell computes its new state by reading its own and its
neighbours' **old** state (`Read`), and writes **only itself** (`Write[i]`).

Never write `Write[neighbourIndex]`. Never mutate `Read` during the gather.

This is what makes results independent of iteration order, reproducible, and safely
parallel. The upstream project does the opposite and its own source comments admit the
results drift visibly over long runs.

If a phase running *before* the gather needs to affect the next state, flag the cell in
`Read` (see `CellFlags.Doomed`) and let the gather honour it. Writing `Write` from an earlier
phase does not work — the gather rebuilds the whole cell from `Read` and will silently
discard it. That bug shipped once already.

## 3. Run the harness before pushing

```sh
dotnet test tests/Tectonic.Tests          # 28 tests, ~6 s
./run-headless.sh --steps 2500 --scenario radial   # invariants + PNG maps
```

The headless run exits non-zero if any invariant fails. Both run in CI.

## Things that look like cleanups but aren't

These have each been measured, and each is documented where it lives. Please read the
comment before "fixing" it:

- **`RockColumn` as a `fixed float[11]`** rather than a `List<RockLayer>`. A crust column
  holds at most one layer per rock type, so the list collapses to a fixed vector. Avoids
  roughly a million allocations per second.
- **`Cell` as a blittable struct**, not a class.
- **No seed copy between the double buffers.** The gather writes every cell, so `Write`
  needs no seeding. Adding a copy back would mask exactly the bug described in rule 2.
- **Plate-vs-plate contact is hand-rolled, not collider-based.** Subduction *requires*
  interpenetration; a contact solver exists to prevent it. Unity's physics is used for
  plate-vs-wall only, enforced by the layer collision matrix.
- **Spatial constants are authored in world units**, with per-cell values derived in
  `TectonicConfig.Rebuild()`. A constant expressed "per cell" silently changes meaning at
  every resolution.

## Style

Follow what is there: explicit names over clever ones, comments that explain *why* rather
than *what*, and a measurement attached to any claim about performance.
