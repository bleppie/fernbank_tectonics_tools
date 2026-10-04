using System;
using Unity.Mathematics;

namespace Tectonic
{
    /// <summary>
    /// A tectonic plate: a set of cells on its own local hex lattice, plus a rigid transform.
    ///
    /// THE UNITY SEAM LIVES HERE. The core never writes <see cref="Position"/>,
    /// <see cref="Rotation"/> or the velocities — Unity's Rigidbody2D owns those and mirrors
    /// them in once per step. The core writes only <see cref="Force"/> and
    /// <see cref="Torque"/>, which Unity reads back and applies. Run the core headlessly and
    /// a simple integrator stands in for Unity; nothing else changes.
    ///
    /// Cells are stored as a dense rectangle in axial space with an occupancy flag, rather
    /// than a hash map. Plate blobs are convex enough that the fill ratio is good, lookup is
    /// a subtraction and an array index, and the two buffers swap by reference each step.
    /// </summary>
    public sealed class Plate
    {
        public int Id;

        /// <summary>
        /// An anchored plate. It participates in everything — collisions, subduction, mountain
        /// building, crust generation — but forces never move it.
        ///
        /// This is what lets small plates a visitor can shove be arranged around an immovable
        /// centre. Mechanically it is just "ignore Force/Torque": the host skips integration
        /// and Unity sets the body kinematic. Nothing in the geology needs to know.
        ///
        /// Note the side effect, which is the point: a static plate is an infinitely heavy
        /// collision partner, so everything that runs into it decelerates hard and builds
        /// relief instead of pushing it away.
        /// </summary>
        public bool IsStatic;

        /// <summary>Stacking rank as well as a physical property: lower density rides on top.
        /// The world keeps its plate list sorted by this, and the sort order is load-bearing
        /// for deciding which plate subducts.</summary>
        public float Density;

        // ---- Written by Unity (Rigidbody2D), read by the core -----------------------
        public float2 Position;
        public float Rotation;
        public float2 LinearVelocity;
        public float AngularVelocity;

        // ---- Written by the core, read by Unity -------------------------------------
        public float2 Force;
        public float Torque;

        // ---- Derived, refreshed each step -------------------------------------------
        public float Mass;
        public float Area;
        public float2 CenterOfMassLocal;
        public float MomentOfInertia;
        public int LiveCellCount;

        /// <summary>
        /// World-space bounding box of the live cells, refreshed once per step.
        ///
        /// This is the collision broad phase. Without it, every cell walks the whole plate
        /// list doing a full inverse-transform and lattice round per candidate, so cost grows
        /// linearly with plate count. With it, a plate that obviously cannot contain the point
        /// is rejected in four float compares — which is what keeps the exhibit free to try
        /// three plates or fifteen without a redesign.
        /// </summary>
        public float2 WorldMin, WorldMax;

        float2 _localMin, _localMax;

        /// <summary>Cheap reject: is this point even inside the plate's bounding box?</summary>
        public bool WorldBoundsContain(float2 p) =>
            p.x >= WorldMin.x && p.x <= WorldMax.x && p.y >= WorldMin.y && p.y <= WorldMax.y;

        /// <summary>
        /// Recompute the world AABB from the cached local one. Plates rotate, so the four
        /// corners of the local box are transformed and re-bounded — conservative, but tight
        /// enough to reject almost everything, and O(1) rather than O(cells).
        /// </summary>
        public void UpdateWorldBounds()
        {
            if (LiveCellCount == 0) { WorldMin = WorldMax = Position; return; }

            float2 lo = new float2(float.MaxValue);
            float2 hi = new float2(float.MinValue);

            for (int i = 0; i < 4; i++)
            {
                float2 corner = new float2(
                    (i & 1) == 0 ? _localMin.x : _localMax.x,
                    (i & 2) == 0 ? _localMin.y : _localMax.y);
                float2 w = LocalToWorld(corner);
                lo = math.min(lo, w);
                hi = math.max(hi, w);
            }

            WorldMin = lo;
            WorldMax = hi;
        }

        /// <summary>Front buffer: the current committed state. Every phase reads this.</summary>
        public Cell[] Read;

        /// <summary>Back buffer: the state being built. The geology gather writes EVERY cell of
        /// it, then <see cref="EndStep"/> swaps the two — so no cell ever observes another
        /// cell's partial update within a step, and no seeding copy is needed.</summary>
        public Cell[] Write;

        /// <summary>Preallocated per-cell scratch for the two-pass geology gather. Reusing
        /// these is what keeps the step allocation-free, which matters far more for an
        /// exhibit running for days than it did for a browser tab.</summary>
        public float[] ScratchElevation;
        public float[] ScratchErodible;
        public float[] ScratchTotalDrop;

        public int2 Origin;
        public int2 Size;

        public Plate(int id, float density, int2 origin, int2 size)
        {
            Id = id;
            Density = density;
            Origin = origin;
            Size = math.max(size, new int2(1, 1));
            int n = Size.x * Size.y;
            Read = new Cell[n];
            Write = new Cell[n];
            for (int i = 0; i < n; i++) { Read[i] = Cell.Empty; Write[i] = Cell.Empty; }
            AllocScratch(n);
        }

        void AllocScratch(int n)
        {
            ScratchElevation = new float[n];
            ScratchErodible = new float[n];
            ScratchTotalDrop = new float[n];
        }

        // ---- Indexing ---------------------------------------------------------------

        /// <summary>Array index for an axial coordinate, or -1 if outside the allocated rect.</summary>
        public int Index(int2 axial)
        {
            int x = axial.x - Origin.x;
            int y = axial.y - Origin.y;
            if ((uint)x >= (uint)Size.x || (uint)y >= (uint)Size.y) return -1;
            return x + y * Size.x;
        }

        public int2 AxialAt(int index) => new int2(
            Origin.x + index % Size.x,
            Origin.y + index / Size.x);

        public bool IsAlive(int2 axial)
        {
            int i = Index(axial);
            return i >= 0 && Read[i].Alive;
        }

        // ---- Transform --------------------------------------------------------------

        public float2 LocalToWorld(float2 local)
        {
            math.sincos(Rotation, out float s, out float c);
            return new float2(local.x * c - local.y * s, local.x * s + local.y * c) + Position;
        }

        public float2 WorldToLocal(float2 world)
        {
            float2 d = world - Position;
            math.sincos(-Rotation, out float s, out float c);
            return new float2(d.x * c - d.y * s, d.x * s + d.y * c);
        }

        public float2 CellWorld(int2 axial, float cellSize) => LocalToWorld(HexLattice.AxialToLocal(axial, cellSize));

        /// <summary>Velocity of a point given in world space: v + ω × r, with the 2D cross
        /// product ω × r = ω · (−r.y, r.x).</summary>
        public float2 VelocityAtWorld(float2 world)
        {
            float2 r = world - Position;
            return LinearVelocity + AngularVelocity * new float2(-r.y, r.x);
        }

        /// <summary>Cell containing a world-space point, or the sentinel if this plate has
        /// nothing there. This is the entire narrow-phase collision test — no colliders, no
        /// spatial index, just an inverse transform and a rounding.</summary>
        public bool TryCellAtWorld(float2 world, float cellSize, out int index, out int2 axial)
        {
            axial = HexLattice.LocalToAxial(WorldToLocal(world), cellSize);
            index = Index(axial);
            if (index < 0 || !Read[index].Alive) { index = -1; return false; }
            return true;
        }

        // ---- Growth -----------------------------------------------------------------

        /// <summary>
        /// Widen the allocated rect so it contains <paramref name="axial"/>, with slack.
        ///
        /// This is the ONLY allocating operation in a step, and it matters more than the byte
        /// count suggests: these cell arrays are megabytes, so they go on the Large Object
        /// Heap, which is only collected with Gen2 and is not compacted by default. Steady
        /// churn there means fragmentation over a multi-day run — the one GC failure mode an
        /// unattended exhibit actually has to care about.
        ///
        /// So the slack grows GEOMETRICALLY with the plate rather than being a fixed ring.
        /// That makes the total number of reallocations logarithmic in final size, so growth
        /// allocation decays to nothing instead of ticking along forever.
        /// </summary>
        public void EnsureContains(int2 axial)
        {
            if (Index(axial) >= 0) return;

            // Try to SHIFT before growing. A plate consumed at one margin and growing at the
            // other keeps the same rect SIZE while its origin drifts, so reallocating is
            // usually unnecessary — and reallocating is what hurts: these arrays are megabytes,
            // so they live on the Large Object Heap, which only Gen2 collects and does not
            // compact by default. Steady LOH churn is the one GC failure mode that actually
            // bites an exhibit running for days.
            if (TryReorigin(axial)) return;

            int slack = math.max(24, math.cmax(Size) / 4);
            int2 lo = math.min(Origin, axial - slack);
            int2 hi = math.max(Origin + Size - 1, axial + slack);
            int2 newOrigin = lo;
            int2 newSize = hi - lo + 1;

            int n = newSize.x * newSize.y;
            var read = new Cell[n];
            var write = new Cell[n];
            for (int i = 0; i < n; i++) { read[i] = Cell.Empty; write[i] = Cell.Empty; }

            for (int y = 0; y < Size.y; y++)
            {
                int srcRow = y * Size.x;
                int dstRow = (Origin.y + y - newOrigin.y) * newSize.x + (Origin.x - newOrigin.x);
                Array.Copy(Read, srcRow, read, dstRow, Size.x);
                Array.Copy(Write, srcRow, write, dstRow, Size.x);
            }

            Read = read;
            Write = write;
            Origin = newOrigin;
            Size = newSize;
            _shift = null;
            AllocScratch(n);
        }

        /// <summary>
        /// Reusable buffer for in-place rect shifting. Allocated once per rect size, then
        /// reused forever — the whole point is to make re-origining allocation-free.
        /// </summary>
        Cell[] _shift;

        /// <summary>
        /// Slide the allocated rect so it covers the live cells plus <paramref name="axial"/>,
        /// keeping the same size and reusing the same arrays. Returns false if the content
        /// genuinely will not fit, in which case the caller reallocates.
        /// </summary>
        bool TryReorigin(int2 axial)
        {
            // Tight bounds of what is actually alive. O(n), but only on the rare step where a
            // plate reaches outside its rect.
            int2 lo = axial;
            int2 hi = axial;
            for (int i = 0; i < Read.Length; i++)
            {
                if (!Read[i].Alive && !Write[i].Alive) continue;
                int2 a = AxialAt(i);
                lo = math.min(lo, a);
                hi = math.max(hi, a);
            }

            const int Pad = 3;
            int2 needed = (hi - lo + 1) + 2 * Pad;
            if (needed.x > Size.x || needed.y > Size.y) return false;

            // Centre the live content in the rect, so the next drift in either direction has
            // room and this stays rare.
            int2 newOrigin = lo - Pad - (Size - needed) / 2;
            if (math.all(newOrigin == Origin)) return false;

            if (_shift == null || _shift.Length != Read.Length) _shift = new Cell[Read.Length];

            ShiftInto(ref Read, newOrigin);
            ShiftInto(ref Write, newOrigin);
            Origin = newOrigin;
            return true;
        }

        void ShiftInto(ref Cell[] src, int2 newOrigin)
        {
            var dst = _shift;
            for (int i = 0; i < dst.Length; i++) dst[i] = Cell.Empty;

            int2 d = Origin - newOrigin;   // a cell at old (x,y) lands at (x + d.x, y + d.y)

            for (int y = 0; y < Size.y; y++)
            {
                int ny = y + d.y;
                if ((uint)ny >= (uint)Size.y) continue;

                int copyLen = Size.x - math.abs(d.x);
                if (copyLen <= 0) continue;

                int srcX = d.x >= 0 ? 0 : -d.x;
                int dstX = d.x >= 0 ? d.x : 0;
                Array.Copy(src, y * Size.x + srcX, dst, ny * Size.x + dstX, copyLen);
            }

            // Swap: the array we just vacated becomes the scratch for next time.
            var old = src;
            src = dst;
            _shift = old;
        }

        /// <summary>
        /// Shrink the allocated rect to hug the live cells.
        ///
        /// Scenario builders size a plate from a bounding box of its intended SHAPE, and for
        /// anything that isn't axis-aligned — a wedge of an annulus, say — that box is mostly
        /// empty. The waste is invisible with two plates and crippling with sixteen: every
        /// per-cell pass — mass properties, zones, boundary flags, collisions, the geology
        /// gather — walks the whole rect, so a half-empty rect is a half-wasted step. This is
        /// what made a six-configuration plate-count sweep run for seventy minutes.
        ///
        /// Called once after construction, so it costs nothing at run time.
        /// </summary>
        public void ShrinkToFit(int pad = 6)
        {
            int2 lo = new int2(int.MaxValue);
            int2 hi = new int2(int.MinValue);
            bool any = false;

            for (int i = 0; i < Read.Length; i++)
            {
                if (!Read[i].Alive) continue;
                int2 a = AxialAt(i);
                lo = any ? math.min(lo, a) : a;
                hi = any ? math.max(hi, a) : a;
                any = true;
            }
            if (!any) return;

            int2 newOrigin = lo - pad;
            int2 newSize = (hi - lo + 1) + 2 * pad;
            if (newSize.x >= Size.x && newSize.y >= Size.y) return;   // already tight

            int n = newSize.x * newSize.y;
            var read = new Cell[n];
            var write = new Cell[n];
            for (int i = 0; i < n; i++) { read[i] = Cell.Empty; write[i] = Cell.Empty; }

            for (int y = 0; y < newSize.y; y++)
            {
                int srcY = newOrigin.y + y - Origin.y;
                if ((uint)srcY >= (uint)Size.y) continue;
                int srcX = newOrigin.x - Origin.x;

                int copyLo = math.max(0, -srcX);
                int copyLen = math.min(newSize.x - copyLo, Size.x - math.max(0, srcX));
                if (copyLen <= 0) continue;

                Array.Copy(Read, srcY * Size.x + srcX + copyLo, read, y * newSize.x + copyLo, copyLen);
                Array.Copy(Write, srcY * Size.x + srcX + copyLo, write, y * newSize.x + copyLo, copyLen);
            }

            Read = read;
            Write = write;
            Origin = newOrigin;
            Size = newSize;
            _shift = null;
            AllocScratch(n);
        }

        // ---- Double buffering -------------------------------------------------------

        public void EndStep()
        {
            var t = Read; Read = Write; Write = t;
        }

        // ---- Bookkeeping ------------------------------------------------------------

        public void RefreshMassProperties(TectonicConfig cfg)
        {
            float cellArea = HexLattice.CellArea(cfg.CellSize);
            float mass = 0f;
            float2 weighted = float2.zero;
            int count = 0;
            float2 localMin = new float2(float.MaxValue);
            float2 localMax = new float2(float.MinValue);

            for (int i = 0; i < Read.Length; i++)
            {
                ref readonly Cell c = ref Read[i];
                if (!c.Alive) continue;
                count++;
                // Continental crust is thicker and lighter but carries more column; the JS
                // model weights it up, and the same weighting keeps collisions feeling right.
                float m = cellArea * (c.Rock.IsContinental ? 3f : 1f);
                mass += m;
                float2 lp = HexLattice.AxialToLocal(AxialAt(i), cfg.CellSize);
                weighted += m * lp;
                localMin = math.min(localMin, lp);
                localMax = math.max(localMax, lp);
            }

            LiveCellCount = count;
            Area = count * cellArea;
            _localMin = localMin;
            _localMax = localMax;
            UpdateWorldBounds();
            Mass = math.max(mass, 1e-4f);
            CenterOfMassLocal = mass > 0f ? weighted / mass : float2.zero;

            float inertia = 0f;
            for (int i = 0; i < Read.Length; i++)
            {
                ref readonly Cell c = ref Read[i];
                if (!c.Alive) continue;
                float m = cellArea * (c.Rock.IsContinental ? 3f : 1f);
                float2 d = HexLattice.AxialToLocal(AxialAt(i), cfg.CellSize) - CenterOfMassLocal;
                inertia += m * math.lengthsq(d);
            }
            MomentOfInertia = math.max(inertia, 1e-4f);
        }
    }
}
