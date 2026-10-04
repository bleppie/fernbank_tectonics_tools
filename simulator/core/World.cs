using System;
using System.Collections.Generic;
using Unity.Mathematics;

namespace Tectonic
{
    /// <summary>
    /// The simulation. Owns plates, runs the step loop, and knows nothing about Unity.
    ///
    /// Step contract with the host:
    ///   1. host writes each plate's Position/Rotation/LinearVelocity/AngularVelocity
    ///      (from Rigidbody2D, or from the headless integrator)
    ///   2. host calls Step()
    ///   3. host reads each plate's Force/Torque and applies them
    ///
    /// Phases that precede RunGeology write transient contact and zone state back into Read.
    /// RunGeology itself is a strict gather: each cell reads Read and writes only its own slot
    /// in Write, so its result does not depend on iteration order and it runs in parallel
    /// unchanged.
    /// </summary>
    public sealed partial class World
    {
        public readonly TectonicConfig Config;
        public readonly List<Plate> Plates = new List<Plate>();

        public int StepIdx;
        public float Time;

        /// <summary>Consecutive steps with negligible relative plate motion. The exhibit's
        /// idle signal — more reliable than "no visitors", since a visitor can stand and
        /// watch while the map settles.</summary>
        public int StaticSteps;

        /// <summary>Diagnostics for the current step.</summary>
        public int CellsCreatedLastStep;

        /// <summary>Total cells created since the world was built.</summary>
        public int CellsCreatedTotal;

        public bool IsSettled => StaticSteps >= Config.StaticStepsBeforeIdle;

        public float TimeInMillionYears => Time * TectonicConfig.ModelTimeToMillionYears;

        public World(TectonicConfig config)
        {
            Config = config ?? new TectonicConfig();
            Config.Rebuild();   // resolution-derived values must be current before anything reads them
        }

        // ---- Deterministic per-cell randomness --------------------------------------

        /// <summary>
        /// Randomness derived from (seed, step, plate, cell) rather than drawn from a shared
        /// stream. A shared stream would make results depend on iteration order, which is
        /// exactly what Jacobi buffering is here to eliminate — and it is what makes the JS
        /// model's output depend on Map insertion order.
        /// </summary>
        public Unity.Mathematics.Random CellRandom(int plateId, int cellIndex, uint salt = 0u)
        {
            uint h = Config.Seed;
            h = (h ^ (uint)StepIdx) * 16777619u;
            h = (h ^ (uint)plateId) * 16777619u;
            h = (h ^ (uint)cellIndex) * 16777619u;
            h = (h ^ salt) * 16777619u;
            h ^= h >> 15;
            return new Unity.Mathematics.Random(h == 0u ? 1u : h);
        }

        // ---- Derived cell queries ---------------------------------------------------

        public float Elevation(in Cell c)
        {
            float e = c.Rock.ThicknessAboveZero(Config.CrustThicknessToElevation) - Config.CrustBelowZeroElevation;

            if (c.BendingProgress > 0f && c.Rock.IsOceanic)
                e += Config.SubductionMinElevation * c.BendingProgress;

            float na = NormalizedAge(c);
            if (na < 1f)
                e += Config.OceanicRidgeElevation * math.sqrt(1f - na);

            return e;
        }

        public float NormalizedAge(in Cell c) =>
            math.min(Config.FreshCrustMaxAge, c.Age / math.max(Config.MatureCrustAge, 1e-5f));

        public float SubductionProgress(in Cell c)
        {
            if (!c.Subducting) return 0f;
            float t = math.saturate(c.SubductionDist / Config.SubductionMaxDist);
            return t * t;
        }

        // ---- World queries ----------------------------------------------------------

        /// <summary>Topmost plate (lowest density) with a live cell at a world point.</summary>
        public int TopPlateAt(float2 world, out int cellIndex)
        {
            for (int p = 0; p < Plates.Count; p++)
            {
                if (Plates[p].TryCellAtWorld(world, Config.CellSize, out cellIndex, out _))
                    return p;
            }
            cellIndex = -1;
            return -1;
        }

        public bool IsInsideVisible(float2 world) =>
            math.all(math.abs(world) <= Config.VisibleHalfExtent);

        public bool IsInsideWalls(float2 world) =>
            math.all(math.abs(world) <= Config.WallHalfExtent);

        /// <summary>Keeps the plate list sorted so index order is stacking order: index 0 rides
        /// on top. Ties broken by Id so the sort is stable across platforms — List.Sort is not
        /// a stable sort, unlike the JS Array.prototype.sort the original relied on.</summary>
        static readonly Comparison<Plate> ByDensity = (a, b) =>
        {
            int c = a.Density.CompareTo(b.Density);
            return c != 0 ? c : a.Id.CompareTo(b.Id);
        };

        public void SortPlatesByDensity() => Plates.Sort(ByDensity);

        // ---- The step ---------------------------------------------------------------

        public void Step()
        {
            float dt = Config.Timestep;

            SortPlatesByDensity();

            for (int i = 0; i < Plates.Count; i++) Plates[i].RefreshMassProperties(Config);

            UpdateZones();          // frozen flags: which cells are outside the visible region
            UpdateBoundaryFlags();  // which cells sit on a plate edge
            DetectCollisions();     // transient contact state, written into Read
            UpdateContinentBuffers();
            AccumulateForces();     // -> Plate.Force / Plate.Torque, for the host to apply
            RunGeology(dt);         // the Jacobi gather pass: Read -> Write
            GenerateNewCrust(dt);   // divergent boundaries and the wall spreading ridges

            for (int i = 0; i < Plates.Count; i++) Plates[i].EndStep();

            RemoveEmptyPlates();

            StepIdx++;
            Time += dt;

            UpdateSettledCounter();
        }

        void RemoveEmptyPlates()
        {
            for (int i = Plates.Count - 1; i >= 0; i--)
            {
                bool any = false;
                var cells = Plates[i].Read;
                for (int c = 0; c < cells.Length; c++) if (cells[c].Alive) { any = true; break; }
                if (!any) Plates.RemoveAt(i);
            }
        }

        void UpdateSettledCounter()
        {
            float maxRel = 0f;
            for (int a = 0; a < Plates.Count; a++)
                for (int b = a + 1; b < Plates.Count; b++)
                    maxRel = math.max(maxRel, math.length(Plates[a].LinearVelocity - Plates[b].LinearVelocity));

            if (maxRel < Config.StaticMotionThreshold) StaticSteps++;
            else StaticSteps = 0;
        }

        // ---- Zones ------------------------------------------------------------------

        /// <summary>
        /// Mark cells outside the visible region as frozen. They keep every byte of their
        /// state and are still carried by the plate transform, but geology skips them.
        /// That is what makes a visitor's push reversible: shove a continent off-screen,
        /// pull it back, and it returns exactly as it left rather than as fresh sea floor.
        /// It also buys back most of the cost of simulating a larger-than-visible world.
        /// </summary>
        void UpdateZones()
        {
            float cs = Config.CellSize;
            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];
                var read = plate.Read;
                for (int i = 0; i < read.Length; i++)
                {
                    if (!read[i].Alive) continue;
                    bool frozen = !IsInsideVisible(plate.CellWorld(plate.AxialAt(i), cs));
                    read[i].Frozen = frozen;
                }
            }
        }

        void UpdateBoundaryFlags()
        {
            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];
                var read = plate.Read;
                for (int i = 0; i < read.Length; i++)
                {
                    if (!read[i].Alive) continue;
                    int2 a = plate.AxialAt(i);
                    bool boundary = false;
                    for (int d = 0; d < HexLattice.NeighborCount; d++)
                    {
                        if (!plate.IsAlive(a + HexLattice.NeighborDirs[d])) { boundary = true; break; }
                    }
                    read[i].Boundary = boundary;
                }
            }
        }
    }
}
