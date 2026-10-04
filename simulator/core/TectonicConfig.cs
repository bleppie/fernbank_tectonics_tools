using Unity.Mathematics;

namespace Tectonic
{
    /// <summary>
    /// All tunable constants in one place. In Unity this is mirrored by a ScriptableObject;
    /// the core never reads global state, so a World always carries its own config.
    ///
    /// Elevation units follow the JS model: 0.5 is sea level, ~1.0 is the highest mountain.
    /// Time units also follow it: 1 unit = 3 million years.
    ///
    /// RESOLUTION INDEPENDENCE
    /// -----------------------
    /// <see cref="CellsAcrossVisibleWidth"/> is the resolution knob, set once at startup.
    /// Everything spatial is authored in WORLD UNITS, never in cells, and the per-cell values
    /// the simulation actually uses are derived in <see cref="Rebuild"/>.
    ///
    /// This distinction is not cosmetic. A constant like "the bending front decays by 0.25 per
    /// neighbour" silently means a different physical distance at every resolution — double
    /// the cell count and the trench gets half as wide. The same applies to subduction zone
    /// width, continent buffers, metamorphic halos, surface transport, and per-cell event
    /// probabilities. Authoring in world units and deriving per-cell values keeps the geology
    /// looking like the same planet at 120 cells across or at 600.
    /// </summary>
    public sealed class TectonicConfig
    {
        // ---- Resolution -------------------------------------------------------------

        /// <summary>
        /// Cells spanning the visible width. THE resolution knob — startup only.
        ///
        /// Cost scales with the square: doubling this quadruples both cell count and step
        /// time. See the benchmark table in the README for measured numbers.
        /// </summary>
        public int CellsAcrossVisibleWidth = 120;

        /// <summary>Hex circumradius in world units. Derived from the resolution.</summary>
        public float CellSize { get; private set; } = 1f;

        /// <summary>Centre-to-centre distance between adjacent cells. Derived.</summary>
        public float CellSpacing { get; private set; } = HexLattice.Sqrt3;

        /// <summary>Area of one cell. Derived.</summary>
        public float CellArea { get; private set; }

        // ---- World extent -----------------------------------------------------------

        /// <summary>Half-extent of the VISIBLE region, in world units.</summary>
        public float2 VisibleHalfExtent = new float2(64f, 36f);   // 128 x 72 world units = 16:9

        /// <summary>
        /// Extra half-extent beyond the visible region before the hard walls. Crust here is
        /// retained but frozen — this is what makes a visitor's actions reversible. It is NOT
        /// destroyed, so pushing a continent off-screen and pulling it back returns it intact.
        /// </summary>
        public float2 MarginHalfExtent = new float2(9.6f, 5.4f);  // 15% of the visible extent

        public float2 WallHalfExtent => VisibleHalfExtent + MarginHalfExtent;

        // ---- Time -------------------------------------------------------------------

        /// <summary>Geology timestep, in model time units. 1 unit = 3 Myr.</summary>
        public float Timestep = 0.04f;

        public const float ModelTimeToMillionYears = 3f;

        // ---- Crust / elevation (thicknesses, so resolution-invariant) ---------------

        public float SeaLevel = 0.5f;
        public float BaseOceanElevation = 0.1f;
        public float BaseContinentElevation = 0.55f;
        public float BaseOceanicCrustThickness = 0.5f;
        public float BaseContinentalCrustThickness = 1.25f;
        public float MaxCrustThicknessBase = 2f;

        /// <summary>
        /// Ceiling on how thick an OCEANIC column can get, from any process.
        ///
        /// Without it, ocean crust in a busy multi-plate jam accumulates arc volcanics and
        /// folding until it breaches sea level and whole ocean plates render as land. Oceanic
        /// lithosphere does not build mountain ranges; it subducts, or it accretes onto a
        /// margin. Kept just below the ~1.17 thickness that would reach sea level.
        ///
        /// LIMITATION: this also prevents genuine volcanic island arcs from emerging. Islands
        /// would need their own mechanism — concentrated uplift at a few cells rather than a
        /// blanket ceiling. Not needed for the first slice.
        /// </summary>
        public float MaxOceanicThickness = 1.1f;

        /// <summary>Fraction of a rock layer's thickness that shows above zero elevation;
        /// the rest is the isostatic root. Oceanic sediment is exempt (it sits on top).</summary>
        public float CrustThicknessToElevation = 0.6f;

        public float CrustBelowZeroElevation => BaseOceanicCrustThickness * CrustThicknessToElevation - BaseOceanElevation;

        public float MaxRegularSedimentThickness = 0.05f;
        public float OceanicRidgeElevation = 0.15f;
        public float SubductionMinElevation = -0.3f;

        // ---- Ages (distance travelled, already world units) -------------------------

        public float MatureCrustAge = 650f / 6371f * 60f;

        /// <summary>
        /// Floor on how fast crust ages, in world-distance-equivalent per unit time.
        ///
        /// Age is "distance travelled" in this model, inherited from the original. That breaks
        /// down for an anchored plate: its speed is zero, so new sea floor at its margins would
        /// stay age-zero forever and render as a permanently bright young ridge. Sea floor
        /// matures with time whether or not the plate happens to be moving.
        /// </summary>
        public float MinAgeRate = 0.15f;
        public float FreshCrustMaxAge = 5.5f;
        public float PreexistingCrustAge = 5.5f;

        // ---- Spatial extents, authored in WORLD UNITS -------------------------------

        /// <summary>How far a slab travels down the zone before it is fully consumed.</summary>
        public float SubductionWidth = 19f;

        /// <summary>How far the "approaching continent" halo reaches into open ocean.</summary>
        public float ContinentBufferWidth = 12f;

        /// <summary>Distance over which the trench depression decays away from the slab.</summary>
        public float BendingFalloffDistance = 7f;

        /// <summary>Distance over which a metamorphic halo decays to ~a third of its intensity.</summary>
        public float MetamorphismHaloDistance = 3f;

        /// <summary>Crust thinning per world unit of continental rifting.</summary>
        public float ContinentalStretchingPerUnit = 0.145f;

        public float TrenchMaxDepth = 0.15f;
        public float MinContinentalCrustThickness = 0.95f;

        /// <summary>Below this thickness a non-subductable layer no longer blocks subduction.</summary>
        public float SubductionBlockingThickness = 0.1f;

        // ---- Geological rates -------------------------------------------------------

        /// <summary>Surface-transport diffusivity, world units squared per unit time.
        /// Resolution-independent; <see cref="ErosionCoefficient"/> is derived from it.</summary>
        public float ErosionDiffusivity = 0.66f;

        /// <summary>Bedrock-to-sediment conversion. A thickness rate, so already invariant.</summary>
        public float WeatheringRate = 0.022f;

        /// <summary>Earthquakes per unit AREA per unit time, so the number visible on screen
        /// doesn't quadruple when the grid gets finer.</summary>
        public float EarthquakeRatePerArea = 0.0077f;

        // Thickness rates: independent of cell size by construction.
        public float SubductionSedimentRate = 0.5f;
        public float FoldRate = 1.1f;
        public float UpliftRate = 0.4f;
        public float VolcanicRockRate = 0.9f;
        public float BasaltGabbroAccretionRate = 0.6f;

        // ---- Forces -----------------------------------------------------------------
        // Force per cell scales with area and so does mass, so acceleration is already
        // resolution-invariant. Nothing to derive here.

        public float DragCoefficient = 0.0060f;
        public float OrogenicDragCoefficient = 0.0300f;

        /// <summary>Outward push from the off-screen spreading ridges at the walls. The
        /// ambient drive that keeps the exhibit alive between visitors.</summary>
        public float RidgePush = 0.0022f;

        /// <summary>Ridge push falls off over this distance from the wall, in world units.</summary>
        public float RidgePushRange = 24f;

        // ---- New crust generation ---------------------------------------------------

        /// <summary>A new cell needs at least this many existing neighbours, so plates grow as
        /// sheets rather than sprouting spindles into open water.</summary>
        public int NewCrustMinNeighbors = 2;

        // ---- Execution --------------------------------------------------------------

        /// <summary>
        /// Run the per-cell geology across threads. Safe by construction: the sweep is a pure
        /// gather, so no cell writes to another and the result is bit-identical either way.
        /// This is the dividend of the Jacobi decision — it costs one line to switch on.
        ///
        /// DEFAULT OFF, and the measurements are worth reading before changing that:
        ///
        ///   - In an isolated geology microbenchmark it is 1.02-1.56x. That number moved twice:
        ///     with the old per-step double-buffer memcpy competing for bandwidth it measured
        ///     0.86-1.10x; removing the memcpy flipped the sign.
        ///   - In a FULL step it is ~1.0x below roughly 300 cells across, because the geology
        ///     is only part of the step and the loop is partly memory-bandwidth-bound.
        ///   - Parallel.For allocates (closures, task infrastructure): measured ~14 KB/step
        ///     against ~600 B/step serial. Gen0 churn, so harmless for an 8-hour run, but it
        ///     is not free and it swamps the allocation figure the harness reports.
        ///
        /// Worth enabling above ~320 cells across, where the geology dominates the step and
        /// the gain reaches ~1.3x. Below that it costs allocation and buys nothing.
        /// </summary>
        public bool ParallelGeology = false;

        // ---- Exhibit behaviour ------------------------------------------------------

        public float StaticMotionThreshold = 0.004f;
        public int StaticStepsBeforeIdle = 600;

        // ---- Determinism ------------------------------------------------------------

        public uint Seed = 0x7EC70127u;

        // ---- Derived per-cell values ------------------------------------------------

        public float SubductionMaxDist { get; private set; }
        public int ContinentBufferCells { get; private set; }
        public float BendingFalloffPerCell { get; private set; }
        public float MetamorphismDecayPerCell { get; private set; }
        public float ErosionCoefficient { get; private set; }
        public float ContinentalStretchPerCell { get; private set; }
        public float EarthquakeProbPerCell { get; private set; }

        /// <summary>
        /// Raised when the derived erosion coefficient had to be clamped for stability. Finer
        /// grids need a proportionally smaller timestep for explicit surface transport; rather
        /// than let an unattended exhibit diverge, the coefficient is capped and this flag is
        /// set so the host can warn.
        /// </summary>
        public bool ErosionClampedForStability { get; private set; }

        /// <summary>
        /// Recompute everything that depends on resolution. The World constructor calls this,
        /// so ordinary use never has to; call it directly only after poking a field.
        /// </summary>
        public void Rebuild()
        {
            // Flat-top hexes step 1.5 * circumradius horizontally.
            float visibleWidth = VisibleHalfExtent.x * 2f;
            CellSize = visibleWidth / (1.5f * math.max(1, CellsAcrossVisibleWidth));
            CellSpacing = HexLattice.CellSpacing(CellSize);
            CellArea = HexLattice.CellArea(CellSize);

            SubductionMaxDist = SubductionWidth;
            ContinentBufferCells = math.max(1, (int)math.round(ContinentBufferWidth / CellSpacing));

            // Decay PER NEIGHBOUR STEP, chosen so the front spans the same world distance
            // regardless of how many cells that happens to be.
            BendingFalloffPerCell = (TrenchMaxDepth + 1f) * (CellSpacing / math.max(BendingFalloffDistance, 1e-4f));
            MetamorphismDecayPerCell = math.exp(-CellSpacing / math.max(MetamorphismHaloDistance, 1e-4f));

            ContinentalStretchPerCell = ContinentalStretchingPerUnit * CellSpacing;

            // Surface transport is diffusion: D = k * spacing^2, so k = D / spacing^2.
            float k = ErosionDiffusivity / (CellSpacing * CellSpacing);

            // Explicit diffusion stability: k * dt * neighbours must stay under ~1. Finer grids
            // hit this first. Clamping keeps the exhibit alive and flags the compromise.
            float kMax = 0.8f / (Timestep * HexLattice.NeighborCount);
            ErosionClampedForStability = k > kMax;
            ErosionCoefficient = math.min(k, kMax);

            EarthquakeProbPerCell = math.saturate(EarthquakeRatePerArea * CellArea);
        }

        public TectonicConfig() => Rebuild();

        public TectonicConfig Clone()
        {
            var c = (TectonicConfig)MemberwiseClone();
            c.Rebuild();
            return c;
        }
    }
}
