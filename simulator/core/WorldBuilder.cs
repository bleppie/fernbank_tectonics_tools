using Unity.Mathematics;

namespace Tectonic
{
    /// <summary>
    /// Builds starting scenarios. In the exhibit this is also the authoring path for the
    /// default state the system slowly resets to when the room is empty.
    /// </summary>
    public static class WorldBuilder
    {
        /// <summary>
        /// The classic convergent-margin setup, and the smallest arrangement that exercises
        /// the whole geology: an ocean plate driving into a continental one. Produces a
        /// trench, a subducting slab, an arc volcano chain on the overriding plate, and
        /// eventually mountains.
        /// </summary>
        public static World DefaultScenario(TectonicConfig cfg = null)
        {
            cfg ??= new TectonicConfig();
            var world = new World(cfg);

            float2 half = cfg.VisibleHalfExtent;

            // Plate 0 — rides on top. Continental block on the left, ocean to its west.
            var west = MakeRectPlate(world, id: 0, density: 0f,
                center: new float2(-half.x * 0.5f, 0f),
                halfExtent: new float2(half.x * 0.5f, half.y),
                fill: (local, w) =>
                {
                    // Continental crust in the eastern half of this plate: the margin the
                    // ocean plate runs into. The margin is deliberately irregular — a
                    // perfectly straight trench reads as a rendering artifact rather than
                    // as geology.
                    float margin = -half.x * 0.30f
                                 + 7f * math.sin(local.y * 0.075f)
                                 + 3f * math.sin(local.y * 0.21f + 1.3f);
                    return local.x > margin ? CrustKind.Continent : CrustKind.Ocean;
                });

            // Plate 1 — denser, so it goes underneath, and oceanic so it can actually subduct.
            var east = MakeRectPlate(world, id: 1, density: 1f,
                center: new float2(half.x * 0.5f, 0f),
                halfExtent: new float2(half.x * 0.5f, half.y),
                fill: (local, w) =>
                {
                    // A second continent riding in from the east. The ocean between the two
                    // has to subduct away first; when the continents finally meet, neither
                    // column can subduct and the collision builds mountains instead.
                    float margin = -half.x * 0.55f
                                 + 6f * math.sin(local.y * 0.06f + 2.1f)
                                 + 3f * math.sin(local.y * 0.18f);
                    return local.x > margin ? CrustKind.Continent : CrustKind.Ocean;
                });

            // Drive them together. Unity's Rigidbody2D takes over from here.
            west.LinearVelocity = new float2(0.05f, 0f);
            east.LinearVelocity = new float2(-0.35f, 0f);

            world.SortPlatesByDensity();
            for (int i = 0; i < world.Plates.Count; i++)
                world.Plates[i].RefreshMassProperties(cfg);

            return world;
        }

        /// <summary>
        /// The exhibit topology: an ANCHORED central continent with a ring of smaller mobile
        /// plates around it, each one something a visitor can shove.
        ///
        /// The anchored centre is what makes this work as an exhibit. It is an infinitely
        /// heavy collision partner, so a plate pushed against it decelerates hard and spends
        /// that energy building relief rather than shoving the continent across the map. The
        /// result of a visitor's push is a mountain range, not a translation.
        /// </summary>
        public static World RadialScenario(TectonicConfig cfg = null, int satellites = 6)
        {
            cfg ??= new TectonicConfig();
            cfg.Rebuild();
            var world = new World(cfg);

            float2 half = cfg.VisibleHalfExtent;
            const float innerR = 0.44f;   // continent edge, as a fraction of the half-extent
            const float outerR = 1.14f;   // satellites reach into the margin

            // --- the anchored centre -------------------------------------------------
            var centre = MakeRectPlate(world, id: 0, density: 0f,
                center: float2.zero,
                halfExtent: half * (innerR + 0.06f),
                fill: (local, w) =>
                {
                    float2 n = local / (half * innerR);
                    // Irregular coastline: a perfect ellipse reads as a logo, not a continent.
                    float wobble = 1f + 0.09f * math.sin(math.atan2(n.y, n.x) * 3f)
                                      + 0.05f * math.sin(math.atan2(n.y, n.x) * 7f + 1.7f);
                    return math.length(n) <= wobble ? CrustKind.Continent : CrustKind.None;
                });
            centre.IsStatic = true;

            // --- the mobile ring -----------------------------------------------------
            // Unequal sectors. Perfectly regular wedges converging on a disc produce a
            // six-pointed flower — an artifact of the initial condition, not of the physics.
            // Real plate boundaries are irregular, so the default scenario should be too.
            var edges = new float[satellites + 1];
            for (int i = 0; i <= satellites; i++)
            {
                float t = i / (float)satellites;
                float jitter = 0.16f * math.sin(i * 2.399f) + 0.09f * math.sin(i * 5.117f + 1.3f);
                edges[i] = (t + (i == 0 || i == satellites ? 0f : jitter / satellites)) * 2f * math.PI;
            }

            for (int i = 0; i < satellites; i++)
            {
                float a0 = edges[i];
                float a1 = edges[i + 1];

                SectorBounds(half, a0, a1, innerR, outerR, out float2 sectorCentre, out float2 sectorHalf);

                int id = i + 1;
                var plate = MakeRectPlate(world, id: id, density: id,
                    center: sectorCentre,
                    halfExtent: sectorHalf,
                    fill: (local, w) =>
                    {
                        float2 worldPos = sectorCentre + local;
                        float2 n = worldPos / half;
                        float r = math.length(n);
                        // Vary how far each plate starts from the continent, so they don't all
                        // arrive at once and the map keeps changing for longer.
                        float startR = innerR * (1.03f + 0.16f * math.abs(math.sin(id * 1.71f)));
                        if (r < startR || r > outerR) return CrustKind.None;

                        float ang = math.atan2(n.y, n.x);
                        if (ang < 0f) ang += 2f * math.PI;
                        // A hair of clearance so neighbouring satellites don't start overlapped.
                        float pad = 0.012f;
                        if (ang < a0 + pad || ang >= a1 - pad) return CrustKind.None;

                        return CrustKind.Ocean;
                    });

                // Drift inward, toward the anchored continent.
                float mid = (a0 + a1) * 0.5f;
                float speed = 0.20f + 0.16f * math.abs(math.sin(id * 2.31f));
                plate.LinearVelocity = -new float2(math.cos(mid), math.sin(mid)) * speed;
            }

            world.SortPlatesByDensity();
            for (int i = 0; i < world.Plates.Count; i++)
                world.Plates[i].RefreshMassProperties(cfg);

            return world;
        }

        /// <summary>
        /// Tight bounding box of an annulus sector, so each satellite allocates a rect that
        /// fits its own wedge instead of the whole world. At high resolution the difference
        /// is tens of megabytes per plate.
        /// </summary>
        static void SectorBounds(float2 half, float a0, float a1, float innerR, float outerR,
                                 out float2 centre, out float2 halfExtent)
        {
            float2 lo = new float2(float.MaxValue);
            float2 hi = new float2(float.MinValue);

            const int Samples = 48;
            for (int k = 0; k <= Samples; k++)
            {
                float a = math.lerp(a0, a1, k / (float)Samples);
                float2 dir = new float2(math.cos(a), math.sin(a));
                for (int j = 0; j < 2; j++)
                {
                    float2 p = dir * (j == 0 ? innerR : outerR) * half;
                    lo = math.min(lo, p);
                    hi = math.max(hi, p);
                }
            }

            centre = (lo + hi) * 0.5f;
            halfExtent = (hi - lo) * 0.5f + 2f;
        }

        public enum CrustKind { None, Ocean, Continent }

        public delegate CrustKind FillFn(float2 localPos, World world);

        public static Plate MakeRectPlate(World world, int id, float density,
                                          float2 center, float2 halfExtent, FillFn fill)
        {
            var cfg = world.Config;
            float cs = cfg.CellSize;

            HexLattice.AxialBoundsForRect(-halfExtent, halfExtent, cs, out int2 lo, out int2 hi);
            int2 size = hi - lo + 1;

            var plate = new Plate(id, density, lo, size) { Position = center, Rotation = 0f };

            for (int i = 0; i < plate.Read.Length; i++)
            {
                int2 axial = plate.AxialAt(i);
                float2 local = HexLattice.AxialToLocal(axial, cs);
                if (math.any(math.abs(local) > halfExtent)) continue;

                var kind = fill(local, world);
                if (kind == CrustKind.None) continue;

                plate.Read[i] = MakeCell(world, plate.Id, i, kind);
                plate.Write[i] = plate.Read[i];
            }

            plate.ShrinkToFit();
            world.Plates.Add(plate);
            return plate;
        }

        static Cell MakeCell(World world, int plateId, int index, CrustKind kind)
        {
            var cfg = world.Config;
            var rng = world.CellRandom(plateId, index, 0x5BD1E995u);

            Cell c = Cell.Empty;
            c.Alive = true;
            c.MaxCrustThickness = cfg.MaxCrustThicknessBase + rng.NextFloat();
            c.UpliftCapacity = 0.6f + 0.7f * rng.NextFloat();

            if (kind == CrustKind.Ocean)
            {
                float t = cfg.BaseOceanicCrustThickness;
                c.Rock[RockProps.LayerOceanicSediment] = cfg.MaxRegularSedimentThickness;
                c.Rock[RockProps.LayerBasalt] = t * 0.3f;
                c.Rock[RockProps.LayerGabbro] = t * 0.7f;
                c.Age = cfg.PreexistingCrustAge * cfg.MatureCrustAge;
            }
            else
            {
                float t = cfg.BaseContinentalCrustThickness;
                c.Rock[RockProps.LayerContinentalSediment] = 0.03f;
                c.Rock[RockProps.LayerSandstone] = t * 0.15f;
                c.Rock[RockProps.LayerGranite] = t * 0.85f;
                c.Age = cfg.PreexistingCrustAge * cfg.MatureCrustAge;
            }

            return c;
        }
    }
}
