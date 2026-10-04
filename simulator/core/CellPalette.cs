using Unity.Mathematics;

namespace Tectonic
{
    /// <summary>
    /// The single elevation ramp and overlay set, shared by the Unity mesh renderer and the
    /// headless PNG renderer.
    ///
    /// This lives in the core, not in either renderer, because the project's whole offline
    /// verification story depends on the two agreeing exactly: what gets checked in a headless
    /// run has to be what the projector shows. It used to be duplicated, and a duplicated
    /// constant elsewhere (FoldRate) had already silently drifted between the two — so this is
    /// a copy worth not having.
    ///
    /// Returns linear 0..1 colour. Each host does its own quantisation, which is the only part
    /// that needs host types.
    /// </summary>
    public static class CellPalette
    {
        // Elevation ramp
        public const float SeaDepthRange = 0.7f;
        public const float LandElevationRange = 0.55f;

        // Overlay thresholds and blend strengths
        public const float VolcanicVisibleIntensity = 0.01f;
        public const float MetamorphicVisibleGrade = 0.55f;
        const float SubductionBlend = 0.45f;
        const float VolcanicBlend = 0.5f;
        const float MetamorphicBlend = 0.16f;

        static readonly float3 ShallowSea = new float3(0.28f, 0.55f, 0.78f);
        static readonly float3 DeepSea = new float3(0.02f, 0.06f, 0.22f);
        static readonly float3 Lowland = new float3(0.22f, 0.48f, 0.24f);
        static readonly float3 Upland = new float3(0.52f, 0.42f, 0.24f);
        static readonly float3 Snow = new float3(0.96f, 0.96f, 0.98f);
        static readonly float3 SubductionTint = new float3(0.45f, 0.10f, 0.14f);
        static readonly float3 VolcanicTint = new float3(1.00f, 0.45f, 0.05f);
        static readonly float3 MetamorphicTint = new float3(0.48f, 0.26f, 0.55f);
        static readonly float3 EarthquakeFlash = new float3(1f, 1f, 0.3f);

        /// <summary>
        /// Colour for one cell: an elevation ramp with process overlays composited on top.
        ///
        /// The overlays (subduction, arc volcanism, earthquakes, metamorphism) are drawn over
        /// the terrain rather than encoded into it because the exhibit is about the process,
        /// not the outcome — a visitor should be able to see WHY the mountains are there.
        /// </summary>
        public static float3 ColorOf(World world, in Cell cell)
        {
            float elevation = world.Elevation(cell);
            float seaLevel = world.Config.SeaLevel;

            float3 color;
            if (elevation < seaLevel)
            {
                float depth = math.saturate((seaLevel - elevation) / SeaDepthRange);
                color = math.lerp(ShallowSea, DeepSea, depth);
            }
            else
            {
                float height = math.saturate((elevation - seaLevel) / LandElevationRange);
                color = height < 0.5f
                    ? math.lerp(Lowland, Upland, height * 2f)
                    : math.lerp(Upland, Snow, (height - 0.5f) * 2f);
            }

            if (cell.Subducting) color = math.lerp(color, SubductionTint, SubductionBlend);
            if (cell.VolcanicIntensity > VolcanicVisibleIntensity) color = math.lerp(color, VolcanicTint, VolcanicBlend);
            if (cell.EarthquakeLifespan > 0f) color = EarthquakeFlash;
            if (cell.Metamorphic > MetamorphicVisibleGrade) color = math.lerp(color, MetamorphicTint, MetamorphicBlend);

            return math.saturate(color);
        }

        /// <summary>
        /// Colour for a hexagon corner: the average of the (up to three) cells meeting there.
        ///
        /// This is what dissolves the visible hex grid. Give the corner vertices this value and
        /// the renderer interpolates across cell boundaries, so a 12-pixel hexagon on a 4K
        /// projector reads as smooth terrain rather than a mosaic — visual resolution decoupled
        /// from simulation resolution, paid for by the GPU rather than by more cells.
        /// </summary>
        public static float3 CornerColor(World world, Plate plate, int2 axial, int corner, float3 centerColor)
        {
            float3 sum = centerColor;
            int count = 1;
            int2 pair = HexLattice.CornerNeighbors[corner];

            for (int side = 0; side < 2; side++)
            {
                int dir = side == 0 ? pair.x : pair.y;
                int neighbor = plate.Index(axial + HexLattice.NeighborDirs[dir]);
                if (neighbor < 0 || !plate.Read[neighbor].Alive) continue;
                sum += ColorOf(world, plate.Read[neighbor]);
                count++;
            }

            return sum / count;
        }
    }
}
