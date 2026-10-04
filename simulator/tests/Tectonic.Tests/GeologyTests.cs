using Xunit;

namespace Tectonic.Tests
{
    /// <summary>
    /// The geological rules that decide what a collision produces, and the invariants a long
    /// unattended run depends on.
    /// </summary>
    public class GeologyTests
    {
        static RockColumn Oceanic(float thickness = 0.5f)
        {
            var c = new RockColumn();
            c[RockProps.LayerOceanicSediment] = 0.05f;
            c[RockProps.LayerBasalt] = thickness * 0.3f;
            c[RockProps.LayerGabbro] = thickness * 0.7f;
            return c;
        }

        static RockColumn Continental(float thickness = 1.25f)
        {
            var c = new RockColumn();
            c[RockProps.LayerContinentalSediment] = 0.03f;
            c[RockProps.LayerSandstone] = thickness * 0.15f;
            c[RockProps.LayerGranite] = thickness * 0.85f;
            return c;
        }

        /// <summary>
        /// The single lithological rule the whole simulation turns on: ocean subducts,
        /// continent does not. Everything downstream — trenches, arcs, mountain belts — is a
        /// consequence of this predicate.
        /// </summary>
        [Fact]
        public void OnlyOceanicColumnsCanSubduct()
        {
            const float blocking = 0.1f;
            Assert.True(Oceanic().CanSubduct(blocking));
            Assert.False(Continental().CanSubduct(blocking));
        }

        /// <summary>
        /// Accumulating enough arc volcanics turns an oceanic column non-subductable, which is
        /// what eventually converts a subduction margin into a collision. Worth pinning: it is
        /// the mechanism by which the map stops being static.
        /// </summary>
        [Fact]
        public void VolcanicAccumulationEventuallyBlocksSubduction()
        {
            var column = Oceanic();
            Assert.True(column.CanSubduct(0.1f));

            column[RockProps.LayerAndesite] = 0.4f;
            Assert.False(column.CanSubduct(0.1f));
        }

        /// <summary>
        /// Rock identity is POSITIONAL — the bottom layer decides. This surprises people, so
        /// it is worth a test: adding continental rocks on top of oceanic crust does not make
        /// it continental.
        /// </summary>
        [Fact]
        public void CrustTypeIsDecidedByTheBottomLayer()
        {
            var column = Oceanic();
            column[RockProps.LayerContinentalSediment] = 0.4f;
            column[RockProps.LayerAndesite] = 0.4f;

            Assert.True(column.IsOceanic);
            Assert.False(column.IsContinental);
        }

        [Fact]
        public void RockLayerOrderIsTopToBottom()
        {
            // Continental sediment sits above granite; gabbro is the floor.
            Assert.True(RockProps.LayerContinentalSediment < RockProps.LayerGranite);
            Assert.Equal(RockProps.LayerCount - 1, RockProps.LayerGabbro);
        }

        /// <summary>
        /// Isostasy is faked with a single ratio: part of a column stands proud, the rest is
        /// the root. Oceanic sediment is exempt because it sits on top rather than in it.
        /// </summary>
        [Fact]
        public void OceanicSedimentIsExemptFromTheIsostaticRatio()
        {
            var bare = new RockColumn();
            bare[RockProps.LayerGabbro] = 1f;

            var draped = new RockColumn();
            draped[RockProps.LayerGabbro] = 1f;
            draped[RockProps.LayerOceanicSediment] = 0.1f;

            // The sediment contributes its full thickness; the gabbro only 0.6 of its own.
            Assert.Equal(0.6f, bare.ThicknessAboveZero(0.6f), 4);
            Assert.Equal(0.7f, draped.ThicknessAboveZero(0.6f), 4);
        }

        /// <summary>
        /// A continent driven into an anchored plate must build relief rather than sail
        /// through it. This is the exhibit's whole premise, so it is asserted end to end.
        /// </summary>
        [Fact]
        public void ConvergenceBuildsReliefAndDrivesSubduction()
        {
            var cfg = new TectonicConfig { CellsAcrossVisibleWidth = 60, ParallelGeology = false };
            cfg.Rebuild();
            var world = WorldBuilder.RadialScenario(cfg, satellites: 4);

            bool sawSubduction = false;
            for (int i = 0; i < 500; i++)
            {
                DeterminismTests.SimpleStep(world);
                if (!sawSubduction)
                    foreach (var p in world.Plates)
                        foreach (var c in p.Read)
                            if (c.Alive && c.Subducting) { sawSubduction = true; break; }
            }

            Assert.True(sawSubduction, "no cell ever entered a subduction zone");

            float min = float.MaxValue, max = float.MinValue;
            foreach (var p in world.Plates)
                foreach (var c in p.Read)
                {
                    if (!c.Alive) continue;
                    float e = world.Elevation(c);
                    Assert.False(float.IsNaN(e) || float.IsInfinity(e), "non-finite elevation");
                    Assert.True(c.Rock.Total > 0f, "live cell with an empty rock column");
                    if (e < min) min = e;
                    if (e > max) max = e;
                }

            Assert.True(max - min > 0.4f, $"insufficient relief developed: {max - min:F2}");
        }

        /// <summary>
        /// An anchored plate participates in everything but never moves. A drifting "immovable"
        /// centre is the single worst failure this exhibit could have, and it would be easy to
        /// introduce by accident.
        /// </summary>
        [Fact]
        public void AnchoredPlatesNeverMove()
        {
            var cfg = new TectonicConfig { CellsAcrossVisibleWidth = 60, ParallelGeology = false };
            cfg.Rebuild();
            var world = WorldBuilder.RadialScenario(cfg, satellites: 4);

            Plate anchored = null;
            foreach (var p in world.Plates) if (p.IsStatic) anchored = p;
            Assert.NotNull(anchored);

            var startPos = anchored.Position;
            float startRot = anchored.Rotation;

            for (int i = 0; i < 300; i++) DeterminismTests.SimpleStep(world);

            Assert.Equal(startPos.x, anchored.Position.x, 5);
            Assert.Equal(startPos.y, anchored.Position.y, 5);
            Assert.Equal(startRot, anchored.Rotation, 5);
        }
    }
}
