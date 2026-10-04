using Unity.Mathematics;
using Xunit;

namespace Tectonic.Tests
{
    /// <summary>
    /// Resolution independence is the design claim that lets the exhibit's grid be retuned
    /// without retuning the geology. It is easy to break by accident — anyone adding a
    /// constant expressed "per cell" instead of "per world unit" breaks it silently, and the
    /// symptom is a trench that looks slightly wrong at one resolution.
    /// </summary>
    public class ConfigTests
    {
        [Theory]
        [InlineData(120)]
        [InlineData(240)]
        [InlineData(480)]
        public void WorldDistancesSurviveResolutionChanges(int cellsAcross)
        {
            var cfg = new TectonicConfig { CellsAcrossVisibleWidth = cellsAcross };
            cfg.Rebuild();

            // The continent buffer is authored in world units; its per-cell form must reproduce
            // that distance to within one cell.
            float bufferWorld = cfg.ContinentBufferCells * cfg.CellSpacing;
            Assert.InRange(bufferWorld,
                cfg.ContinentBufferWidth - cfg.CellSpacing,
                cfg.ContinentBufferWidth + cfg.CellSpacing);

            // Likewise the bending falloff: decaying by BendingFalloffPerCell each neighbour
            // step must span BendingFalloffDistance in world units.
            float falloffWorld = (cfg.TrenchMaxDepth + 1f) / cfg.BendingFalloffPerCell * cfg.CellSpacing;
            Assert.InRange(falloffWorld,
                cfg.BendingFalloffDistance * 0.9f,
                cfg.BendingFalloffDistance * 1.1f);
        }

        [Fact]
        public void CellsAcrossActuallySpansTheVisibleWidth()
        {
            foreach (int across in new[] { 60, 120, 320, 560 })
            {
                var cfg = new TectonicConfig { CellsAcrossVisibleWidth = across };
                cfg.Rebuild();
                float spanned = across * 1.5f * cfg.CellSize;       // flat-top column pitch
                Assert.Equal(cfg.VisibleHalfExtent.x * 2f, spanned, 3);
            }
        }

        [Fact]
        public void RebuildIsIdempotent()
        {
            var cfg = new TectonicConfig { CellsAcrossVisibleWidth = 240 };
            cfg.Rebuild();
            float size = cfg.CellSize, erosion = cfg.ErosionCoefficient;
            int buffer = cfg.ContinentBufferCells;

            cfg.Rebuild();

            Assert.Equal(size, cfg.CellSize);
            Assert.Equal(erosion, cfg.ErosionCoefficient);
            Assert.Equal(buffer, cfg.ContinentBufferCells);
        }

        /// <summary>
        /// Explicit surface transport is only stable while k*dt*neighbours stays below ~1.
        /// Finer grids hit the limit first, so the coefficient is clamped rather than allowed
        /// to diverge — an unattended exhibit must never blow up.
        /// </summary>
        [Fact]
        public void ErosionCoefficientNeverExceedsTheStabilityLimit()
        {
            foreach (int across in new[] { 60, 120, 240, 320, 420, 560, 800 })
            {
                var cfg = new TectonicConfig { CellsAcrossVisibleWidth = across };
                cfg.Rebuild();

                Assert.True(cfg.ErosionCoefficient * cfg.Timestep * HexLattice.NeighborCount <= 0.81f,
                    $"unstable erosion coefficient at {across} cells across");

                float unclamped = cfg.ErosionDiffusivity / (cfg.CellSpacing * cfg.CellSpacing);
                float limit = 0.8f / (cfg.Timestep * HexLattice.NeighborCount);
                Assert.Equal(unclamped > limit, cfg.ErosionClampedForStability);
            }
        }
    }

    /// <summary>
    /// Plate cell storage: a dense rect in axial space that grows, shrinks and slides beneath
    /// the simulation. It holds the hairiest index arithmetic in the project and is the only
    /// code in a step that allocates, so it is worth pinning down.
    /// </summary>
    public class PlateStorageTests
    {
        static Plate MakePlate()
        {
            var plate = new Plate(id: 0, density: 0f, origin: new int2(-10, -10), size: new int2(24, 24));
            for (int q = -4; q <= 4; q++)
            for (int r = -4; r <= 4; r++)
            {
                int i = plate.Index(new int2(q, r));
                var cell = Cell.Empty;
                cell.Alive = true;
                cell.Age = q * 100f + r;                       // a checksummable marker
                cell.Rock[RockProps.LayerGabbro] = 0.5f;
                plate.Read[i] = cell;
                plate.Write[i] = cell;
            }
            return plate;
        }

        static double Checksum(Plate p)
        {
            double sum = 0;
            for (int i = 0; i < p.Read.Length; i++)
            {
                if (!p.Read[i].Alive) continue;
                int2 a = p.AxialAt(i);
                sum += a.x * 31.0 + a.y * 17.0 + p.Read[i].Age * 0.5;
            }
            return sum;
        }

        [Fact]
        public void IndexAndAxialAreInverses()
        {
            var p = MakePlate();
            for (int i = 0; i < p.Read.Length; i++)
                Assert.Equal(i, p.Index(p.AxialAt(i)));
        }

        [Fact]
        public void GrowingPreservesEveryLiveCell()
        {
            var p = MakePlate();
            double before = Checksum(p);

            p.EnsureContains(new int2(200, 200));   // far outside: forces a real reallocation

            Assert.True(p.Index(new int2(200, 200)) >= 0);
            Assert.Equal(before, Checksum(p), 3);
        }

        /// <summary>
        /// The re-origin path is the 37x allocation win the design leans on: when a plate is
        /// consumed at one margin and grows at the other, its rect only needs to SLIDE. If
        /// this silently regressed to reallocating, nothing else would notice — the results
        /// would be identical and only the allocation rate would change.
        /// </summary>
        [Fact]
        public void DriftSlidesTheRectInsteadOfReallocating()
        {
            var p = MakePlate();

            // Reproduce the situation the optimisation exists for: a plate consumed at one
            // margin and growing at the other. Its live content moves across the lattice while
            // its SIZE stays put, so the rect only needs to slide. (Simply asking for a distant
            // cell would not exercise this — the content would no longer fit, and the fallback
            // reallocation is the correct answer there.)
            for (int i = 0; i < p.Read.Length; i++)
            {
                if (p.AxialAt(i).x >= 0) continue;
                p.Read[i] = Cell.Empty;
                p.Write[i] = Cell.Empty;
            }

            double before = Checksum(p);
            var originalRead = p.Read;

            p.EnsureContains(new int2(16, 0));

            Assert.True(p.Index(new int2(16, 0)) >= 0);

            // The arrays were reused, not reallocated: after sliding both buffers through the
            // single scratch, the array Read used to occupy is now Write.
            Assert.Same(originalRead, p.Write);
            Assert.Equal(before, Checksum(p), 3);
        }

        [Fact]
        public void ShrinkToFitKeepsEveryLiveCell()
        {
            var p = MakePlate();
            double before = Checksum(p);
            int sizeBefore = p.Size.x * p.Size.y;

            p.ShrinkToFit();

            Assert.True(p.Size.x * p.Size.y < sizeBefore);
            Assert.Equal(before, Checksum(p), 3);
        }
    }
}
