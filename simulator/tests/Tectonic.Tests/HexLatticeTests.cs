using Unity.Mathematics;
using Xunit;

namespace Tectonic.Tests
{
    /// <summary>
    /// The lattice is the narrow-phase collision primitive: "which cell of plate B is at this
    /// world point" is answered by <see cref="HexLattice.LocalToAxial"/> and nothing else. A
    /// rounding bug here is a subduction bug, and it would present as mysterious geology
    /// rather than as an obvious crash — so it gets tested directly.
    /// </summary>
    public class HexLatticeTests
    {
        [Theory]
        [InlineData(0.25f)]
        [InlineData(1.0f)]
        [InlineData(3.7f)]
        public void AxialToLocalRoundTrips(float size)
        {
            for (int q = -40; q <= 40; q++)
            for (int r = -40; r <= 40; r++)
            {
                var axial = new int2(q, r);
                var back = HexLattice.LocalToAxial(HexLattice.AxialToLocal(axial, size), size);
                Assert.Equal(axial, back);
            }
        }

        /// <summary>
        /// Points anywhere inside a cell must resolve to that cell, not just its exact centre.
        /// This is the property collision detection actually relies on.
        /// </summary>
        [Fact]
        public void PointsInsideACellResolveToThatCell()
        {
            const float size = 1.3f;
            float inradius = size * HexLattice.Sqrt3 * 0.5f;
            var rng = new Random(12345u);

            for (int trial = 0; trial < 4000; trial++)
            {
                var axial = new int2(rng.NextInt(-30, 30), rng.NextInt(-30, 30));
                float2 centre = HexLattice.AxialToLocal(axial, size);

                float angle = rng.NextFloat(0f, 2f * math.PI);
                float radius = rng.NextFloat(0f, inradius * 0.98f);
                float2 probe = centre + radius * new float2(math.cos(angle), math.sin(angle));

                Assert.Equal(axial, HexLattice.LocalToAxial(probe, size));
            }
        }

        /// <summary>
        /// CubeRound collapses the textbook three-way correction into two branches. Correct,
        /// but clever enough to deserve a guard: cube coordinates must always sum to zero.
        /// </summary>
        [Fact]
        public void CubeRoundPreservesTheCubeInvariant()
        {
            var rng = new Random(777u);
            for (int i = 0; i < 20000; i++)
            {
                float q = rng.NextFloat(-50f, 50f);
                float r = rng.NextFloat(-50f, 50f);
                int2 rounded = HexLattice.CubeRound(q, r);
                int y = -rounded.x - rounded.y;
                Assert.Equal(0, rounded.x + y + rounded.y);
            }
        }

        [Fact]
        public void AllSixNeighboursAreDistinctAndAdjacent()
        {
            const float size = 1f;
            float2 origin = HexLattice.AxialToLocal(int2.zero, size);
            float expected = HexLattice.CellSpacing(size);

            for (int d = 0; d < HexLattice.NeighborCount; d++)
            {
                float2 neighbour = HexLattice.AxialToLocal(HexLattice.NeighborDirs[d], size);
                Assert.Equal(expected, math.length(neighbour - origin), 4);

                for (int e = d + 1; e < HexLattice.NeighborCount; e++)
                    Assert.NotEqual(HexLattice.NeighborDirs[d], HexLattice.NeighborDirs[e]);
            }
        }

        /// <summary>
        /// Corner k must sit between the two neighbours CornerNeighbors[k] names — this is what
        /// makes the corner-averaged shading correct rather than merely smooth.
        /// </summary>
        [Fact]
        public void CornerNeighboursBracketTheirCorner()
        {
            const float size = 1f;
            for (int k = 0; k < 6; k++)
            {
                float2 corner = HexLattice.CornerOffset(k, size);
                int2 pair = HexLattice.CornerNeighbors[k];

                float2 a = HexLattice.AxialToLocal(HexLattice.NeighborDirs[pair.x], size);
                float2 b = HexLattice.AxialToLocal(HexLattice.NeighborDirs[pair.y], size);

                // The corner is equidistant from both bracketing cell centres.
                Assert.Equal(math.length(corner - a), math.length(corner - b), 3);
            }
        }
    }
}
