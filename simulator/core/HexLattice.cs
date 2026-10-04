using Unity.Mathematics;

namespace Tectonic
{
    /// <summary>
    /// Flat-top hexagonal lattice in axial (q, r) coordinates.
    ///
    /// Replaces the JS model's icosahedral geodesic sphere (`peels/`, ~1,657 LOC including
    /// the Voronoi raster and kd-tree). A flat hex grid gives what the geodesic sphere was
    /// approximating — six roughly-equidistant neighbours — with none of the machinery:
    /// no pentagons, no peel seams, no nearest-neighbour acceleration structure.
    ///
    /// Hex rather than square is deliberate. Most of the geology is diffusion-shaped
    /// (sediment spread, erosion, bending propagation, metamorphism). On a 4-neighbour
    /// square grid those produce visible diamond artifacts; on 8-neighbour, square ones.
    /// Hex is isotropic enough that the artifacts don't read.
    ///
    /// Coordinates are always PLATE-LOCAL. Each plate carries its own transform, so two
    /// plates' lattices do not align in world space — exactly as in the JS model, where
    /// every plate held its own quaternion over the shared grid.
    /// </summary>
    public static class HexLattice
    {
        public const int NeighborCount = 6;
        public const float Sqrt3 = 1.7320508075688772f;

        /// <summary>Axial offsets of the six neighbours, counter-clockwise from +q.</summary>
        public static readonly int2[] NeighborDirs =
        {
            new int2( 1,  0),
            new int2( 1, -1),
            new int2( 0, -1),
            new int2(-1,  0),
            new int2(-1,  1),
            new int2( 0,  1),
        };

        public static int2 Neighbor(int2 axial, int dir) => axial + NeighborDirs[dir];

        /// <summary>
        /// The two neighbours sharing corner k of a cell, as indices
        /// into <see cref="NeighborDirs"/>. Corner k sits at angle 60k degrees; the neighbour
        /// directions sit at 30, -30, -90, -150, 150 and 90 degrees, so each corner falls
        /// between exactly two of them.
        ///
        /// Used to give a corner vertex the average of the three cells that meet there, which
        /// makes the mesh interpolate smoothly and dissolves the visible hex grid — visual
        /// smoothness without paying for more cells.
        /// </summary>
        public static readonly int2[] CornerNeighbors =
        {
            new int2(0, 1),
            new int2(0, 5),
            new int2(5, 4),
            new int2(4, 3),
            new int2(3, 2),
            new int2(2, 1),
        };

        /// <summary>Local-space offset of hex corner k, for a cell of the given size.</summary>
        public static float2 CornerOffset(int corner, float size)
        {
            float a = math.radians(60f * corner);
            return new float2(math.cos(a), math.sin(a)) * size;
        }

        /// <summary>Plate-local position of a cell centre. <paramref name="size"/> is the circumradius.</summary>
        public static float2 AxialToLocal(int2 a, float size) => new float2(
            size * 1.5f * a.x,
            size * Sqrt3 * (a.y + a.x * 0.5f));

        public static float2 AxialToLocal(float2 a, float size) => new float2(
            size * 1.5f * a.x,
            size * Sqrt3 * (a.y + a.x * 0.5f));

        /// <summary>
        /// Plate-local position to the containing cell. This is the narrow-phase collision
        /// primitive — the whole reason the sphere needed a 205 KB Voronoi raster, reduced
        /// here to a division and a rounding. O(1), allocation-free.
        /// </summary>
        public static int2 LocalToAxial(float2 p, float size)
        {
            float q = (2f / 3f * p.x) / size;
            float r = (-1f / 3f * p.x + Sqrt3 / 3f * p.y) / size;
            return CubeRound(q, r);
        }

        /// <summary>Round fractional axial coordinates to the nearest cell (via cube coords).</summary>
        public static int2 CubeRound(float q, float r)
        {
            float y = -q - r;

            float rq = math.round(q);
            float ry = math.round(y);
            float rr = math.round(r);

            float dq = math.abs(rq - q);
            float dy = math.abs(ry - y);
            float dr = math.abs(rr - r);

            // Re-derive whichever component drifted most, so q + y + r == 0 still holds.
            if (dq > dy && dq > dr) rq = -ry - rr;
            else if (dy <= dr) rr = -rq - ry;

            return new int2((int)rq, (int)rr);
        }

        /// <summary>Hex (Manhattan-on-cube) distance in cells.</summary>
        public static int Distance(int2 a, int2 b)
        {
            int dq = a.x - b.x;
            int dr = a.y - b.y;
            return (math.abs(dq) + math.abs(dq + dr) + math.abs(dr)) / 2;
        }

        /// <summary>Centre-to-centre spacing of adjacent cells.</summary>
        public static float CellSpacing(float size) => size * Sqrt3;

        /// <summary>Area of one hexagonal cell.</summary>
        public static float CellArea(float size) => 1.5f * Sqrt3 * size * size;

        /// <summary>
        /// Axial bounds that cover an axis-aligned local-space rectangle, with two cells of
        /// slack so edge cells whose centres fall just outside are still included.
        /// </summary>
        public static void AxialBoundsForRect(float2 min, float2 max, float size, out int2 loOut, out int2 hiOut)
        {
            int2 lo = new int2(int.MaxValue, int.MaxValue);
            int2 hi = new int2(int.MinValue, int.MinValue);

            // The axial basis is sheared, so corners of the rect are not the extremes of q/r.
            // Sampling all four corners and padding by two is exact enough and stays cheap.
            for (int i = 0; i < 4; i++)
            {
                float2 corner = new float2((i & 1) == 0 ? min.x : max.x, (i & 2) == 0 ? min.y : max.y);
                int2 a = LocalToAxial(corner, size);
                lo = math.min(lo, a);
                hi = math.max(hi, a);
            }

            loOut = lo - 2;
            hiOut = hi + 2;
        }
    }
}
