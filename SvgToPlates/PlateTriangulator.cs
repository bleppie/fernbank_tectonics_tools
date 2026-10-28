using System;
using System.Collections.Generic;
using andywiecko.BurstTriangulator;   // https://github.com/andywiecko/BurstTriangulator (package: com.andywiecko.burst.triangulator)
using Unity.Collections;
using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>
    /// The only file that talks to the triangulation package. If the package API changes (or you swap in another
    /// triangulator) this is the one place to edit. Written against BurstTriangulator 2.5+/3.x, where positions are
    /// double2; on older 2.x versions change double2 to float2 below.
    /// </summary>
    public static class PlateTriangulator
    {
        /// <summary>Vertices and triangles returned by the triangulator (refinement can add vertices).</summary>
        public struct Result
        {
            public List<float2> positions;   // the input points first, in the same order, then any added points
            public int[] triangles;          // three indices per triangle
        }

        /// <summary>
        /// Constrained Delaunay triangulation of the convex hull of the points, optionally refined.
        /// constraintEdgePairs holds index pairs (a0, b0, a1, b1, ...) into points; those edges appear in the result
        /// (refinement may split them). Without refinement, triangles outside the outline are still present and are
        /// removed afterwards by KeepInsideTrianglesFacingUp. With refinement, RestoreBoundary is enabled so the
        /// outside is removed before refining and no vertices are wasted there.
        /// allowRefinement = false ignores settings.refineMesh (plain Delaunay of the points, nothing added).
        /// </summary>
        public static Result Triangulate(IReadOnlyList<float2> points, IReadOnlyList<int> constraintEdgePairs, PlateMeshSettings settings,
                                         bool allowRefinement = true)
        {
            bool refine = allowRefinement && settings.refineMesh;

            // Fill managed arrays first: a 'using var' variable is read-only, so its elements cannot be assigned.
            var positionData = new double2[points.Count];
            for (int i = 0; i < positionData.Length; i++)
                positionData[i] = points[i];

            var constraintData = new int[constraintEdgePairs.Count];
            for (int i = 0; i < constraintData.Length; i++)
                constraintData[i] = constraintEdgePairs[i];

            using var positions = new NativeArray<double2>(positionData, Allocator.Persistent);
            using var constraints = new NativeArray<int>(constraintData, Allocator.Persistent);

            using var triangulator = new Triangulator(Allocator.Persistent)
            {
                Input = { Positions = positions, ConstraintEdges = constraints },
                Settings =
                {
                    RefineMesh = refine,
                    RestoreBoundary = refine,
                    RefinementThresholds =
                    {
                        Angle = math.radians(settings.minimumAngleDegrees),
                        Area = settings.maximumTriangleArea
                    },
                    ConcentricShellsParameter = settings.concentricShellsParameter
                }
            };
            if (constraintData.Length == 0)
                triangulator.Input.ConstraintEdges = default;   // no constraints: plain Delaunay
            triangulator.Run();

            var result = new Result
            {
                positions = new List<float2>(triangulator.Output.Positions.Length),
                triangles = triangulator.Output.Triangles.AsArray().ToArray()
            };
            for (int i = 0; i < triangulator.Output.Positions.Length; i++)
            {
                double2 p = triangulator.Output.Positions[i];
                result.positions.Add(new float2((float)p.x, (float)p.y));
            }

            // The attribute code relies on the first points keeping their index.
            for (int i = 0; i < points.Count; i++)
            {
                if (math.distance(result.positions[i], points[i]) > 1e-4f)
                    throw new InvalidOperationException("The triangulator reordered or moved the input points (is a preprocessor enabled?).");
            }
            return result;
        }

        /// <summary>
        /// Keeps only the triangles whose centroid lies inside the outline, and flips them if needed so every
        /// triangle faces +y (Unity front faces are clockwise seen from above; with x to the right and z up that is a
        /// negative 2D signed area). Because the outline edges are constraints, a triangle is either wholly inside or
        /// wholly outside, so testing the centroid is enough.
        /// </summary>
        public static List<int> KeepInsideTrianglesFacingUp(IReadOnlyList<float2> points, int[] triangles, float2[] outline)
        {
            var kept = new List<int>(triangles.Length);
            for (int t = 0; t < triangles.Length; t += 3)
            {
                int a = triangles[t], b = triangles[t + 1], c = triangles[t + 2];
                float2 pa = points[a], pb = points[b], pc = points[c];

                float2 centroid = (pa + pb + pc) / 3f;
                if (!PlateGeometry.IsInside(centroid, outline)) continue;

                if (PlateGeometry.Cross(pb - pa, pc - pa) > 0f)
                    (b, c) = (c, b);

                kept.Add(a);
                kept.Add(b);
                kept.Add(c);
            }
            return kept;
        }
    }
}
