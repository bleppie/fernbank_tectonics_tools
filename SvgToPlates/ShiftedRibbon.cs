using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>
    /// Builds a ribbon as a grid: copies of the shared path shifted to each side, connected column by column.
    ///
    /// Every path vertex (plus extra points on segments longer than interiorSpacing) starts one column. The column's
    /// vertices are the path vertex moved by -ribbonExtent .. +ribbonExtent along a shift direction, so
    ///     uv.x (distance along the path) is exactly the arc length of the path vertex, the same all the way across, and
    ///     uv.y (offset) is exactly the signed shift, the same all the way along.
    /// Both vary linearly over every triangle, so the gradients are clean everywhere, including the ends.
    ///
    /// The shift direction is the path's vertex normal, smoothed along the path over ribbonNormalSmoothing so that
    /// neighbouring columns cannot cross where the path wiggles (ShiftedPath), or the single average normal of the whole
    /// path, perpendicular to the line from its first to its last point (ShiftedAlongAverageNormal).
    /// </summary>
    public static class ShiftedRibbon
    {
        public struct Result
        {
            public float2[] positions;
            public Vector2[] distanceUvs;   // (arc length along the path, signed offset)
            public Vector2[] normalUvs;     // unit normal (x, z), pointing to the +offset side
            public List<int> triangles;
            public int foldedTriangles;     // triangles turned over because columns crossed (should be 0)
        }

        /// <summary>
        /// path runs clockwise around the lower-numbered plate, so its left side (positive offset) is away from that plate,
        /// matching the Outline method.
        /// </summary>
        public static Result Build(float2[] path, PlateMeshSettings settings)
        {
            float2[] columns = Subdivide(path, settings.interiorSpacing);
            int columnCount = columns.Length;
            float[] arcLength = PathGeometry.CumulativeLength(columns);

            float2[] normals = SmoothedNormals(columns, arcLength, settings.ribbonNormalSmoothing);
            var shift = new float2[columnCount];
            if (settings.ribbonMethod == RibbonMethod.ShiftedAlongAverageNormal)
            {
                float2 average = PlateGeometry.Perp(math.normalize(columns[columnCount - 1] - columns[0]));
                for (int i = 0; i < columnCount; i++) shift[i] = average;
            }
            else
            {
                normals.CopyTo(shift, 0);
            }

            // Rows: -extent .. +extent, with the path itself as the middle row.
            int perSide = settings.interiorSpacing > 0f
                ? math.max(1, (int)math.ceil(settings.ribbonExtent / settings.interiorSpacing))
                : 1;
            int rowCount = 2 * perSide + 1;

            var result = new Result
            {
                positions = new float2[columnCount * rowCount],
                distanceUvs = new Vector2[columnCount * rowCount],
                normalUvs = new Vector2[columnCount * rowCount],
                triangles = new List<int>(6 * (columnCount - 1) * (rowCount - 1))
            };

            for (int i = 0; i < columnCount; i++)
            {
                for (int k = 0; k < rowCount; k++)
                {
                    float offset = settings.ribbonExtent * (k - perSide) / perSide;
                    int v = i * rowCount + k;
                    result.positions[v] = columns[i] + shift[i] * offset;
                    result.distanceUvs[v] = new Vector2(arcLength[i], offset);
                    result.normalUvs[v] = new Vector2(normals[i].x, normals[i].y);
                }
            }

            for (int i = 0; i < columnCount - 1; i++)
            {
                for (int k = 0; k < rowCount - 1; k++)
                {
                    int a = i * rowCount + k;      // this column, this row
                    int b = a + rowCount;          // next column, same row
                    int c = b + 1;                 // next column, next row
                    int d = a + 1;                 // this column, next row
                    result.triangles.AddRange(new[] { a, b, c, a, c, d });
                }
            }

            OrientFacingUp(result.positions, result.triangles, out result.foldedTriangles);
            return result;
        }

        /// <summary>The path's vertices, with extra points on segments longer than spacing (spacing 0: none).</summary>
        static float2[] Subdivide(float2[] path, float spacing)
        {
            var points = new List<float2> { path[0] };
            for (int i = 0; i < path.Length - 1; i++)
            {
                float length = math.distance(path[i], path[i + 1]);
                int pieces = spacing > 0f ? math.max(1, (int)math.ceil(length / spacing)) : 1;
                for (int j = 1; j <= pieces; j++)
                    points.Add(math.lerp(path[i], path[i + 1], j / (float)pieces));
            }
            return points.ToArray();
        }

        /// <summary>
        /// Unit left normal at every vertex (the tangent from TangentField rotated 90 degrees counter-clockwise),
        /// Gaussian-averaged along the path with standard deviation 'smoothing' (world units, measured along the path).
        /// </summary>
        static float2[] SmoothedNormals(float2[] path, float[] arcLength, float smoothing)
        {
            float2[] tangents = TangentField.ComputePathTangents(path);
            int n = path.Length;
            var normals = new float2[n];
            for (int i = 0; i < n; i++) normals[i] = PlateGeometry.Perp(tangents[i]);
            if (smoothing <= 0f) return normals;

            var smoothed = new float2[n];
            for (int i = 0; i < n; i++)
            {
                double2 sum = double2.zero;
                for (int j = 0; j < n; j++)
                {
                    double d = (arcLength[j] - arcLength[i]) / smoothing;
                    if (d * d > 36.0) continue;   // beyond 6 standard deviations
                    sum += math.exp(-0.5 * d * d) * (double2)normals[j];
                }
                smoothed[i] = math.normalizesafe((float2)sum, normals[i]);
            }
            return smoothed;
        }

        /// <summary>
        /// Flips every triangle if needed so the ribbon faces +y (Unity front faces are clockwise seen from above, i.e. a
        /// negative 2D signed area). Triangles left with the other orientation are folded over: columns crossed.
        /// </summary>
        static void OrientFacingUp(float2[] positions, List<int> triangles, out int folded)
        {
            double total = 0.0;
            for (int t = 0; t < triangles.Count; t += 3)
                total += PlateGeometry.Cross(positions[triangles[t + 1]] - positions[triangles[t]],
                                             positions[triangles[t + 2]] - positions[triangles[t]]);
            if (total > 0.0)
            {
                for (int t = 0; t < triangles.Count; t += 3)
                    (triangles[t + 1], triangles[t + 2]) = (triangles[t + 2], triangles[t + 1]);
            }

            folded = 0;
            for (int t = 0; t < triangles.Count; t += 3)
            {
                if (PlateGeometry.Cross(positions[triangles[t + 1]] - positions[triangles[t]],
                                        positions[triangles[t + 2]] - positions[triangles[t]]) >= 0f)
                    folded++;
            }
        }
    }
}
