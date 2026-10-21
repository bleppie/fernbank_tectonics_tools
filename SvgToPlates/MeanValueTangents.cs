using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>
    /// Mean value coordinates (Floater 2003; Hormann and Floater 2006, "Mean value coordinates for arbitrary planar
    /// polygons"). For a point inside a simple polygon they give weights, one per polygon vertex, that sum to 1,
    /// reproduce linear data exactly, vary smoothly, and interpolate boundary values (linearly along each edge).
    /// They do not depend on any triangulation. They are only meaningful for points inside the polygon.
    /// </summary>
    public static class MeanValueTangents
    {
        /// <summary>
        /// Weighted sum of the vertex values using the mean value coordinates of p. The result is not normalised.
        /// p must lie inside the polygon and not on its boundary.
        /// </summary>
        public static float2 Interpolate(float2 p, float2[] polygon, float2[] vertexValues)
        {
            int n = polygon.Length;
            var toVertex = new double2[n];
            var distance = new double[n];
            for (int i = 0; i < n; i++)
            {
                toVertex[i] = polygon[i] - p;
                distance[i] = math.length(toVertex[i]);
            }

            // tanHalf[i] = tan(alpha_i / 2), where alpha_i is the angle at p spanned by edge i -> i+1.
            var tanHalf = new double[n];
            for (int i = 0; i < n; i++)
            {
                int next = (i + 1) % n;
                double cross = toVertex[i].x * toVertex[next].y - toVertex[i].y * toVertex[next].x;
                double dot = math.dot(toVertex[i], toVertex[next]);
                tanHalf[i] = cross / (distance[i] * distance[next] + dot);
            }

            double2 weightedSum = double2.zero;
            double weightSum = 0.0;
            for (int i = 0; i < n; i++)
            {
                int previous = (i + n - 1) % n;
                double weight = (tanHalf[previous] + tanHalf[i]) / distance[i];
                weightedSum += weight * (double2)vertexValues[i];
                weightSum += weight;
            }

            double2 result = weightedSum / weightSum;
            return new float2((float)result.x, (float)result.y);
        }
    }
}
