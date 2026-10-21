using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>
    /// Answers "what is the closest point on the original boundary, and what is the tangent there?".
    /// The tangent is interpolated linearly along the closest edge between the tangents of its two end vertices,
    /// so it equals the exact vertex tangent at the vertices and is continuous along the whole boundary.
    /// </summary>
    public sealed class BoundaryLookup
    {
        readonly float2[] ring;
        readonly float2[] tangents;   // unit tangent at every ring vertex

        public BoundaryLookup(float2[] ring, float2[] vertexTangents)
        {
            this.ring = ring;
            tangents = vertexTangents;
        }

        /// <summary>Distance from p to the boundary; also returns the unit tangent at the closest boundary point.</summary>
        public float Closest(float2 p, out float2 tangent)
        {
            float bestDistance = float.MaxValue;
            int bestEdge = 0;
            float bestParameter = 0f;

            for (int i = 0; i < ring.Length; i++)
            {
                float2 a = ring[i];
                float2 ab = ring[(i + 1) % ring.Length] - a;
                float t = math.clamp(math.dot(p - a, ab) / math.dot(ab, ab), 0f, 1f);
                float distance = math.distance(p, a + ab * t);
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    bestEdge = i;
                    bestParameter = t;
                }
            }

            float2 blended = math.lerp(tangents[bestEdge], tangents[(bestEdge + 1) % ring.Length], bestParameter);
            tangent = math.normalizesafe(blended, tangents[bestEdge]);
            return bestDistance;
        }
    }
}
