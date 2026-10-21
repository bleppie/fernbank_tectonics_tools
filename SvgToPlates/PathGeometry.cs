using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>Result of projecting a point onto an open polyline.</summary>
    public struct PathProjection
    {
        public float distance;    // distance to the closest point on the polyline
        public float arcLength;   // distance along the polyline from its start to that closest point
        public float side;        // +1 if the point is to the left of the travel direction, -1 if to the right
        public int segment;       // index of the segment holding the closest point (segment i runs from point i to i + 1)
        public float parameter;   // 0..1 position of the closest point along that segment
    }

    /// <summary>
    /// Helpers for open polylines (paths). "Left" and "right" refer to the travel direction from the first to the last
    /// point, seen from above with x to the right and z up.
    /// </summary>
    public static class PathGeometry
    {
        /// <summary>Cumulative length along the path: element i is the distance from the start to point i.</summary>
        public static float[] CumulativeLength(float2[] path)
        {
            var result = new float[path.Length];
            for (int i = 1; i < path.Length; i++)
                result[i] = result[i - 1] + math.distance(path[i], path[i - 1]);
            return result;
        }

        static float2 LeftNormal(float2 direction) =>
            math.normalizesafe(new float2(-direction.y, direction.x));

        /// <summary>
        /// Offsets an open polyline sideways: positive distance moves it to the left, negative to the right.
        /// Interior vertices move along the mitred (averaged) normal, limited to maxMiterScale * |distance|.
        /// The two end vertices move perpendicular to their only segment, which gives a square cut at each end.
        /// </summary>
        public static float2[] OffsetOpenPolyline(float2[] path, float distance, float maxMiterScale)
        {
            int n = path.Length;
            var result = new float2[n];
            for (int i = 0; i < n; i++)
            {
                float2 normalIn = i > 0 ? LeftNormal(path[i] - path[i - 1]) : float2.zero;
                float2 normalOut = i < n - 1 ? LeftNormal(path[i + 1] - path[i]) : float2.zero;
                if (i == 0) normalIn = normalOut;
                if (i == n - 1) normalOut = normalIn;

                float2 vertexNormal = math.normalizesafe(normalIn + normalOut, normalOut);
                float cosHalfTurn = math.max(math.dot(vertexNormal, normalOut), 1f / maxMiterScale);
                result[i] = path[i] + vertexNormal * (distance / cosHalfTurn);
            }
            return result;
        }

        /// <summary>
        /// Closest point on the polyline to p: its distance, arc length position, and which side of the path p is on.
        /// Where the closest point is a vertex, the side is judged against the average direction of the two
        /// neighbouring segments, so it stays consistent around bends.
        /// </summary>
        public static PathProjection Project(float2 p, float2[] path, float[] cumulativeLength)
        {
            float bestDistance = float.MaxValue;
            int bestSegment = 0;
            float bestParameter = 0f;

            for (int i = 0; i < path.Length - 1; i++)
            {
                float2 a = path[i];
                float2 ab = path[i + 1] - a;
                float t = math.clamp(math.dot(p - a, ab) / math.dot(ab, ab), 0f, 1f);
                float distance = math.distance(p, a + ab * t);
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    bestSegment = i;
                    bestParameter = t;
                }
            }

            float2 segment = path[bestSegment + 1] - path[bestSegment];
            float2 tangent = math.normalize(segment);
            if (bestParameter == 0f && bestSegment > 0)
                tangent = math.normalizesafe(tangent + math.normalize(path[bestSegment] - path[bestSegment - 1]), tangent);
            else if (bestParameter == 1f && bestSegment < path.Length - 2)
                tangent = math.normalizesafe(tangent + math.normalize(path[bestSegment + 2] - path[bestSegment + 1]), tangent);

            float2 closest = path[bestSegment] + segment * bestParameter;
            return new PathProjection
            {
                distance = bestDistance,
                arcLength = cumulativeLength[bestSegment] + bestParameter * math.length(segment),
                side = PlateGeometry.Cross(tangent, p - closest) >= 0f ? 1f : -1f,
                segment = bestSegment,
                parameter = bestParameter
            };
        }

        /// <summary>
        /// Unit tangent at the closest point of a projection: the vertex tangents of the segment's two end points,
        /// interpolated linearly and normalised. Equals the exact vertex tangent at the vertices.
        /// </summary>
        public static float2 TangentAt(PathProjection projection, float2[] vertexTangents)
        {
            float2 blended = math.lerp(vertexTangents[projection.segment], vertexTangents[projection.segment + 1], projection.parameter);
            return math.normalizesafe(blended, vertexTangents[projection.segment]);
        }
    }
}
