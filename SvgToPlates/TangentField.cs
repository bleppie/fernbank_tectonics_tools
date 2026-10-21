using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>Tangent vectors at the vertices of the original boundary.</summary>
    public static class TangentField
    {
        /// <summary>
        /// Tangent at every vertex of a counter-clockwise ring: the weighted average of the unit tangents of the
        /// incoming and outgoing edge. Weights are the edge lengths, which makes the result the normalised sum of the
        /// two edge vectors. Averaging unit vectors (instead of angles) avoids trouble at the +-pi wrap-around.
        /// </summary>
        public static float2[] ComputeBoundaryTangents(float2[] ring)
        {
            int n = ring.Length;
            var result = new float2[n];
            for (int i = 0; i < n; i++)
            {
                float2 incomingEdge = ring[i] - ring[(i + n - 1) % n];
                float2 outgoingEdge = ring[(i + 1) % n] - ring[i];

                float incomingLength = math.length(incomingEdge);
                float outgoingLength = math.length(outgoingEdge);
                float2 incomingTangent = incomingEdge / incomingLength;   // = (cos, sin) of atan2(dz, dx)
                float2 outgoingTangent = outgoingEdge / outgoingLength;

                float incomingWeight = incomingLength;    // use 1f for both to get a plain average
                float outgoingWeight = outgoingLength;
                float2 blended = incomingWeight * incomingTangent + outgoingWeight * outgoingTangent;
                result[i] = math.normalizesafe(blended, outgoingTangent);
            }
            return result;
        }

        /// <summary>
        /// Tangent at every vertex of an open path, in travel direction. Interior vertices use the same edge-length
        /// weighted average as the plate outlines; the two end vertices use the direction of their only segment.
        /// </summary>
        public static float2[] ComputePathTangents(float2[] path)
        {
            int n = path.Length;
            var result = new float2[n];
            for (int i = 0; i < n; i++)
            {
                float2 incomingEdge = i > 0 ? path[i] - path[i - 1] : path[i + 1] - path[i];
                float2 outgoingEdge = i < n - 1 ? path[i + 1] - path[i] : path[i] - path[i - 1];

                float2 incomingTangent = math.normalize(incomingEdge);
                float2 outgoingTangent = math.normalize(outgoingEdge);
                float2 blended = math.length(incomingEdge) * incomingTangent + math.length(outgoingEdge) * outgoingTangent;
                result[i] = math.normalizesafe(blended, outgoingTangent);
            }
            return result;
        }
    }
}
