using System;
using System.Collections.Generic;
using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>Oriented rectangle: two perpendicular unit axes, a centre, and the side lengths along those axes.</summary>
    public struct OrientedBox
    {
        public float2 axisX;
        public float2 axisZ;
        public float2 center;
        public float2 size;
    }

    /// <summary>
    /// Plain 2D geometry helpers. Points are float2 where .x is world x and .y is world z.
    /// A "ring" is a closed polygon outline stored without repeating its first point.
    /// Rings are kept counter-clockwise (positive signed area), so the interior is to the left of each edge.
    /// </summary>
    public static class PlateGeometry
    {
        /// <summary>Points closer than this (world units) to a ring count as lying on it.</summary>
        public const float OnBoundaryTolerance = 1e-5f;

        public static float Cross(float2 a, float2 b) => a.x * b.y - a.y * b.x;
        public static float2 Perp(float2 v) => new float2(-v.y, v.x);

        // ================================================================ ring clean-up

        public static float SignedArea(float2[] ring)
        {
            float twiceArea = 0f;
            for (int i = 0; i < ring.Length; i++)
                twiceArea += Cross(ring[i], ring[(i + 1) % ring.Length]);
            return 0.5f * twiceArea;
        }

        public static float2[] EnsureCounterClockwise(float2[] ring)
        {
            if (SignedArea(ring) >= 0f) return ring;
            float2[] reversed = (float2[])ring.Clone();
            Array.Reverse(reversed);
            return reversed;
        }

        /// <summary>Drops points that are within minDistance of the previously kept point (and of the first point at the end).</summary>
        public static float2[] RemoveClosePoints(float2[] ring, float minDistance)
        {
            var kept = new List<float2>(ring.Length);
            foreach (float2 p in ring)
            {
                if (kept.Count == 0 || math.distance(p, kept[kept.Count - 1]) > minDistance)
                    kept.Add(p);
            }
            while (kept.Count > 1 && math.distance(kept[0], kept[kept.Count - 1]) <= minDistance)
                kept.RemoveAt(kept.Count - 1);
            return kept.ToArray();
        }

        // ================================================================ convex hull + rotating calipers

        /// <summary>Andrew's monotone chain. Returns the hull counter-clockwise without collinear points.</summary>
        public static float2[] ConvexHull(IEnumerable<float2> input)
        {
            var points = new List<float2>(input);
            points.Sort((a, b) => a.x != b.x ? a.x.CompareTo(b.x) : a.y.CompareTo(b.y));

            var hull = new List<float2>();
            for (int pass = 0; pass < 2; pass++)   // pass 0 builds the lower hull, pass 1 the upper hull
            {
                int start = hull.Count;
                for (int i = 0; i < points.Count; i++)
                {
                    float2 p = points[pass == 0 ? i : points.Count - 1 - i];
                    while (hull.Count >= start + 2 &&
                           Cross(hull[hull.Count - 1] - hull[hull.Count - 2], p - hull[hull.Count - 2]) <= 0f)
                        hull.RemoveAt(hull.Count - 1);
                    hull.Add(p);
                }
                hull.RemoveAt(hull.Count - 1);     // last point of each pass is the first point of the next
            }
            return hull.ToArray();
        }

        /// <summary>
        /// Smallest-area enclosing rectangle of a convex polygon (counter-clockwise) using rotating calipers.
        /// One side of the optimal rectangle always lies on a hull edge, so each hull edge is tried in turn while
        /// three calipers (right, top, left) only ever move forward around the hull.
        /// </summary>
        public static OrientedBox FindMinimumAreaBox(float2[] hull)
        {
            int n = hull.Length;
            if (n < 3) throw new ArgumentException("Need at least three hull points.");

            var best = new OrientedBox();
            float bestArea = float.MaxValue;

            // Start the calipers at the correct places for edge 0 (brute force), then walk them.
            float2 u = EdgeDirection(hull, 0);
            float2 v = Perp(u);
            int right = IndexOfExtreme(hull, u, +1f);
            int top = IndexOfExtreme(hull, v, +1f);
            int left = IndexOfExtreme(hull, u, -1f);

            for (int edge = 0; edge < n; edge++)
            {
                u = EdgeDirection(hull, edge);   // direction along the box side that lies on this hull edge
                v = Perp(u);                     // points into the hull

                right = AdvanceWhileImproving(hull, right, u, +1f);
                top = AdvanceWhileImproving(hull, top, v, +1f);
                left = AdvanceWhileImproving(hull, left, u, -1f);

                float minU = math.dot(hull[left], u);
                float maxU = math.dot(hull[right], u);
                float minV = math.dot(hull[edge], v);
                float maxV = math.dot(hull[top], v);

                float width = maxU - minU;
                float height = maxV - minV;
                float area = width * height;
                if (area < bestArea)
                {
                    bestArea = area;
                    best.axisX = u;
                    best.axisZ = v;
                    best.size = new float2(width, height);
                    best.center = u * (0.5f * (minU + maxU)) + v * (0.5f * (minV + maxV));
                }
            }
            return best;
        }

        /// <summary>Rotates the box frame by 90 degrees if needed so that size.x >= size.y.</summary>
        public static OrientedBox WithLongSideAlongX(OrientedBox box)
        {
            if (box.size.y <= box.size.x) return box;
            return new OrientedBox
            {
                axisX = box.axisZ,
                axisZ = -box.axisX,
                center = box.center,
                size = box.size.yx
            };
        }

        /// <summary>Expresses a point in the box frame: origin at the box centre, x along axisX, z along axisZ.</summary>
        public static float2 ToBoxSpace(float2 p, OrientedBox box)
        {
            float2 d = p - box.center;
            return new float2(math.dot(d, box.axisX), math.dot(d, box.axisZ));
        }

        static float2 EdgeDirection(float2[] hull, int edge) =>
            math.normalize(hull[(edge + 1) % hull.Length] - hull[edge]);

        static int IndexOfExtreme(float2[] points, float2 axis, float sign)
        {
            int best = 0;
            for (int i = 1; i < points.Length; i++)
                if (sign * math.dot(points[i], axis) > sign * math.dot(points[best], axis)) best = i;
            return best;
        }

        static int AdvanceWhileImproving(float2[] hull, int index, float2 axis, float sign)
        {
            int n = hull.Length;
            while (sign * math.dot(hull[(index + 1) % n], axis) > sign * math.dot(hull[index], axis))
                index = (index + 1) % n;
            return index;
        }

        // ================================================================ offsetting

        /// <summary>
        /// Moves every vertex of a counter-clockwise ring outward along its averaged (mitred) normal so that
        /// the adjacent edges move by 'distance'. The mitre is limited to maxMiterScale * distance.
        /// Short edges at concave corners can flip and create small loops; run RemoveSelfIntersectionLoops afterwards.
        /// </summary>
        public static float2[] OffsetOutward(float2[] ring, float distance, float maxMiterScale)
        {
            int n = ring.Length;
            var result = new float2[n];
            for (int i = 0; i < n; i++)
            {
                float2 previous = ring[(i + n - 1) % n];
                float2 next = ring[(i + 1) % n];
                float2 incomingNormal = OutwardNormal(ring[i] - previous);
                float2 outgoingNormal = OutwardNormal(next - ring[i]);

                float2 vertexNormal = math.normalizesafe(incomingNormal + outgoingNormal, outgoingNormal);
                float cosHalfTurn = math.max(math.dot(vertexNormal, outgoingNormal), 1f / maxMiterScale);
                result[i] = ring[i] + vertexNormal * (distance / cosHalfTurn);
            }
            return result;
        }

        static float2 OutwardNormal(float2 edge) => math.normalizesafe(new float2(edge.y, -edge.x));

        /// <summary>
        /// Repeatedly finds a pair of crossing edges. The crossing splits the ring into two closed loops; the loop
        /// that is turned inside out (clockwise, negative area) is the artefact of offsetting and is discarded.
        /// If both loops have the same orientation, the smaller one is discarded.
        /// </summary>
        public static float2[] RemoveSelfIntersectionLoops(float2[] ring)
        {
            var points = new List<float2>(ring);
            for (int iteration = 0; iteration < 10000; iteration++)
            {
                if (!TryFindSelfIntersection(points, out int i, out int j, out float2 crossing))
                    return points.ToArray();

                // Edges i and j cross at 'crossing'. Loop A is the crossing followed by vertices i+1 .. j,
                // loop B is the crossing followed by vertices j+1 .. end, 0 .. i.
                var loopA = new List<float2> { crossing };
                loopA.AddRange(points.GetRange(i + 1, j - i));

                var loopB = new List<float2> { crossing };
                loopB.AddRange(points.GetRange(j + 1, points.Count - j - 1));
                loopB.AddRange(points.GetRange(0, i + 1));

                float areaA = SignedArea(loopA.ToArray());
                float areaB = SignedArea(loopB.ToArray());
                bool keepA = (areaA >= 0f) != (areaB >= 0f) ? areaA >= 0f : areaA > areaB;
                points = keepA ? loopA : loopB;
            }
            throw new InvalidOperationException("Could not remove self-intersections from the offset outline.");
        }

        static bool TryFindSelfIntersection(List<float2> ring, out int edgeI, out int edgeJ, out float2 crossing)
        {
            int n = ring.Count;
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 2; j < n; j++)
                {
                    if (i == 0 && j == n - 1) continue;   // these two edges are neighbours
                    if (TryIntersectSegments(ring[i], ring[(i + 1) % n], ring[j], ring[(j + 1) % n], out crossing))
                    {
                        edgeI = i;
                        edgeJ = j;
                        return true;
                    }
                }
            }
            edgeI = edgeJ = -1;
            crossing = default;
            return false;
        }

        // ================================================================ segment tests

        /// <summary>True if segments ab and cd cross strictly inside both (shared end points do not count).</summary>
        public static bool TryIntersectSegments(float2 a, float2 b, float2 c, float2 d, out float2 point)
        {
            point = default;
            float2 r = b - a;
            float2 s = d - c;
            float denominator = Cross(r, s);
            if (math.abs(denominator) < 1e-12f) return false;   // parallel

            float2 ac = c - a;
            float t = Cross(ac, s) / denominator;
            float u = Cross(ac, r) / denominator;
            if (t <= 0f || t >= 1f || u <= 0f || u >= 1f) return false;

            point = a + r * t;
            return true;
        }

        /// <summary>Checks the constraint edges (pairs of indices into points) for crossings.</summary>
        public static bool HasCrossingEdges(List<float2> points, List<int> edgePairs)
        {
            int edgeCount = edgePairs.Count / 2;
            for (int e = 0; e < edgeCount; e++)
            {
                for (int f = e + 1; f < edgeCount; f++)
                {
                    int a = edgePairs[2 * e], b = edgePairs[2 * e + 1];
                    int c = edgePairs[2 * f], d = edgePairs[2 * f + 1];
                    if (a == c || a == d || b == c || b == d) continue;   // share an end point
                    if (TryIntersectSegments(points[a], points[b], points[c], points[d], out _)) return true;
                }
            }
            return false;
        }

        // ================================================================ distance to a ring

        public static float DistanceToRing(float2 p, float2[] ring)
        {
            float best = float.MaxValue;
            for (int i = 0; i < ring.Length; i++)
                best = math.min(best, DistanceToSegment(p, ring[i], ring[(i + 1) % ring.Length]));
            return best;
        }

        /// <summary>Distance to the ring: negative inside, positive outside, exactly 0 on the ring.</summary>
        public static float SignedDistanceToRing(float2 p, float2[] ring)
        {
            float distance = DistanceToRing(p, ring);
            if (distance < OnBoundaryTolerance) return 0f;
            return IsInside(p, ring) ? -distance : distance;
        }

        /// <summary>Crossing-number point-in-polygon test.</summary>
        public static bool IsInside(float2 p, float2[] ring)
        {
            bool inside = false;
            for (int i = 0, j = ring.Length - 1; i < ring.Length; j = i++)
            {
                float2 a = ring[i];
                float2 b = ring[j];
                if ((a.y > p.y) != (b.y > p.y))
                {
                    float xAtP = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
                    if (p.x < xAtP) inside = !inside;
                }
            }
            return inside;
        }

        static float DistanceToSegment(float2 p, float2 a, float2 b)
        {
            float2 ab = b - a;
            float t = math.clamp(math.dot(p - a, ab) / math.dot(ab, ab), 0f, 1f);
            return math.distance(p, a + ab * t);
        }
    }
}
