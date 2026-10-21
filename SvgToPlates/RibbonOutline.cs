using System;
using System.Collections.Generic;
using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>
    /// Outline of a ribbon: the region within a given distance of an open path, with square cuts at the two ends.
    ///
    /// Offsetting each vertex along its own normal breaks down at sharp corners when a neighbouring segment is shorter
    /// than the offset: the offset edges then cross the path itself. This class builds the outline robustly instead:
    ///   1. Every segment is shifted sideways by the extent (a strip). At each corner the shifted segments are joined
    ///      by a round arc, a bevel or a mitre on the outside, and on the inside they are joined through the path vertex
    ///      itself, so that every strip is fully enclosed by the curve. The resulting closed curve may cross itself.
    ///   2. All crossings (and places where a vertex touches an edge) are found, which turns the curve into a planar graph.
    ///   3. The outer boundary of that graph is traced. It is the outline of the union of all strips and joins,
    ///      correct at any corner angle and any extent, and it never crosses the path.
    ///
    /// All computation is done in double precision.
    /// </summary>
    public static class RibbonOutline
    {
        const double Tolerance = 1e-5;   // world units: points closer than this are the same point; a vertex this close to an edge touches it

        /// <summary>Counter-clockwise outline around the path. The path's end vertices lie on the end cuts if those are on the outline.</summary>
        public static float2[] Build(float2[] path, float extent, RibbonJoin join, float curveTolerance, float maxMiterScale)
        {
            var points = new double2[path.Length];
            for (int i = 0; i < points.Length; i++) points[i] = path[i];

            List<double2> raw = BuildRawOutline(points, extent, join, curveTolerance, maxMiterScale);
            List<double2> contour = TraceOuterContour(raw);

            var result = new float2[contour.Count];
            for (int i = 0; i < result.Length; i++)
                result[i] = new float2((float)contour[i].x, (float)contour[i].y);
            return result;
        }

        // ================================================================ step 1: the raw, self-crossing curve

        /// <summary>Right side forward, the end cut through the last vertex, left side backward, the start cut through the first vertex.</summary>
        static List<double2> BuildRawOutline(double2[] path, double extent, RibbonJoin join, double curveTolerance, double miterLimit)
        {
            int last = path.Length - 1;
            List<double2> right = BuildSide(path, -1.0, extent, join, curveTolerance, miterLimit);
            List<double2> left = BuildSide(path, +1.0, extent, join, curveTolerance, miterLimit);
            left.Reverse();

            var raw = new List<double2>();
            foreach (double2 p in right) AddDistinct(raw, p);
            AddDistinct(raw, path[last]);
            foreach (double2 p in left) AddDistinct(raw, p);
            AddDistinct(raw, path[0]);
            while (raw.Count > 1 && math.distance(raw[0], raw[raw.Count - 1]) <= 1e-9)
                raw.RemoveAt(raw.Count - 1);
            return raw;
        }

        static void AddDistinct(List<double2> list, double2 p)
        {
            if (list.Count == 0 || math.distance(p, list[list.Count - 1]) > 1e-9)
                list.Add(p);
        }

        /// <summary>
        /// One side of the shifted strips, in path direction. side = +1 is the left of the travel direction, -1 the right.
        /// </summary>
        static List<double2> BuildSide(double2[] path, double side, double extent, RibbonJoin join, double curveTolerance, double miterLimit)
        {
            int segments = path.Length - 1;
            var direction = new double2[segments];
            var normal = new double2[segments];
            for (int i = 0; i < segments; i++)
            {
                direction[i] = math.normalize(path[i + 1] - path[i]);
                normal[i] = side * new double2(-direction[i].y, direction[i].x);
            }

            var points = new List<double2>();
            for (int i = 0; i < segments; i++)
            {
                if (i > 0)
                    AddCornerJoin(points, path[i], direction[i - 1], direction[i], normal[i - 1], normal[i], side, extent, join, curveTolerance, miterLimit);
                points.Add(path[i] + extent * normal[i]);
                points.Add(path[i + 1] + extent * normal[i]);
            }
            return points;
        }

        static void AddCornerJoin(List<double2> points, double2 corner, double2 directionBefore, double2 directionAfter,
                                  double2 normalBefore, double2 normalAfter, double side, double extent,
                                  RibbonJoin join, double curveTolerance, double miterLimit)
        {
            double turn = Cross(directionBefore, directionAfter);   // > 0: the path turns left
            if (Math.Abs(turn) < 1e-12) return;                     // straight on: the shifted segments already meet

            bool isOutsideOfCorner = turn * side < 0.0;
            if (!isOutsideOfCorner)
            {
                // Inside of the corner: the shifted segments overlap. Go through the path vertex so both strips stay enclosed.
                points.Add(corner);
                return;
            }

            double cosAngle = math.dot(normalBefore, normalAfter);
            switch (join)
            {
                case RibbonJoin.Round:
                {
                    double angle = Math.Atan2(Cross(normalBefore, normalAfter), cosAngle);   // signed rotation from one normal to the other
                    double step = extent > curveTolerance ? 2.0 * Math.Acos(Math.Max(-1.0, Math.Min(1.0, 1.0 - curveTolerance / extent))) : Math.PI;
                    int count = (int)Math.Ceiling(Math.Abs(angle) / Math.Max(step, 1e-3));
                    for (int k = 1; k < count; k++)
                        points.Add(corner + extent * Rotate(normalBefore, angle * k / count));
                    break;
                }
                case RibbonJoin.Miter:
                {
                    double cosHalfAngle = Math.Sqrt(Math.Max(0.0, (1.0 + cosAngle) * 0.5));
                    if (cosHalfAngle > 1e-9 && 1.0 / cosHalfAngle <= miterLimit)
                        points.Add(corner + extent * (normalBefore + normalAfter) / (1.0 + cosAngle));
                    break;   // beyond the limit it stays a bevel
                }
                // RibbonJoin.Bevel: the two shifted end points are simply connected.
            }
        }

        static double2 Rotate(double2 v, double angle)
        {
            double c = Math.Cos(angle), s = Math.Sin(angle);
            return new double2(c * v.x - s * v.y, s * v.x + c * v.y);
        }

        static double Cross(double2 a, double2 b) => a.x * b.y - a.y * b.x;

        // ================================================================ steps 2 and 3: planar graph and its outer boundary

        static List<double2> TraceOuterContour(List<double2> raw)
        {
            int m = raw.Count;

            // One graph node per distinct position; raw vertices at the same position share a node.
            var positions = new List<double2>();
            var nodeOfVertex = new int[m];
            for (int v = 0; v < m; v++)
                nodeOfVertex[v] = NodeFor(positions, raw[v]);

            // cuts[e]: extra nodes lying on edge e (edge e runs from raw[e] to raw[e + 1]), with their position along it.
            var cuts = new List<(double t, int node)>[m];
            for (int e = 0; e < m; e++) cuts[e] = new List<(double t, int node)>();
            AddTouchingVertexCuts(raw, nodeOfVertex, cuts);
            AddCrossingCuts(raw, positions, cuts);

            HashSet<int>[] adjacency = BuildAdjacency(m, nodeOfVertex, cuts, positions.Count);
            List<int> loop = WalkOuterFace(positions, adjacency);

            var contour = new List<double2>(loop.Count);
            foreach (int node in loop) contour.Add(positions[node]);
            RemoveSpikes(contour);
            if (contour.Count < 3) throw new InvalidOperationException("The ribbon outline collapsed.");
            return contour;
        }

        static int NodeFor(List<double2> positions, double2 p)
        {
            for (int k = 0; k < positions.Count; k++)
                if (math.distance(positions[k], p) < Tolerance) return k;
            positions.Add(p);
            return positions.Count - 1;
        }

        /// <summary>A vertex that lies on another edge (a T-junction, e.g. at exact right angles) splits that edge.</summary>
        static void AddTouchingVertexCuts(List<double2> raw, int[] nodeOfVertex, List<(double t, int node)>[] cuts)
        {
            int m = raw.Count;
            for (int e = 0; e < m; e++)
            {
                double2 a = raw[e];
                double2 b = raw[(e + 1) % m];
                double2 ab = b - a;
                double lengthSquared = math.dot(ab, ab);
                if (lengthSquared < 1e-24) continue;

                for (int v = 0; v < m; v++)
                {
                    if (v == e || v == (e + 1) % m) continue;
                    double t = math.dot(raw[v] - a, ab) / lengthSquared;
                    if (t <= 0.0 || t >= 1.0) continue;
                    if (math.distance(raw[v], a + ab * t) >= Tolerance) continue;
                    if (math.distance(raw[v], a) <= Tolerance || math.distance(raw[v], b) <= Tolerance) continue;
                    cuts[e].Add((t, nodeOfVertex[v]));
                }
            }
        }

        /// <summary>Two edges that cross strictly inside both get a new node at the crossing.</summary>
        static void AddCrossingCuts(List<double2> raw, List<double2> positions, List<(double t, int node)>[] cuts)
        {
            int m = raw.Count;
            for (int i = 0; i < m; i++)
            {
                double2 a = raw[i];
                double2 b = raw[(i + 1) % m];
                for (int j = i + 2; j < m; j++)
                {
                    if (i == 0 && j == m - 1) continue;   // neighbouring edges
                    double2 c = raw[j];
                    double2 d = raw[(j + 1) % m];

                    double2 r = b - a;
                    double2 s = d - c;
                    double denominator = Cross(r, s);
                    if (Math.Abs(denominator) < 1e-14) continue;   // parallel

                    double2 ac = c - a;
                    double t = Cross(ac, s) / denominator;
                    double u = Cross(ac, r) / denominator;
                    if (t <= 0.0 || t >= 1.0 || u <= 0.0 || u >= 1.0) continue;

                    int node = NodeFor(positions, a + r * t);
                    cuts[i].Add((t, node));
                    cuts[j].Add((u, node));
                }
            }
        }

        static HashSet<int>[] BuildAdjacency(int m, int[] nodeOfVertex, List<(double t, int node)>[] cuts, int nodeCount)
        {
            var adjacency = new HashSet<int>[nodeCount];
            for (int k = 0; k < nodeCount; k++) adjacency[k] = new HashSet<int>();

            for (int e = 0; e < m; e++)
            {
                int first = nodeOfVertex[e];
                int last = nodeOfVertex[(e + 1) % m];
                cuts[e].Sort((x, y) => x.t.CompareTo(y.t));

                int previous = first;
                foreach ((double t, int node) cut in cuts[e])
                {
                    if (cut.node == first || cut.node == last || cut.node == previous) continue;
                    Link(adjacency, previous, cut.node);
                    previous = cut.node;
                }
                Link(adjacency, previous, last);
            }
            return adjacency;
        }

        static void Link(HashSet<int>[] adjacency, int a, int b)
        {
            if (a == b) return;
            adjacency[a].Add(b);
            adjacency[b].Add(a);
        }

        /// <summary>
        /// Walks the boundary of the unbounded face, counter-clockwise. Starting at the leftmost node (which is certainly
        /// on it), at every node the next edge is the first one met when turning counter-clockwise from the edge we came by.
        /// </summary>
        static List<int> WalkOuterFace(List<double2> positions, HashSet<int>[] adjacency)
        {
            int start = -1;
            for (int k = 0; k < adjacency.Length; k++)
            {
                if (adjacency[k].Count == 0) continue;
                if (start < 0 || positions[k].x < positions[start].x ||
                    (positions[k].x == positions[start].x && positions[k].y < positions[start].y))
                    start = k;
            }
            if (start < 0) throw new InvalidOperationException("The ribbon outline is empty.");

            var loop = new List<int> { start };
            int current = start;
            double2 back = new double2(-1.0, 0.0);   // as if we arrived from the left
            int maxSteps = 20 * adjacency.Length + 100;
            for (int step = 0; step < maxSteps; step++)
            {
                double backAngle = Math.Atan2(back.y, back.x);
                int next = -1;
                double bestTurn = double.MaxValue;
                foreach (int neighbour in adjacency[current])
                {
                    double2 toNeighbour = positions[neighbour] - positions[current];
                    double turn = Math.Atan2(toNeighbour.y, toNeighbour.x) - backAngle;
                    turn -= 2.0 * Math.PI * Math.Floor(turn / (2.0 * Math.PI));   // 0 .. 2 pi
                    if (turn < 1e-12) turn = 2.0 * Math.PI;                       // the edge we came by is the last resort
                    if (turn < bestTurn)
                    {
                        bestTurn = turn;
                        next = neighbour;
                    }
                }

                back = positions[current] - positions[next];
                current = next;
                if (current == start) return loop;
                loop.Add(current);
            }
            throw new InvalidOperationException("Tracing the ribbon outline did not close.");
        }

        /// <summary>Removes out-and-back spikes (A, B, A) that can appear where the curve doubles back on itself.</summary>
        static void RemoveSpikes(List<double2> contour)
        {
            bool changed = true;
            while (changed && contour.Count > 3)
            {
                changed = false;
                for (int i = 0; i < contour.Count; i++)
                {
                    int previous = (i + contour.Count - 1) % contour.Count;
                    int next = (i + 1) % contour.Count;
                    if (math.distance(contour[previous], contour[next]) >= Tolerance) continue;

                    // Drop the tip and the repeated point, highest index first.
                    int a = Math.Max(i, next), b = Math.Min(i, next);
                    contour.RemoveAt(a);
                    contour.RemoveAt(b);
                    changed = true;
                    break;
                }
            }
        }
    }
}
