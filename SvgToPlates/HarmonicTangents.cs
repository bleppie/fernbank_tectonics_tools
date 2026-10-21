using System;
using System.Collections.Generic;
using Unity.Mathematics;

namespace PlateMeshes
{
    /// <summary>
    /// Fills in the tangent vector at the free vertices of a triangle mesh so that it varies as smoothly as possible
    /// (minimum Dirichlet energy) given fixed values at the boundary vertices.
    ///
    /// This is linear (barycentric) interpolation inside every triangle, with the unknown vertex values chosen by
    /// solving the discrete Laplace equation: every free vertex ends up as a weighted average of its neighbours.
    /// Weights are cotangent weights, clamped to be positive, which guarantees that no value overshoots the
    /// boundary values. Free vertices on the outer offset ring have no fixed neighbours outside, so the field simply
    /// continues smoothly outward from the original boundary.
    /// </summary>
    public static class HarmonicTangents
    {
        const double MinimumWeight = 1e-3;     // keeps every mesh edge connected
        const double MaximumCotangent = 50.0;  // tames sliver triangles
        const int MaximumIterations = 5000;
        const double RelativeTolerance = 1e-10;

        struct Link
        {
            public int to;
            public double weight;
        }

        /// <summary>
        /// Returns a tangent vector for every vertex. Vertices flagged in isFixed keep fixedValues[i] (unit length);
        /// the others get the solved vector, which is NOT normalised: it is an average of unit vectors, so its length
        /// is at most 1 and shrinks where the boundary tangents point in opposing directions (the middle of the plate).
        /// </summary>
        public static float2[] Solve(IReadOnlyList<float2> positions, IReadOnlyList<int> triangles, bool[] isFixed, float2[] fixedValues)
        {
            int vertexCount = positions.Count;
            List<Link>[] graph = BuildCotangentGraph(positions, triangles);

            var freeIndex = new int[vertexCount];   // position of each free vertex among the unknowns, -1 if fixed
            int freeCount = 0;
            for (int v = 0; v < vertexCount; v++)
                freeIndex[v] = isFixed[v] ? -1 : freeCount++;

            var result = (float2[])fixedValues.Clone();
            if (freeCount == 0) return result;

            var fixedX = new double[vertexCount];
            var fixedZ = new double[vertexCount];
            for (int v = 0; v < vertexCount; v++)
            {
                fixedX[v] = fixedValues[v].x;
                fixedZ[v] = fixedValues[v].y;
            }

            double[] solutionX = SolveComponent(graph, isFixed, freeIndex, freeCount, fixedX);
            double[] solutionZ = SolveComponent(graph, isFixed, freeIndex, freeCount, fixedZ);

            for (int v = 0; v < vertexCount; v++)
            {
                if (freeIndex[v] < 0) continue;
                result[v] = new float2((float)solutionX[freeIndex[v]], (float)solutionZ[freeIndex[v]]);
            }
            return result;
        }

        /// <summary>
        /// Scalar version: free vertices get the smoothest (harmonic) extension of the values at the fixed vertices.
        /// fixedValues must hold a value for every vertex; only those flagged in isFixed are kept as they are.
        /// </summary>
        public static float[] SolveScalar(IReadOnlyList<float2> positions, IReadOnlyList<int> triangles, bool[] isFixed, float[] fixedValues)
        {
            int vertexCount = positions.Count;
            List<Link>[] graph = BuildCotangentGraph(positions, triangles);

            var freeIndex = new int[vertexCount];
            int freeCount = 0;
            for (int v = 0; v < vertexCount; v++)
                freeIndex[v] = isFixed[v] ? -1 : freeCount++;

            var result = (float[])fixedValues.Clone();
            if (freeCount == 0) return result;

            var fixedDouble = new double[vertexCount];
            for (int v = 0; v < vertexCount; v++) fixedDouble[v] = fixedValues[v];

            double[] solution = SolveComponent(graph, isFixed, freeIndex, freeCount, fixedDouble);
            for (int v = 0; v < vertexCount; v++)
                if (freeIndex[v] >= 0) result[v] = (float)solution[freeIndex[v]];
            return result;
        }

        // ================================================================ mesh graph

        static List<Link>[] BuildCotangentGraph(IReadOnlyList<float2> positions, IReadOnlyList<int> triangles)
        {
            int n = positions.Count;
            var edgeWeights = new Dictionary<long, double>();
            for (int t = 0; t < triangles.Count; t += 3)
            {
                int a = triangles[t], b = triangles[t + 1], c = triangles[t + 2];
                AddCotangentWeight(edgeWeights, n, positions, a, b, c);   // angle at c, opposite edge a-b
                AddCotangentWeight(edgeWeights, n, positions, b, c, a);
                AddCotangentWeight(edgeWeights, n, positions, c, a, b);
            }

            var graph = new List<Link>[n];
            for (int i = 0; i < n; i++) graph[i] = new List<Link>();
            foreach (KeyValuePair<long, double> entry in edgeWeights)
            {
                int i = (int)(entry.Key / n);
                int j = (int)(entry.Key % n);
                double weight = Math.Max(entry.Value, MinimumWeight);
                graph[i].Add(new Link { to = j, weight = weight });
                graph[j].Add(new Link { to = i, weight = weight });
            }
            return graph;
        }

        /// <summary>Adds 0.5 * cot(angle at apex) to the weight of edge (i, j). Interior edges receive two such terms.</summary>
        static void AddCotangentWeight(Dictionary<long, double> weights, int n, IReadOnlyList<float2> positions, int i, int j, int apex)
        {
            double2 u = positions[i] - positions[apex];
            double2 v = positions[j] - positions[apex];
            double cross = Math.Abs(u.x * v.y - u.y * v.x);
            double cotangent = Math.Min(math.dot(u, v) / Math.Max(cross, 1e-12), MaximumCotangent);

            long key = i < j ? (long)i * n + j : (long)j * n + i;
            weights.TryGetValue(key, out double sum);
            weights[key] = sum + 0.5 * cotangent;
        }

        // ================================================================ linear solve (preconditioned conjugate gradients)

        /// <summary>
        /// Solves A x = b for one vector component, where A is the graph Laplacian restricted to the free vertices
        /// and b collects the contribution of the fixed neighbours. A is symmetric positive definite.
        /// </summary>
        static double[] SolveComponent(List<Link>[] graph, bool[] isFixed, int[] freeIndex, int freeCount, double[] fixedValue)
        {
            var diagonal = new double[freeCount];
            var b = new double[freeCount];
            for (int v = 0; v < graph.Length; v++)
            {
                int row = freeIndex[v];
                if (row < 0) continue;
                foreach (Link link in graph[v])
                {
                    diagonal[row] += link.weight;
                    if (isFixed[link.to]) b[row] += link.weight * fixedValue[link.to];
                }
            }

            var x = new double[freeCount];
            var r = (double[])b.Clone();              // residual b - A x, with x = 0
            var z = new double[freeCount];            // preconditioned residual (Jacobi)
            for (int i = 0; i < freeCount; i++) z[i] = r[i] / diagonal[i];
            var direction = (double[])z.Clone();
            var product = new double[freeCount];

            double rz = Dot(r, z);
            double stopAt = RelativeTolerance * Math.Sqrt(Dot(b, b));
            for (int iteration = 0; iteration < MaximumIterations; iteration++)
            {
                if (Math.Sqrt(Dot(r, r)) <= stopAt) break;

                MultiplyLaplacian(graph, isFixed, freeIndex, diagonal, direction, product);
                double alpha = rz / Dot(direction, product);
                for (int i = 0; i < freeCount; i++)
                {
                    x[i] += alpha * direction[i];
                    r[i] -= alpha * product[i];
                    z[i] = r[i] / diagonal[i];
                }

                double rzNext = Dot(r, z);
                double beta = rzNext / rz;
                rz = rzNext;
                for (int i = 0; i < freeCount; i++) direction[i] = z[i] + beta * direction[i];
            }
            return x;
        }

        static void MultiplyLaplacian(List<Link>[] graph, bool[] isFixed, int[] freeIndex, double[] diagonal, double[] vector, double[] result)
        {
            for (int v = 0; v < graph.Length; v++)
            {
                int row = freeIndex[v];
                if (row < 0) continue;
                double sum = diagonal[row] * vector[row];
                foreach (Link link in graph[v])
                    if (!isFixed[link.to]) sum -= link.weight * vector[freeIndex[link.to]];
                result[row] = sum;
            }
        }

        static double Dot(double[] a, double[] b)
        {
            double sum = 0.0;
            for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
            return sum;
        }
    }
}
