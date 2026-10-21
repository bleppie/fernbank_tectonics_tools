using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>A path along which two plates abut.</summary>
    public struct SharedPath
    {
        public int plateA;        // lower plate number
        public int plateB;        // higher plate number
        public string name;       // e.g. "ribbon_00_01"
        public float2[] points;   // open polyline, ordered clockwise around plateA (seen from above)
    }

    /// <summary>Finds the paths where plates abut by comparing every pair of plate outlines.</summary>
    public static class SharedPathFinder
    {
        /// <summary>
        /// rings are the counter-clockwise plate outlines; the list index is the plate number. An outline edge counts as
        /// shared if both end points and its midpoint lie within 'tolerance' of the other plate's outline. Each maximal
        /// run of shared edges becomes one path, taken from the lower-numbered plate and reversed so it runs clockwise
        /// around that plate (the plate is then on the right of the travel direction).
        /// </summary>
        public static List<SharedPath> Find(List<float2[]> rings, float tolerance)
        {
            var paths = new List<SharedPath>();
            for (int a = 0; a < rings.Count; a++)
                for (int b = a + 1; b < rings.Count; b++)
                    AddPathsBetween(paths, rings, a, b, tolerance);
            return paths;
        }

        static void AddPathsBetween(List<SharedPath> paths, List<float2[]> rings, int a, int b, float tolerance)
        {
            float2[] ringA = rings[a];
            float2[] ringB = rings[b];
            int count = ringA.Length;

            var isShared = new bool[count];
            int sharedCount = 0;
            for (int i = 0; i < count; i++)
            {
                float2 start = ringA[i];
                float2 end = ringA[(i + 1) % count];
                isShared[i] = IsOnRing(start, ringB, tolerance) &&
                              IsOnRing(end, ringB, tolerance) &&
                              IsOnRing(0.5f * (start + end), ringB, tolerance);
                if (isShared[i]) sharedCount++;
            }

            if (sharedCount == 0) return;
            if (sharedCount == count)
            {
                Debug.LogWarning($"plate_{a:00} and plate_{b:00} share their entire outline (a closed path). No ribbon created.");
                return;
            }

            // Start walking at an edge that is not shared, so no run wraps around the end of the array.
            int first = 0;
            while (isShared[first]) first++;

            var run = new List<int>();
            int part = 0;
            for (int k = 1; k <= count; k++)
            {
                int edge = (first + k) % count;
                if (isShared[edge])
                {
                    run.Add(edge);
                }
                else if (run.Count > 0)
                {
                    paths.Add(MakePath(ringA, run, a, b, part++));
                    run.Clear();
                }
            }
        }

        static SharedPath MakePath(float2[] ringA, List<int> edges, int a, int b, int part)
        {
            // Edges are consecutive; the polyline visits the start of each edge and then the end of the last one.
            var points = new List<float2>(edges.Count + 1);
            foreach (int edge in edges)
                points.Add(ringA[edge]);
            points.Add(ringA[(edges[edges.Count - 1] + 1) % ringA.Length]);

            // The ring runs counter-clockwise; reverse it so the path runs clockwise around plate a.
            points.Reverse();

            // The first path between two plates has no suffix; further ones get _b, _c, ...
            string suffix = part == 0 ? "" : "_" + (char)('a' + part);
            return new SharedPath
            {
                plateA = a,
                plateB = b,
                name = $"ribbon_{a:00}_{b:00}{suffix}",
                points = points.ToArray()
            };
        }

        static bool IsOnRing(float2 p, float2[] ring, float tolerance) =>
            PlateGeometry.DistanceToRing(p, ring) <= tolerance;
    }
}
