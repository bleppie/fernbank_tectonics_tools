using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>
    /// Builds one ribbon mesh per path where two plates abut.
    ///
    /// The ribbon is the region within ribbonExtent of the path, closed by a square cut at each end (see RibbonOutline). The original path is
    /// kept as constrained edges inside it, exactly as the original outline is kept inside the offset plates, and the
    /// same triangulation settings apply (interior grid, refinement, weld distance, miter limit).
    ///
    /// Three UV sets (Unity channels 0, 1 and 2, i.e. mesh.uv, mesh.uv2 and mesh.uv3), all Vector2:
    ///
    /// UV0 (mesh.uv):  texture mapping, identical to the plates': the shared bounding box of all plates maps so that its
    ///                 long side spans 0..1 and its short side is centred on 0.5, so a texture lines up across the plates
    ///                 and the ribbons.
    ///
    /// UV1 (mesh.uv2): u = distance along the path (0 at its start, path length at its end), running clockwise around the
    ///                 lower-numbered plate; v = signed distance from the path, 0 on it, positive on the side away
    ///                 from the lower-numbered plate, negative towards it.
    ///
    /// UV2 (mesh.uv3): the unit normal (x and z components), computed like the plates' normal field: the path's tangent
    ///                 field is extended smoothly across the ribbon and rotated 90 degrees, then normalised. The path
    ///                 runs clockwise around the lower-numbered plate, so the rotation is counter-clockwise and the
    ///                 normal points away from that plate (the direction in which v grows). Set to (0, 0) where the
    ///                 field is shorter than minNormalLength.
    /// </summary>
    public static class RibbonBuilder
    {
        const float PathVertexMatchTolerance = 2e-5f;   // slightly above the outline tracer's merge tolerance

        public static List<Mesh> BuildAll(List<float2[]> plateRings, PlateMeshSettings settings, float2 boxSize)
        {
            var meshes = new List<Mesh>();
            if (settings.ribbonExtent <= 2f * settings.weldDistance)
            {
                Debug.LogWarning("ribbonExtent is too small compared to weldDistance. No ribbons created.");
                return meshes;
            }

            foreach (SharedPath path in SharedPathFinder.Find(plateRings, settings.weldDistance))
            {
                Mesh mesh = BuildRibbon(path, settings, boxSize);
                if (mesh != null) meshes.Add(mesh);
            }
            return meshes;
        }

        static Mesh BuildRibbon(SharedPath shared, PlateMeshSettings settings, float2 boxSize)
        {
            float2[] path = shared.points;
            int last = path.Length - 1;
            float[] arcLength = PathGeometry.CumulativeLength(path);

            float2[] outline = RibbonOutline.Build(path, settings.ribbonExtent, settings.ribbonJoin, settings.curveTolerance, settings.maxMiterScale);
            outline = PlateGeometry.EnsureCounterClockwise(PlateGeometry.RemoveClosePoints(outline, settings.weldDistance));
            if (outline.Length < 3)
            {
                Debug.LogError($"{shared.name}: the ribbon outline collapsed. Ribbon skipped.");
                return null;
            }

            // Vertex layout: [path vertices][outline vertices that are not path vertices][interior points].
            // An outline vertex that coincides with a path vertex (for example where an end cut passes through the path's end
            // point) shares its index, so the two sets of constraint edges connect there instead of creating near-duplicate points.
            var positions = new List<float2>(path);
            var outlineIndex = new List<int>();
            foreach (float2 p in outline)
            {
                int pathVertex = IndexOfPathVertex(p, path);
                if (pathVertex >= 0)
                {
                    outlineIndex.Add(pathVertex);
                }
                else
                {
                    positions.Add(p);
                    outlineIndex.Add(positions.Count - 1);
                }
            }

            positions.AddRange(SampleInteriorPoints(path, outline, settings.interiorSpacing));

            var constraintEdges = new List<int>();
            for (int i = 0; i < last; i++)
            {
                constraintEdges.Add(i);
                constraintEdges.Add(i + 1);
            }
            for (int k = 0; k < outlineIndex.Count; k++)
            {
                constraintEdges.Add(outlineIndex[k]);
                constraintEdges.Add(outlineIndex[(k + 1) % outlineIndex.Count]);
            }

            if (PlateGeometry.HasCrossingEdges(positions, constraintEdges))
            {
                Debug.LogError($"{shared.name}: the ribbon outline crosses the path. Try a smaller ribbonExtent. Ribbon skipped.");
                return null;
            }

            PlateTriangulator.Result triangulation = PlateTriangulator.Triangulate(positions, constraintEdges, settings);
            List<int> triangles = PlateTriangulator.KeepInsideTrianglesFacingUp(triangulation.positions, triangulation.triangles, outline);
            if (triangles.Count == 0)
            {
                Debug.LogError($"{shared.name}: triangulation produced no triangles. Ribbon skipped.");
                return null;
            }

            RibbonAttributes attributes = ComputeAttributes(triangulation.positions, triangles, path, arcLength, settings);
            int[] sourceIndex = PlateMeshPipeline.CompactTriangles(triangles, triangulation.positions.Count);
            float2[] finalPositions = PlateMeshPipeline.Gather(triangulation.positions, sourceIndex);
            Vector2[] finalDistanceUvs = PlateMeshPipeline.Gather(attributes.distanceUvs, sourceIndex);
            Vector2[] finalNormalUvs = PlateMeshPipeline.Gather(attributes.normalUvs, sourceIndex);

            // Texture mapping from the shared bounding box, exactly as for the plates.
            float textureScale = 1f / math.max(boxSize.x, boxSize.y);
            var textureUvs = new Vector2[finalPositions.Length];
            for (int i = 0; i < textureUvs.Length; i++)
                textureUvs[i] = PlateMeshPipeline.TextureCoordinate(finalPositions[i], textureScale);

            Mesh mesh = CreateMesh(shared.name, finalPositions, textureUvs, finalDistanceUvs, finalNormalUvs, triangles);
            Debug.Log($"{shared.name}: path length {arcLength[last]:F3}, {finalPositions.Length} vertices, {triangles.Count / 3} triangles.");
            return mesh;
        }

        static int IndexOfPathVertex(float2 p, float2[] path)
        {
            for (int i = 0; i < path.Length; i++)
                if (math.distance(p, path[i]) < PathVertexMatchTolerance) return i;
            return -1;
        }

        /// <summary>Grid points inside the outline and at least half a cell away from both the path and the outline.</summary>
        static List<float2> SampleInteriorPoints(float2[] path, float2[] outline, float spacing)
        {
            var points = new List<float2>();
            if (spacing <= 0f) return points;

            float2 min = outline[0];
            float2 max = outline[0];
            foreach (float2 p in outline)
            {
                min = math.min(min, p);
                max = math.max(max, p);
            }

            float[] arcLength = PathGeometry.CumulativeLength(path);
            float minGap = 0.5f * spacing;
            int firstX = (int)math.ceil(min.x / spacing), lastX = (int)math.floor(max.x / spacing);
            int firstZ = (int)math.ceil(min.y / spacing), lastZ = (int)math.floor(max.y / spacing);
            for (int ix = firstX; ix <= lastX; ix++)
            {
                for (int iz = firstZ; iz <= lastZ; iz++)
                {
                    var p = new float2(ix, iz) * spacing;
                    if (!PlateGeometry.IsInside(p, outline)) continue;
                    if (PathGeometry.Project(p, path, arcLength).distance < minGap) continue;
                    if (PlateGeometry.DistanceToRing(p, outline) < minGap) continue;
                    points.Add(p);
                }
            }
            return points;
        }

        struct RibbonAttributes
        {
            public Vector2[] distanceUvs;   // UV channel 0: (u, v)
            public Vector2[] normalUvs;     // UV channel 1: unit normal (x, z)
        }

        /// <summary>
        /// v is the exact signed distance to the path. u is the arc length of the closest path point for vertices on the
        /// path (the path vertices and any points refinement adds on it); every other vertex gets the smooth (harmonic)
        /// extension of u across the ribbon, which avoids the jumps the closest-point arc length has inside bends.
        /// The normal comes from the path's tangent, extended across the ribbon the same way (see HarmonicTangents).
        /// </summary>
        static RibbonAttributes ComputeAttributes(List<float2> positions, List<int> triangles, float2[] path, float[] arcLength, PlateMeshSettings settings)
        {
            float2[] pathTangents = TangentField.ComputePathTangents(path);

            int count = positions.Count;
            var used = new bool[count];
            foreach (int index in triangles) used[index] = true;

            var v = new float[count];
            var uOnPath = new float[count];             // arc length of the closest path point, valid for every vertex
            var tangentOnPath = new float2[count];      // tangent at the closest path point, valid for every vertex
            var isFixed = new bool[count];
            for (int i = 0; i < count; i++)
            {
                PathProjection projection = PathGeometry.Project(positions[i], path, arcLength);
                bool onPath = projection.distance < PlateGeometry.OnBoundaryTolerance;
                v[i] = onPath ? 0f : projection.side * projection.distance;
                uOnPath[i] = projection.arcLength;
                tangentOnPath[i] = PathGeometry.TangentAt(projection, pathTangents);
                isFixed[i] = onPath || !used[i];
            }

            float[] u = HarmonicTangents.SolveScalar(positions, triangles, isFixed, uOnPath);
            float2[] tangentField = HarmonicTangents.Solve(positions, triangles, isFixed, tangentOnPath);

            var attributes = new RibbonAttributes { distanceUvs = new Vector2[count], normalUvs = new Vector2[count] };
            for (int i = 0; i < count; i++)
            {
                attributes.distanceUvs[i] = new Vector2(u[i], v[i]);
                attributes.normalUvs[i] = NormalOf(tangentField[i], settings.minNormalLength);
            }
            return attributes;
        }

        /// <summary>
        /// Unit normal from a tangent vector. The path runs clockwise around the lower-numbered plate, so rotating the
        /// tangent 90 degrees counter-clockwise (seen from above, x right, z up) points away from that plate, which is the
        /// direction of increasing v. Too-short vectors give (0, 0).
        /// </summary>
        static Vector2 NormalOf(float2 tangent, float minNormalLength)
        {
            float length = math.length(tangent);
            if (length <= 0f || length < minNormalLength) return Vector2.zero;
            return new Vector2(-tangent.y, tangent.x) / length;
        }

        static Mesh CreateMesh(string name, float2[] positions, Vector2[] textureUvs, Vector2[] distanceUvs, Vector2[] normalUvs, List<int> triangles)
        {
            var vertices = new Vector3[positions.Length];
            var normals = new Vector3[positions.Length];
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i] = new Vector3(positions[i].x, 0f, positions[i].y);
                normals[i] = Vector3.up;
            }

            var mesh = new Mesh { name = name };
            if (vertices.Length > 65535)
                mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
            mesh.vertices = vertices;
            mesh.normals = normals;
            mesh.SetUVs(0, new List<Vector2>(textureUvs));    // mesh.uv:  texture mapping
            mesh.SetUVs(1, new List<Vector2>(distanceUvs));   // mesh.uv2: (u, v) along / across the path
            mesh.SetUVs(2, new List<Vector2>(normalUvs));     // mesh.uv3: unit normal (x, z)
            mesh.SetTriangles(triangles, 0);
            mesh.RecalculateBounds();
            return mesh;
        }
    }
}
