using System;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>
    /// SVG text in, one Mesh per outline out. Each stage is its own method so it is easy to follow and to change.
    ///
    /// Mesh layout: vertices lie in the xz plane (y = 0), normals are +y.
    /// UV0 (mesh.uv, Vector2) is for texture mapping: the shared bounding box of all plates maps so that its long
    ///     side spans 0..1 and its short side is centred on 0.5 (u follows x, v follows z, square texels).
    /// UV1 (Vector4, set with SetUVs(1, ...)):
    ///     x, y = outward normal (x and z components): the tangent field rotated 90 degrees clockwise and
    ///            normalised, or (0, 0) where the field is shorter than minNormalLength
    ///     z    = length of the field vector before normalising (1 on the boundary, smaller towards the middle)
    ///     w    = signed distance from the ORIGINAL boundary (negative inside, 0 on it, positive in the offset band)
    /// </summary>
    public static class PlateMeshPipeline
    {
        // Output order: plate_00.., then plate_hull_00.. (convex hull of each plate's original outline), then ribbons.
        public static List<Mesh> Build(string svgText, PlateMeshSettings settings)
        {
            // The curve tolerance is given in world units but flattening happens on raw SVG coordinates.
            float curveToleranceInSvgUnits = settings.curveTolerance / settings.svgUnitsToWorld;
            List<float2[]> svgOutlines = SvgPolygonReader.Read(svgText, curveToleranceInSvgUnits);

            List<float2[]> outlines = PrepareOutlines(svgOutlines, settings);
            OrientedBox box = FindSharedBox(outlines, settings.longSideAlongX);
            List<float2[]> alignedOutlines = AlignToBox(outlines, box);

            var meshes = new List<Mesh>();
            for (int i = 0; i < alignedOutlines.Count; i++)
            {
                Mesh mesh = BuildPlateMesh(alignedOutlines[i], settings, box.size, $"plate_{i:00}");
                if (mesh != null) meshes.Add(mesh);
            }

            if (settings.generateHulls)
            {
                for (int i = 0; i < alignedOutlines.Count; i++)
                {
                    Mesh hull = BuildHullMesh(alignedOutlines[i], settings, box.size, $"plate_hull_{i:00}");
                    if (hull != null) meshes.Add(hull);
                }
            }

            if (settings.generateRibbons)
                meshes.AddRange(RibbonBuilder.BuildAll(alignedOutlines, settings, box.size));
            return meshes;
        }

        // ================================================================ stage 1: outlines

        /// <summary>SVG units -> world units, y flipped (SVG y points down, world z points up), welded, counter-clockwise.</summary>
        static List<float2[]> PrepareOutlines(List<float2[]> svgOutlines, PlateMeshSettings settings)
        {
            var result = new List<float2[]>();
            foreach (float2[] svgOutline in svgOutlines)
            {
                var world = new float2[svgOutline.Length];
                for (int i = 0; i < world.Length; i++)
                    world[i] = new float2(svgOutline[i].x, -svgOutline[i].y) * settings.svgUnitsToWorld;

                float2[] cleaned = PlateGeometry.RemoveClosePoints(world, settings.weldDistance);
                if (cleaned.Length < 3)
                {
                    Debug.LogWarning("Skipped an outline with fewer than three distinct points.");
                    continue;
                }
                result.Add(PlateGeometry.EnsureCounterClockwise(cleaned));
            }
            return result;
        }

        // ================================================================ stage 2: shared rotation + centring

        /// <summary>Smallest-area rectangle (rotating calipers) around the convex hull of every outline.</summary>
        static OrientedBox FindSharedBox(List<float2[]> outlines, bool longSideAlongX)
        {
            var allPoints = new List<float2>();
            foreach (float2[] outline in outlines)
                allPoints.AddRange(outline);

            float2[] hull = PlateGeometry.ConvexHull(allPoints);
            OrientedBox box = PlateGeometry.FindMinimumAreaBox(hull);
            if (longSideAlongX) box = PlateGeometry.WithLongSideAlongX(box);

            float angle = math.degrees(math.atan2(box.axisX.y, box.axisX.x));
            Debug.Log($"Shared bounding box {box.size.x:F3} x {box.size.y:F3}; its x axis pointed {angle:F2} degrees before alignment.");
            return box;
        }

        static List<float2[]> AlignToBox(List<float2[]> outlines, OrientedBox box)
        {
            var result = new List<float2[]>(outlines.Count);
            foreach (float2[] outline in outlines)
            {
                var aligned = new float2[outline.Length];
                for (int i = 0; i < aligned.Length; i++)
                    aligned[i] = PlateGeometry.ToBoxSpace(outline[i], box);
                result.Add(aligned);
            }
            return result;
        }

        // ================================================================ stage 3: one plate

        static Mesh BuildPlateMesh(float2[] ring, PlateMeshSettings settings, float2 boxSize, string name)
        {
            // Outer outline: the offset ring, or the original ring itself when there is no offset.
            float2[] outerRing = ring;
            if (settings.offset > 2f * settings.weldDistance)
                outerRing = BuildOffsetRing(ring, settings);
            bool hasOffsetBand = !ReferenceEquals(outerRing, ring);

            if (hasOffsetBand && !Encloses(outerRing, ring))
            {
                Debug.LogError($"{name}: the offset outline does not enclose the original outline (offset {settings.offset} is too large for this shape). Plate skipped.");
                return null;
            }

            List<float2> interiorPoints = SampleInteriorPoints(ring, outerRing, settings.interiorSpacing);

            // Vertex layout: [original ring][offset ring][interior points]
            var positions = new List<float2>(ring);
            if (hasOffsetBand) positions.AddRange(outerRing);
            positions.AddRange(interiorPoints);

            var constraintEdges = new List<int>();
            AddRingConstraints(constraintEdges, 0, ring.Length);
            if (hasOffsetBand) AddRingConstraints(constraintEdges, ring.Length, outerRing.Length);

            if (PlateGeometry.HasCrossingEdges(positions, constraintEdges))
            {
                Debug.LogError($"{name}: the offset outline crosses the original outline. Try a smaller offset. Plate skipped.");
                return null;
            }

            PlateTriangulator.Result triangulation = PlateTriangulator.Triangulate(positions, constraintEdges, settings);
            List<int> triangles = PlateTriangulator.KeepInsideTrianglesFacingUp(triangulation.positions, triangulation.triangles, outerRing);
            if (triangles.Count == 0)
            {
                Debug.LogError($"{name}: triangulation produced no triangles. Plate skipped.");
                return null;
            }

            // Attributes are computed for every vertex first, then vertices no triangle uses are dropped.
            VertexAttributes attributes = ComputeVertexAttributes(triangulation.positions, triangles, ring, boxSize, settings);
            int[] sourceIndex = CompactTriangles(triangles, triangulation.positions.Count);
            float2[] finalPositions = Gather(triangulation.positions, sourceIndex);
            Vector2[] textureUvs = Gather(attributes.textureUvs, sourceIndex);
            Vector4[] fieldUvs = Gather(attributes.fieldUvs, sourceIndex);

            Mesh mesh = CreateMesh(name, finalPositions, textureUvs, fieldUvs, triangles);
            int addedByRefinement = triangulation.positions.Count - positions.Count;
            Debug.Log($"{name}: offset {settings.offset}, {finalPositions.Length} vertices ({ring.Length} boundary, {(hasOffsetBand ? outerRing.Length : 0)} offset, {interiorPoints.Count} interior, {addedByRefinement} added by refinement), {triangles.Count / 3} triangles.");
            return mesh;
        }

        /// <summary>
        /// The convex hull of the original (not offset) outline, Delaunay-triangulated from the hull corners alone
        /// (no interior points, no refinement). Same layout and UV0 as the plates; no UV1.
        /// </summary>
        static Mesh BuildHullMesh(float2[] ring, PlateMeshSettings settings, float2 boxSize, string name)
        {
            float2[] hull = PlateGeometry.ConvexHull(ring);
            if (hull.Length < 3 || PlateGeometry.SignedArea(hull) <= 0f)
            {
                Debug.LogError($"{name}: the outline's convex hull has no area. Hull skipped.");
                return null;
            }

            PlateTriangulator.Result triangulation = PlateTriangulator.Triangulate(hull, new List<int>(), settings, allowRefinement: false);
            List<int> triangles = PlateTriangulator.KeepInsideTrianglesFacingUp(triangulation.positions, triangulation.triangles, hull);
            if (triangles.Count == 0)
            {
                Debug.LogError($"{name}: triangulation produced no triangles. Hull skipped.");
                return null;
            }

            int[] sourceIndex = CompactTriangles(triangles, triangulation.positions.Count);
            float2[] positions = Gather(triangulation.positions, sourceIndex);

            float textureScale = 1f / math.max(boxSize.x, boxSize.y);
            var textureUvs = new Vector2[positions.Length];
            for (int i = 0; i < positions.Length; i++)
                textureUvs[i] = TextureCoordinate(positions[i], textureScale);

            Mesh mesh = CreateMesh(name, positions, textureUvs, null, triangles);
            Debug.Log($"{name}: {positions.Length} hull vertices (outline had {ring.Length}), {triangles.Count / 3} triangles.");
            return mesh;
        }

        static bool Encloses(float2[] outer, float2[] inner)
        {
            foreach (float2 p in inner)
                if (!PlateGeometry.IsInside(p, outer)) return false;
            return true;
        }

        /// <summary>Grows the ring outward, removes small loops caused by short edges at concave corners, merges close points.</summary>
        static float2[] BuildOffsetRing(float2[] ring, PlateMeshSettings settings)
        {
            float2[] offsetRing = PlateGeometry.OffsetOutward(ring, settings.offset, settings.maxMiterScale);
            offsetRing = PlateGeometry.RemoveSelfIntersectionLoops(offsetRing);
            return PlateGeometry.RemoveClosePoints(offsetRing, settings.weldDistance);
        }

        /// <summary>
        /// Grid points (aligned to the world axes) that lie inside the outer outline and at least half a grid cell
        /// away from both outlines, so they do not create slivers against the boundary vertices.
        /// </summary>
        static List<float2> SampleInteriorPoints(float2[] innerRing, float2[] outerRing, float spacing)
        {
            var points = new List<float2>();
            if (spacing <= 0f) return points;

            float2 min = outerRing[0];
            float2 max = outerRing[0];
            foreach (float2 p in outerRing)
            {
                min = math.min(min, p);
                max = math.max(max, p);
            }

            float minGap = 0.5f * spacing;
            bool hasOffsetBand = !ReferenceEquals(innerRing, outerRing);

            int firstX = (int)math.ceil(min.x / spacing), lastX = (int)math.floor(max.x / spacing);
            int firstZ = (int)math.ceil(min.y / spacing), lastZ = (int)math.floor(max.y / spacing);
            for (int ix = firstX; ix <= lastX; ix++)
            {
                for (int iz = firstZ; iz <= lastZ; iz++)
                {
                    var p = new float2(ix, iz) * spacing;
                    if (!PlateGeometry.IsInside(p, outerRing)) continue;
                    if (PlateGeometry.DistanceToRing(p, innerRing) < minGap) continue;
                    if (hasOffsetBand && PlateGeometry.DistanceToRing(p, outerRing) < minGap) continue;
                    points.Add(p);
                }
            }
            return points;
        }

        static void AddRingConstraints(List<int> edgePairs, int firstVertex, int vertexCount)
        {
            for (int i = 0; i < vertexCount; i++)
            {
                edgePairs.Add(firstVertex + i);
                edgePairs.Add(firstVertex + (i + 1) % vertexCount);
            }
        }

        // ================================================================ stage 4: UVs

        struct VertexAttributes
        {
            public Vector2[] textureUvs;   // UV channel 0 (mesh.uv)
            public Vector4[] fieldUvs;     // UV channel 1 (mesh.uv2)
        }

        /// <summary>
        /// Computes both UV sets for every vertex.
        /// Any vertex lying on the original boundary (the original vertices, and points refinement may add on it)
        /// gets signed distance 0 and the boundary tangent at its position. Every other vertex gets the exact signed
        /// distance to the original boundary and a tangent vector from the selected interpolation method.
        /// </summary>
        static VertexAttributes ComputeVertexAttributes(List<float2> positions, List<int> triangles, float2[] ring, float2 boxSize, PlateMeshSettings settings)
        {
            float2[] ringTangents = TangentField.ComputeBoundaryTangents(ring);
            var boundary = new BoundaryLookup(ring, ringTangents);

            int count = positions.Count;
            var used = new bool[count];
            foreach (int index in triangles) used[index] = true;

            var signedDistance = new float[count];
            var boundaryTangent = new float2[count];   // tangent at the closest boundary point, valid for every vertex
            var isFixed = new bool[count];             // vertices whose tangent is known (not solved for)
            for (int i = 0; i < count; i++)
            {
                float distance = boundary.Closest(positions[i], out boundaryTangent[i]);
                bool onBoundary = distance < PlateGeometry.OnBoundaryTolerance;
                bool inside = !onBoundary && PlateGeometry.IsInside(positions[i], ring);
                signedDistance[i] = onBoundary ? 0f : (inside ? -distance : distance);
                isFixed[i] = onBoundary || !used[i];
            }

            float2[] tangents = InterpolateTangents(positions, triangles, ring, ringTangents, boundaryTangent, isFixed, settings);

            float textureScale = 1f / math.max(boxSize.x, boxSize.y);
            var attributes = new VertexAttributes { textureUvs = new Vector2[count], fieldUvs = new Vector4[count] };
            for (int i = 0; i < count; i++)
            {
                attributes.textureUvs[i] = TextureCoordinate(positions[i], textureScale);
                attributes.fieldUvs[i] = FieldCoordinate(tangents[i], signedDistance[i], settings.minNormalLength);
            }
            return attributes;
        }

        /// <summary>
        /// UV0. The shared bounding box is centred on the origin, so scaling by 1 / (long side) puts the long
        /// coordinate in 0..1 and centres the short one on 0.5 (spanning 0.5 -+ shortSide / longSide / 2).
        /// </summary>
        internal static Vector2 TextureCoordinate(float2 position, float scale)
        {
            return new Vector2(0.5f + position.x * scale, 0.5f + position.y * scale);
        }

        /// <summary>
        /// UV1: (normal.x, normal.z, field length, signed distance). The outward normal is the tangent rotated 90 degrees
        /// clockwise (seen from above, x right, z up), which points outward on the counter-clockwise outlines.
        /// Where the field is too short to give a reliable direction the normal is set to (0, 0).
        /// </summary>
        static Vector4 FieldCoordinate(float2 tangent, float signedDistance, float minNormalLength)
        {
            float length = math.length(tangent);
            float2 normal = float2.zero;
            if (length > 0f && length >= minNormalLength)
                normal = new float2(tangent.y, -tangent.x) / length;
            return new Vector4(normal.x, normal.y, length, signedDistance);
        }

        /// <summary>Returns a tangent vector for every vertex using the method chosen in the settings. Free vertices are not normalised.</summary>
        static float2[] InterpolateTangents(List<float2> positions, List<int> triangles, float2[] ring, float2[] ringTangents,
                                            float2[] boundaryTangent, bool[] isFixed, PlateMeshSettings settings)
        {
            switch (settings.tangentMethod)
            {
                case TangentMethod.MeshHarmonic:
                    return HarmonicTangents.Solve(positions, triangles, isFixed, boundaryTangent);

                case TangentMethod.MeanValue:
                {
                    // Inside the outline: mean value coordinates. Offset band: tangent at the closest boundary point.
                    var result = (float2[])boundaryTangent.Clone();
                    for (int i = 0; i < positions.Count; i++)
                    {
                        if (isFixed[i] || !PlateGeometry.IsInside(positions[i], ring)) continue;
                        result[i] = MeanValueTangents.Interpolate(positions[i], ring, ringTangents);
                    }
                    return result;
                }

                default:
                    throw new ArgumentOutOfRangeException(nameof(settings.tangentMethod));
            }
        }

        // ================================================================ stage 5: Mesh

        /// <summary>
        /// Rewrites the triangle indices so they only refer to vertices that some kept triangle uses (numbered in order
        /// of first use). Returns, for every new vertex, the index it had before.
        /// </summary>
        internal static int[] CompactTriangles(List<int> triangles, int vertexCount)
        {
            var newIndex = new int[vertexCount];
            for (int i = 0; i < newIndex.Length; i++) newIndex[i] = -1;

            var sourceIndex = new List<int>();
            for (int t = 0; t < triangles.Count; t++)
            {
                int old = triangles[t];
                if (newIndex[old] < 0)
                {
                    newIndex[old] = sourceIndex.Count;
                    sourceIndex.Add(old);
                }
                triangles[t] = newIndex[old];
            }
            return sourceIndex.ToArray();
        }

        internal static T[] Gather<T>(IReadOnlyList<T> source, int[] sourceIndex)
        {
            var result = new T[sourceIndex.Length];
            for (int i = 0; i < result.Length; i++) result[i] = source[sourceIndex[i]];
            return result;
        }

        static Mesh CreateMesh(string name, float2[] positions, Vector2[] textureUvs, Vector4[] fieldUvs, List<int> triangles)
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
            mesh.SetUVs(0, new List<Vector2>(textureUvs));   // mesh.uv
            if (fieldUvs != null)
                mesh.SetUVs(1, new List<Vector4>(fieldUvs)); // second UV channel, four components
            mesh.SetTriangles(triangles, 0);
            mesh.RecalculateBounds();
            return mesh;
        }
    }
}
