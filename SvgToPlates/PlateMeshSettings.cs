using System;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>How the boundary tangent vectors are extended to the other vertices.</summary>
    public enum TangentMethod
    {
        /// <summary>Barycentric (P1) interpolation over the final triangles; unknown vertex values come from a discrete Laplace solve. Smoothest; also covers the offset band.</summary>
        MeshHarmonic,

        /// <summary>Hormann-Floater mean value coordinates inside the original outline; closest-boundary-point tangent in the offset band.</summary>
        MeanValue
    }

    /// <summary>How the ribbon meshes along the shared paths are built.</summary>
    public enum RibbonMethod
    {
        /// <summary>Copies of the path shifted to each side along its vertex normals (smoothed over ribbonNormalSmoothing), connected column by column. UV1.x is exactly the arc length of each column's path vertex, UV1.y exactly the offset.</summary>
        ShiftedPath,

        /// <summary>Like ShiftedPath, but every vertex moves along the same direction: the path's average normal (perpendicular to the line from its first to its last point). The ribbon narrows where the path bends away from that direction.</summary>
        ShiftedAlongAverageNormal,

        /// <summary>The original method: the region within ribbonExtent of the path (see ribbonJoin), freely triangulated; UV1.x away from the path is a smooth fill-in, which bends near the path's ends.</summary>
        Outline
    }

    /// <summary>How a ribbon's edge turns around the outside of a corner of its path (Outline method only).</summary>
    public enum RibbonJoin
    {
        /// <summary>A circular arc: every point of the edge is exactly ribbonExtent from the path, so there are no spikes.</summary>
        Round,

        /// <summary>The two edges are connected by a straight line, cutting the corner.</summary>
        Bevel,

        /// <summary>The edges are extended to a point. Corners sharper than maxMiterScale allows fall back to a bevel.</summary>
        Miter
    }

    /// <summary>All user-adjustable parameters of the SVG -> plate mesh conversion.</summary>
    [Serializable]
    public class PlateMeshSettings
    {
        [Header("Units")]
        [Tooltip("World units per SVG unit. Adobe Illustrator writes SVG in points (1/72 inch), so 0.000352778 (= 0.0254 / 72) gives metres. Use 0.001 for an SVG drawn in millimetres.")]
        [Min(1e-9f)] public float svgUnitsToWorld = 0.000352778f;

        [Tooltip("Maximum distance (world units) between a Bezier curve and its straight-segment approximation.")]
        [Min(1e-6f)] public float curveTolerance = 0.001f;

        [Tooltip("Consecutive outline points closer than this (world units) are merged.")]
        [Min(0f)] public float weldDistance = 0.0005f;

        [Header("Layout")]
        [Tooltip("After the minimum-area box is found, make its long side the x axis.")]
        public bool longSideAlongX = true;

        [Header("Offset")]
        [Tooltip("Distance (world units) the polygons are grown outward. 0 = no offset band.")]
        [Min(0f)] public float offset = 0.2f;

        [Tooltip("Plate offsets: corners are mitred, but the offset vertex never travels further than offset * this value. Also the limit for ribbon Miter joins, where sharper corners become a bevel.")]
        [Min(1f)] public float maxMiterScale = 2f;

        [Header("Triangulation")]
        [Tooltip("Spacing (world units) of extra interior points on a world-aligned grid. 0 = boundary vertices only.")]
        [Min(0f)] public float interiorSpacing = 0f;

        [Header("Mesh refinement (BurstTriangulator)")]
        [Tooltip("Insert extra vertices until every triangle meets the angle and area limits below. " +
                 "Also enables RestoreBoundary so only the inside of each plate is refined. " +
                 "With refinement on, interiorSpacing can usually be set to 0.")]
        public bool refineMesh = true;

        [Tooltip("Refinement removes triangles with an angle smaller than this (degrees). Values above about 33 may not terminate.")]
        [Range(0f, 34f)] public float minimumAngleDegrees = 20f;

        [Tooltip("Triangles with an area larger than this (world units squared) are split.")]
        [Min(1e-6f)] public float maximumTriangleArea = 0.01f;

        [Tooltip("Concentric-shells parameter used when refinement splits constrained (boundary) segments. The package default is 0.001.")]
        [Min(1e-6f)] public float concentricShellsParameter = 0.001f;

        [Header("Convex hulls")]
        [Tooltip("Also create plate_hull_AA for every plate: the convex hull of its original (not offset) outline, triangulated from the hull corners only. Same layout and UV0 as the plates.")]
        public bool generateHulls = true;

        [Tooltip("Thickness (world units) of the hull meshes. Above 0 each hull is a prism from y = 0 down to y = -thickness, so a MeshCollider marked Convex follows it (PhysX cannot build a convex collider from flat geometry and falls back to a box). 0 = flat.")]
        [Min(0f)] public float hullThickness = 0.01f;

        [Header("Ribbons (where plates abut)")]
        [Tooltip("Also create a ribbon mesh along every path where two plates abut, named ribbon_AA_BB after the two plates.")]
        public bool generateRibbons = true;

        [Tooltip("ShiftedPath: the path shifted to each side along its smoothed normals, so UV1 is exact (recommended). ShiftedAlongAverageNormal: shifted along one direction for the whole path. Outline: the original free-form outline.")]
        public RibbonMethod ribbonMethod = RibbonMethod.ShiftedPath;

        [Tooltip("ShiftedPath: how far along the path (world units, Gaussian standard deviation) the vertex normals are averaged to get the direction the ribbon is shifted in. Larger keeps neighbouring columns from crossing at tight bends but lets the width vary a little more. 0 = each vertex's own normal. Geometry only: the normal written to UV2 is always the path vertex's own normal.")]
        [Min(0f)] public float ribbonNormalSmoothing = 0.05f;

        [Tooltip("Outline method only: how the ribbon edge turns around the outside of a corner of its path. Round keeps every edge point exactly ribbonExtent from the path.")]
        public RibbonJoin ribbonJoin = RibbonJoin.Round;

        [Tooltip("How far (world units) each ribbon extends to each side of the shared path. The total width is twice this.")]
        [Min(0f)] public float ribbonExtent = 0.2f;

        [Header("Normal field (UV1)")]
        public TangentMethod tangentMethod = TangentMethod.MeshHarmonic;

        [Tooltip("Where the interpolated vector is shorter than this (0 to 1), the normal is set to (0, 0). " +
                 "On the boundary and in the offset band the length is close to 1; it shrinks towards the middle of the plate.")]
        [Range(0f, 1f)] public float minNormalLength = 0.001f;
    }
}
