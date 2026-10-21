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

    /// <summary>How a ribbon's edge turns around the outside of a corner of its path.</summary>
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
        [Tooltip("World units per SVG unit. 0.001 means the SVG is in millimetres and the meshes come out in metres.")]
        [Min(1e-9f)] public float svgUnitsToWorld = 0.001f;

        [Tooltip("Maximum distance (world units) between a Bezier curve and its straight-segment approximation.")]
        [Min(1e-6f)] public float curveTolerance = 0.001f;

        [Tooltip("Consecutive outline points closer than this (world units) are merged.")]
        [Min(0f)] public float weldDistance = 0.0005f;

        [Header("Layout")]
        [Tooltip("After the minimum-area box is found, make its long side the x axis.")]
        public bool longSideAlongX = true;

        [Header("Offset")]
        [Tooltip("Distance (world units) the polygons are grown outward. 0 = no offset band.")]
        [Min(0f)] public float offset = 0.1f;

        [Tooltip("Plate offsets: corners are mitred, but the offset vertex never travels further than offset * this value. Also the limit for ribbon Miter joins, where sharper corners become a bevel.")]
        [Min(1f)] public float maxMiterScale = 2f;

        [Header("Triangulation")]
        [Tooltip("Spacing (world units) of extra interior points on a world-aligned grid. 0 = boundary vertices only.")]
        [Min(0f)] public float interiorSpacing = 0.2f;

        [Header("Mesh refinement (BurstTriangulator)")]
        [Tooltip("Insert extra vertices until every triangle meets the angle and area limits below. " +
                 "Also enables RestoreBoundary so only the inside of each plate is refined. " +
                 "With refinement on, interiorSpacing can usually be set to 0.")]
        public bool refineMesh = false;

        [Tooltip("Refinement removes triangles with an angle smaller than this (degrees). Values above about 33 may not terminate.")]
        [Range(0f, 34f)] public float minimumAngleDegrees = 20f;

        [Tooltip("Triangles with an area larger than this (world units squared) are split.")]
        [Min(1e-6f)] public float maximumTriangleArea = 0.01f;

        [Tooltip("Concentric-shells parameter used when refinement splits constrained (boundary) segments. The package default is 0.001.")]
        [Min(1e-6f)] public float concentricShellsParameter = 0.001f;

        [Header("Ribbons (where plates abut)")]
        [Tooltip("Also create a ribbon mesh along every path where two plates abut, named ribbon_AA_BB after the two plates.")]
        public bool generateRibbons = true;

        [Tooltip("How the ribbon edge turns around the outside of a corner of its path. Round keeps every edge point exactly ribbonExtent from the path.")]
        public RibbonJoin ribbonJoin = RibbonJoin.Round;

        [Tooltip("How far (world units) each ribbon extends to each side of the shared path. The total width is twice this.")]
        [Min(0f)] public float ribbonExtent = 0.1f;

        [Header("Normal field (UV1)")]
        public TangentMethod tangentMethod = TangentMethod.MeshHarmonic;

        [Tooltip("Where the interpolated vector is shorter than this (0 to 1), the normal is set to (0, 0). " +
                 "On the boundary and in the offset band the length is close to 1; it shrinks towards the middle of the plate.")]
        [Range(0f, 1f)] public float minNormalLength = 0.2f;
    }
}
