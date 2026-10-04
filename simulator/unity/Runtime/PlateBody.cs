using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace Tectonic.UnityBridge
{
    /// <summary>
    /// One GameObject per plate: a Rigidbody2D plus a collider that exists ONLY to bump into
    /// the walls.
    ///
    /// The collider is a convex hull of the plate's cells, not its true outline. Plates are
    /// concave and change shape every step, and Unity would have to decompose a concave
    /// PolygonCollider2D into convex pieces every rebuild. Since this collider never touches
    /// another plate — the layer matrix forbids it — and the walls sit off-screen where
    /// contact geometry is invisible, a hull refreshed every ~30 steps is more than enough.
    /// </summary>
    [RequireComponent(typeof(Rigidbody2D))]
    public sealed class PlateBody : MonoBehaviour
    {
        public Plate Plate { get; private set; }
        public Rigidbody2D Body { get; private set; }
        public PolygonCollider2D WallCollider { get; private set; }

        /// <summary>Set when the core has produced new cell data the renderer hasn't drawn yet.</summary>
        public bool Dirty { get; private set; }

        static readonly List<Vector2> HullScratch = new List<Vector2>(4096);
        static readonly List<Vector2> PointScratch = new List<Vector2>(16384);

        public static PlateBody Create(TectonicRunner runner, Plate plate)
        {
            var go = new GameObject($"Plate {plate.Id}");
            go.transform.SetParent(runner.transform, false);
            if (TectonicSetup.PlateLayer >= 0) go.layer = TectonicSetup.PlateLayer;

            var body = go.AddComponent<PlateBody>();
            body.Plate = plate;

            var rb = go.GetComponent<Rigidbody2D>();
            rb.gravityScale = 0f;

            // Kinematic rather than Static: an anchored plate must never respond to forces,
            // but it still has to be movable by script for the slow reset animation. Unity
            // still resolves dynamic-vs-kinematic contacts, so the mobile plates stop against
            // it exactly as they should.
            rb.bodyType = plate.IsStatic ? RigidbodyType2D.Kinematic : RigidbodyType2D.Dynamic;
            rb.interpolation = RigidbodyInterpolation2D.Interpolate;
            rb.collisionDetectionMode = CollisionDetectionMode2D.Continuous;

            // Drag lives in the core (Plate.Force), so Unity's own damping stays off —
            // doing it in both places damps the plates to a standstill.
            rb.linearDamping = 0f;
            rb.angularDamping = 0f;

            // Rigidbody2D.position only reaches the Transform on the next physics step, so set
            // the transform too or every plate renders stacked at the origin for one frame.
            go.transform.SetLocalPositionAndRotation(
                new Vector3(plate.Position.x, plate.Position.y, 0f),
                Quaternion.Euler(0f, 0f, plate.Rotation * Mathf.Rad2Deg));

            rb.useAutoMass = false;   // otherwise mass/centreOfMass writes are silently ignored
            rb.position = new Vector2(plate.Position.x, plate.Position.y);
            rb.rotation = plate.Rotation * Mathf.Rad2Deg;
            // Anchored plates start at rest: Unity integrates a kinematic body's velocity, so
            // seeding one here would make the immovable centre drift forever, unopposed.
            rb.linearVelocity = plate.IsStatic
                ? Vector2.zero
                : new Vector2(plate.LinearVelocity.x, plate.LinearVelocity.y);

            body.Body = rb;
            body.WallCollider = go.AddComponent<PolygonCollider2D>();
            // Non-trigger: real contacts, but only ever against walls (see the layer matrix).
            body.RebuildWallCollider(runner.World.Config.CellSize);

            go.AddComponent<PlateMeshRenderer>().Bind(body, runner);
            return body;
        }

        public void MarkDirty() => Dirty = true;
        public void ClearDirty() => Dirty = false;

        public void RebuildWallCollider(float cellSize)
        {
            if (Plate == null || WallCollider == null) return;

            PointScratch.Clear();
            var cells = Plate.Read;
            for (int i = 0; i < cells.Length; i++)
            {
                if (!cells[i].Alive) continue;
                float2 local = HexLattice.AxialToLocal(Plate.AxialAt(i), cellSize);
                PointScratch.Add(new Vector2(local.x, local.y));
            }

            // A plate worn down to a sliver would otherwise keep its old, oversized hull and
            // go on bumping the walls.
            if (PointScratch.Count < 3) { WallCollider.enabled = false; return; }
            WallCollider.enabled = true;

            ConvexHull(PointScratch, HullScratch);
            if (HullScratch.Count >= 3) WallCollider.points = HullScratch.ToArray();
        }

        /// <summary>Andrew's monotone chain. O(n log n), no allocation beyond the output list.</summary>
        static void ConvexHull(List<Vector2> points, List<Vector2> hull)
        {
            points.Sort((a, b) => a.x != b.x ? a.x.CompareTo(b.x) : a.y.CompareTo(b.y));
            hull.Clear();

            int n = points.Count;
            var tmp = new Vector2[2 * n];
            int k = 0;

            for (int i = 0; i < n; i++)
            {
                while (k >= 2 && Cross(tmp[k - 2], tmp[k - 1], points[i]) <= 0) k--;
                tmp[k++] = points[i];
            }

            int lower = k + 1;
            for (int i = n - 2; i >= 0; i--)
            {
                while (k >= lower && Cross(tmp[k - 2], tmp[k - 1], points[i]) <= 0) k--;
                tmp[k++] = points[i];
            }

            for (int i = 0; i < k - 1; i++) hull.Add(tmp[i]);
        }

        static float Cross(Vector2 o, Vector2 a, Vector2 b) =>
            (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    }
}
