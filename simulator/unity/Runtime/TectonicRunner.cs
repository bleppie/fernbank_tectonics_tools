using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace Tectonic.UnityBridge
{
    /// <summary>
    /// THE SEAM. The only place where the simulation and Unity meet.
    ///
    /// Per fixed step:
    ///   1. mirror each Rigidbody2D's transform and velocity INTO the core
    ///   2. run the geology (at a lower rate than physics — see GeologyEveryNFixedSteps)
    ///   3. apply the core's force and torque back onto the Rigidbody2D
    ///
    /// The core never writes a transform and never reads a MonoBehaviour. Everything below
    /// this file could be deleted and the simulation would still run headlessly, which is
    /// exactly how it is tested.
    ///
    /// Note the deliberate division of labour with Unity's physics:
    ///   - Unity owns integration and PLATE-VS-WALL contacts (non-penetration is correct there)
    ///   - the core owns PLATE-VS-PLATE overlap (interpenetration is REQUIRED — subduction is
    ///     one plate sliding under another, and a contact solver exists to prevent that)
    /// The layer collision matrix is what keeps those two apart. See TectonicSetup.
    /// </summary>
    [DefaultExecutionOrder(-100)]
    public sealed class TectonicRunner : MonoBehaviour
    {
        [Header("Simulation")]
        public TectonicSettings Settings;

        [Tooltip("Run the geology once every N fixed steps. Physics wants ~50 Hz; plate " +
                 "tectonics does not. 5 gives a 10 Hz geology tick against a 50 Hz solver.")]
        [Min(1)] public int GeologyEveryNFixedSteps = 5;

        [Tooltip("Rebuild each plate's wall collider every N geology steps. The exact outline " +
                 "is geologically meaningless off-screen, so a convex hull refreshed " +
                 "occasionally is plenty.")]
        [Min(1)] public int ColliderRebuildInterval = 30;

        public World World { get; private set; }
        public IReadOnlyList<PlateBody> Bodies => _bodies;

        readonly List<PlateBody> _bodies = new List<PlateBody>();
        int _fixedCounter;
        int _geologyCounter;

        void Awake()
        {
            if (transform.localToWorldMatrix != Matrix4x4.identity)
                Debug.LogError("[Tectonic] TectonicRunner must sit at the origin with identity " +
                               "rotation and scale — the core works in world units, and the walls " +
                               "are built in this transform's local space.");

            var cfg = Settings != null ? Settings.ToConfig() : new TectonicConfig();

            World = Settings == null || Settings.RadialScenario
                ? WorldBuilder.RadialScenario(cfg, Settings != null ? Settings.Satellites : 6)
                : WorldBuilder.DefaultScenario(cfg);

            TectonicSetup.BuildWalls(transform, cfg);
            SyncBodies();

            // Before the first ApplyForces. Otherwise plates spend the first few fixed steps at
            // Unity's default mass of 1 while the core pushes forces scaled for their real mass,
            // which is a large spurious kick at start-up.
            PushMassProperties();
        }

        void FixedUpdate()
        {
            if (World == null) return;

            MirrorTransformsIntoCore();

            if (++_fixedCounter >= GeologyEveryNFixedSteps)
            {
                _fixedCounter = 0;
                World.Step();
                _geologyCounter++;

                SyncBodies();

                // Rebuild colliders BEFORE pushing mass properties: changing a body's collider
                // makes Unity recompute its derived inertia, which would silently discard the
                // values we just set.
                if (_geologyCounter % ColliderRebuildInterval == 0)
                    foreach (var b in _bodies) b.RebuildWallCollider(World.Config.CellSize);

                PushMassProperties();

                foreach (var b in _bodies) b.MarkDirty();
            }

            ApplyForces();
        }

        /// <summary>Rigidbody2D -> core. The core reads these and never writes them.</summary>
        void MirrorTransformsIntoCore()
        {
            foreach (var body in _bodies)
            {
                if (body.Plate == null) continue;
                var rb = body.Body;
                body.Plate.Position = (float2)(Vector2)rb.position;
                body.Plate.Rotation = rb.rotation * Mathf.Deg2Rad;

                // An anchored plate reports zero velocity regardless of what the body says.
                // Unity DOES integrate a kinematic body's velocity, so a stray non-zero value
                // would let the "immovable" centre drift — matching HeadlessHost, which zeroes
                // these explicitly, so offline verification stays representative.
                if (body.Plate.IsStatic)
                {
                    body.Plate.LinearVelocity = float2.zero;
                    body.Plate.AngularVelocity = 0f;
                }
                else
                {
                    body.Plate.LinearVelocity = (float2)(Vector2)rb.linearVelocity;
                    body.Plate.AngularVelocity = rb.angularVelocity * Mathf.Deg2Rad;
                }
            }
        }

        /// <summary>core -> Rigidbody2D. The core's entire physical output.</summary>
        void ApplyForces()
        {
            foreach (var body in _bodies)
            {
                if (body.Plate == null || body.Plate.IsStatic) continue;
                var p = body.Plate;
                body.Body.AddForce(new Vector2(p.Force.x, p.Force.y), ForceMode2D.Force);
                body.Body.AddTorque(p.Torque, ForceMode2D.Force);
            }
        }

        /// <summary>
        /// Unity derives mass and inertia from colliders, but ours is a stand-in hull that
        /// says nothing about crust thickness. The core computes the real values from the
        /// rock columns, so push them across explicitly.
        /// </summary>
        void PushMassProperties()
        {
            foreach (var body in _bodies)
            {
                if (body.Plate == null || body.Plate.IsStatic) continue;
                var rb = body.Body;
                rb.mass = body.Plate.Mass;
                rb.centerOfMass = new Vector2(body.Plate.CenterOfMassLocal.x, body.Plate.CenterOfMassLocal.y);
                rb.inertia = body.Plate.MomentOfInertia;
            }
        }

        /// <summary>
        /// Keep one PlateBody per core plate. Plates can appear (rifting) and disappear
        /// (fully subducted), so this reconciles rather than assuming a fixed set.
        /// </summary>
        void SyncBodies()
        {
            for (int i = _bodies.Count - 1; i >= 0; i--)
            {
                if (_bodies[i] == null) { _bodies.RemoveAt(i); continue; }
                if (!World.Plates.Contains(_bodies[i].Plate))
                {
                    Destroy(_bodies[i].gameObject);
                    _bodies.RemoveAt(i);
                }
            }

            // Index loop rather than List.Exists: the closure would allocate per plate per tick.
            foreach (var plate in World.Plates)
            {
                bool found = false;
                for (int i = 0; i < _bodies.Count && !found; i++) found = _bodies[i].Plate == plate;
                if (!found) _bodies.Add(PlateBody.Create(this, plate));
            }
        }

        // ---- Exhibit control --------------------------------------------------------

        /// <summary>
        /// Visitor input: a shove applied at a world point. Goes straight to the rigid body —
        /// the geology finds out about it through the resulting motion, which is exactly how
        /// the ambient ridge push works too.
        /// </summary>
        public void PushAt(Vector2 worldPoint, Vector2 force)
        {
            foreach (var body in _bodies)
            {
                if (body.Plate == null || body.Plate.IsStatic) continue;   // anchored plates ignore shoves
                if (body.Plate.TryCellAtWorld((float2)worldPoint, World.Config.CellSize, out _, out _))
                {
                    body.Body.AddForceAtPosition(force, worldPoint, ForceMode2D.Force);
                    return;
                }
            }
        }

        /// <summary>True when the map has stopped evolving on its own. The exhibit's cue to
        /// begin the slow reset — more reliable than a presence sensor alone, since a visitor
        /// can stand and watch while everything settles.</summary>
        public bool IsSettled => World != null && World.IsSettled;
    }
}
