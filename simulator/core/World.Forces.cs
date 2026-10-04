using Unity.Mathematics;

namespace Tectonic
{
    public sealed partial class World
    {
        /// <summary>
        /// Sum per-cell forces into a single force and torque per plate. The host applies
        /// them — Rigidbody2D.AddForce / AddTorque in Unity, or the headless integrator here.
        /// The core never touches a transform.
        ///
        /// Three contributions:
        ///   basic drag     — resists absolute motion. In Unity most of this can instead be
        ///                    expressed as Rigidbody2D.linearDamping, since it is exactly
        ///                    proportional to velocity; it is kept explicit here so the
        ///                    headless run matches.
        ///   orogenic drag  — resists motion RELATIVE to a plate we are in contact with,
        ///                    amplified by how far the slab has descended. This is the only
        ///                    path by which geology feeds back into the rigid-body motion.
        ///   ridge push     — outward from the off-screen spreading ridges at the walls.
        ///                    The ambient drive that keeps the map evolving between visitors.
        /// </summary>
        void AccumulateForces()
        {
            float cs = Config.CellSize;
            float cellArea = HexLattice.CellArea(cs);

            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];

                // An anchored plate ignores forces, so computing them is wasted work. Zero
                // them so nothing downstream can accidentally apply a stale value.
                if (plate.IsStatic)
                {
                    plate.Force = float2.zero;
                    plate.Torque = 0f;
                    continue;
                }

                var cells = plate.Read;

                float2 totalForce = float2.zero;
                float totalTorque = 0f;
                float2 com = plate.LocalToWorld(plate.CenterOfMassLocal);

                for (int i = 0; i < cells.Length; i++)
                {
                    ref readonly Cell c = ref cells[i];
                    if (!c.Alive) continue;

                    float2 world = plate.CellWorld(plate.AxialAt(i), cs);
                    float2 vel = plate.VelocityAtWorld(world);

                    float2 f = -Config.DragCoefficient * cellArea * vel;

                    if (c.DraggingPlate >= 0 && c.DraggingPlate < Plates.Count)
                    {
                        var other = Plates[c.DraggingPlate];
                        float2 rel = vel - other.VelocityAtWorld(world);
                        float relLen = math.length(rel);
                        if (relLen > 1e-6f)
                        {
                            float slabPull = 1f + SubductionProgress(c) * 20f;
                            f -= (rel / relLen) * math.sqrt(relLen) * slabPull
                                 * Config.OrogenicDragCoefficient * cellArea;
                        }
                    }

                    f += RidgePushAt(world) * cellArea;

                    totalForce += f;
                    float2 r = world - com;
                    totalTorque += r.x * f.y - r.y * f.x;
                }

                plate.Force = totalForce;
                plate.Torque = totalTorque;
            }
        }

        /// <summary>
        /// Inward push from the walls, which stand in for spreading ridges just outside the
        /// region. Falls off with distance, so it drives boundary crust away from the edges
        /// without dominating the interior.
        /// </summary>
        public float2 RidgePushAt(float2 world)
        {
            float2 h = Config.WallHalfExtent;
            float range = math.max(Config.RidgePushRange, 1e-3f);
            float k = Config.RidgePush;
            float2 f = float2.zero;

            float dRight = h.x - world.x;
            if (dRight < range) f.x -= k * (1f - math.max(dRight, 0f) / range);

            float dLeft = world.x + h.x;
            if (dLeft < range) f.x += k * (1f - math.max(dLeft, 0f) / range);

            float dTop = h.y - world.y;
            if (dTop < range) f.y -= k * (1f - math.max(dTop, 0f) / range);

            float dBottom = world.y + h.y;
            if (dBottom < range) f.y += k * (1f - math.max(dBottom, 0f) / range);

            return f;
        }
    }
}
