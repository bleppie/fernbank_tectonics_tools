using System.Collections.Generic;
using Tectonic;
using Unity.Mathematics;

namespace Tectonic.Headless
{
    /// <summary>
    /// Stands in for Unity. Does exactly what the Unity bridge will do:
    ///   1. mirror transforms into the core
    ///   2. Step()
    ///   3. read back Force/Torque and integrate
    ///   4. keep plates inside the walls
    ///
    /// In Unity, steps 3 and 4 are Rigidbody2D.AddForce/AddTorque and the physics solver
    /// resolving contacts against static wall colliders. Here they are twelve lines of
    /// semi-implicit Euler and an AABB clamp — deliberately simple, because the point of
    /// this host is to exercise the geology, not to be a physics engine.
    /// </summary>
    public sealed class HeadlessHost
    {
        public readonly World World;

        public HeadlessHost(World world) { World = world; }

        public void Step()
        {
            World.Step();

            float dt = World.Config.Timestep;

            for (int i = 0; i < World.Plates.Count; i++)
            {
                var p = World.Plates[i];
                if (p.IsStatic) { p.LinearVelocity = float2.zero; p.AngularVelocity = 0f; continue; }

                float2 acc = p.Force / p.Mass;
                float angAcc = p.Torque / p.MomentOfInertia;

                p.LinearVelocity += acc * dt;
                p.AngularVelocity += angAcc * dt;

                p.Position += p.LinearVelocity * dt;
                p.Rotation += p.AngularVelocity * dt;
            }

            ConstrainToWalls();
        }

        /// <summary>
        /// The hard walls. Unity's solver does this properly with contacts and impulses; the
        /// headless stand-in clamps the plate's world AABB and kills the offending velocity
        /// component, which is behaviourally close enough to validate the geology.
        /// </summary>
        void ConstrainToWalls()
        {
            float2 wall = World.Config.WallHalfExtent;
            float cs = World.Config.CellSize;

            for (int i = 0; i < World.Plates.Count; i++)
            {
                var p = World.Plates[i];
                if (p.IsStatic) continue;
                if (!TryWorldBounds(p, cs, out float2 min, out float2 max)) continue;

                float2 push = float2.zero;
                if (max.x > wall.x) push.x -= max.x - wall.x;
                if (min.x < -wall.x) push.x += -wall.x - min.x;
                if (max.y > wall.y) push.y -= max.y - wall.y;
                if (min.y < -wall.y) push.y += -wall.y - min.y;

                if (math.any(push != float2.zero))
                {
                    p.Position += push;
                    if (push.x != 0f) p.LinearVelocity.x = 0f;
                    if (push.y != 0f) p.LinearVelocity.y = 0f;
                }
            }
        }

        public static bool TryWorldBounds(Plate p, float cellSize, out float2 min, out float2 max)
        {
            min = new float2(float.MaxValue);
            max = new float2(float.MinValue);
            bool any = false;

            for (int i = 0; i < p.Read.Length; i++)
            {
                if (!p.Read[i].Alive) continue;
                float2 w = p.CellWorld(p.AxialAt(i), cellSize);
                min = math.min(min, w);
                max = math.max(max, w);
                any = true;
            }
            return any;
        }
    }
}
