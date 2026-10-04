using System;
using Tectonic;
using Unity.Mathematics;

namespace Tectonic.Headless
{
    public struct Stats
    {
        public int Total, Live, Frozen, Subducting, Volcanic;
        public float MinElev, MaxElev, TotalThickness;
        public bool SawNaN, SawEmptyColumn, InsideWalls;
        public bool EverSubducted;

        public static Stats Collect(World world)
        {
            var s = new Stats
            {
                MinElev = float.MaxValue,
                MaxElev = float.MinValue,
                InsideWalls = true,
            };

            float cs = world.Config.CellSize;
            float2 wall = world.Config.WallHalfExtent + 1.5f;

            foreach (var p in world.Plates)
            {
                for (int i = 0; i < p.Read.Length; i++)
                {
                    ref readonly Cell c = ref p.Read[i];
                    if (!c.Alive) continue;

                    s.Total++;
                    if (c.Frozen) s.Frozen++; else s.Live++;
                    if (c.Subducting) s.Subducting++;
                    if (c.VolcanicIntensity > 0.01f) s.Volcanic++;

                    float e = world.Elevation(c);
                    if (float.IsNaN(e) || float.IsInfinity(e)) s.SawNaN = true;
                    else { s.MinElev = math.min(s.MinElev, e); s.MaxElev = math.max(s.MaxElev, e); }

                    float t = c.Rock.Total;
                    if (t <= 0f) s.SawEmptyColumn = true;
                    s.TotalThickness += t;

                    float2 wpos = p.CellWorld(p.AxialAt(i), cs);
                    if (math.any(math.abs(wpos) > wall)) s.InsideWalls = false;
                }
            }

            s.EverSubducted = s.Subducting > 0;
            return s;
        }

        /// <summary>Cheap order-sensitive hash of the whole world, for the determinism check.</summary>
        public static ulong Fingerprint(World world)
        {
            ulong h = 1469598103934665603UL;

            void Mix(float v)
            {
                uint bits = (uint)BitConverter.SingleToInt32Bits(v);
                h = (h ^ bits) * 1099511628211UL;
            }

            foreach (var p in world.Plates)
            {
                Mix(p.Position.x); Mix(p.Position.y); Mix(p.Rotation);
                for (int i = 0; i < p.Read.Length; i++)
                {
                    ref readonly Cell c = ref p.Read[i];
                    if (!c.Alive) continue;
                    Mix(i);
                    Mix(c.Age);
                    Mix(c.Rock.Total);
                    Mix(c.Metamorphic);
                    Mix(c.SubductionDist);
                }
            }
            return h;
        }
    }
}
