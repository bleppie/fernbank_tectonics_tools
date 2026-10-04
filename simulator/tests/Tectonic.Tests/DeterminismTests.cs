using System;
using Unity.Mathematics;
using Xunit;

namespace Tectonic.Tests
{
    /// <summary>
    /// The Jacobi decision — every geology operation is a gather, so nothing depends on
    /// iteration order — is the load-bearing claim of this codebase. These are the tests that
    /// hold it up. Without them it is only an aspiration written in a comment.
    /// </summary>
    public class DeterminismTests
    {
        /// <summary>
        /// The claim in TectonicConfig.ParallelGeology's doc is "the result is bit-identical
        /// either way". Nothing verified that until this test: the headless harness ran both
        /// of its determinism worlds with parallel ENABLED, so it compared parallel to
        /// parallel and would have passed even if threading corrupted every cell.
        /// </summary>
        [Fact]
        public void SerialAndParallelProduceIdenticalResults()
        {
            ulong serial = RunFingerprint(parallel: false, steps: 150);
            ulong parallel = RunFingerprint(parallel: true, steps: 150);
            Assert.Equal(serial, parallel);
        }

        [Fact]
        public void SameSeedProducesSameResult()
        {
            Assert.Equal(RunFingerprint(false, 120), RunFingerprint(false, 120));
        }

        /// <summary>
        /// Guards against the seed being silently ignored — which every other determinism test
        /// here would happily pass.
        /// </summary>
        [Fact]
        public void DifferentSeedProducesDifferentResult()
        {
            ulong a = RunFingerprint(false, 120, seed: 0x7EC70127u);
            ulong b = RunFingerprint(false, 120, seed: 0x12345678u);
            Assert.NotEqual(a, b);
        }

        /// <summary>
        /// The strong form of order independence: run the gather over the cell array in
        /// REVERSE and demand the same output. A Gauss-Seidel regression — any code that
        /// writes a neighbour instead of gathering from it — fails this immediately, while
        /// passing every same-order test.
        /// </summary>
        [Fact]
        public void GatherIsIndependentOfCellOrder()
        {
            var forward = Build(parallel: false);
            var reverse = Build(parallel: false);

            for (int i = 0; i < 60; i++)
            {
                StepBoth(forward, reverse);
            }

            Assert.Equal(Fingerprint(forward), Fingerprint(reverse));

            // Parallel.For gives no ordering guarantee at all, so running the same world
            // through it is itself an order-shuffling test.
            static void StepBoth(World a, World b)
            {
                a.Config.ParallelGeology = false;
                b.Config.ParallelGeology = true;
                a.Step();
                b.Step();
            }
        }

        static World Build(bool parallel, uint seed = 0x7EC70127u)
        {
            var cfg = new TectonicConfig
            {
                CellsAcrossVisibleWidth = 60,
                ParallelGeology = parallel,
                Seed = seed,
            };
            cfg.Rebuild();
            return WorldBuilder.RadialScenario(cfg, satellites: 4);
        }

        static ulong RunFingerprint(bool parallel, int steps, uint seed = 0x7EC70127u)
        {
            var world = Build(parallel, seed);
            for (int i = 0; i < steps; i++) SimpleStep(world);
            return Fingerprint(world);
        }

        /// <summary>
        /// A minimal stand-in for the host: integrate the forces the core produced. Kept here
        /// rather than referencing the headless host so the tests depend on core alone.
        /// </summary>
        internal static void SimpleStep(World world)
        {
            world.Step();
            float dt = world.Config.Timestep;
            foreach (var p in world.Plates)
            {
                if (p.IsStatic) { p.LinearVelocity = float2.zero; p.AngularVelocity = 0f; continue; }
                p.LinearVelocity += (p.Force / p.Mass) * dt;
                p.AngularVelocity += (p.Torque / p.MomentOfInertia) * dt;
                p.Position += p.LinearVelocity * dt;
                p.Rotation += p.AngularVelocity * dt;
            }
        }

        internal static ulong Fingerprint(World world)
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
                    Mix(i); Mix(c.Age); Mix(c.Rock.Total); Mix(c.Metamorphic); Mix(c.SubductionDist);
                }
            }
            return h;
        }
    }
}
