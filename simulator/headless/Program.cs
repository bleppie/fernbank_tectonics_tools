using System;
using System.Diagnostics;
using System.Globalization;
using Tectonic;
using Unity.Mathematics;

namespace Tectonic.Headless
{
    static class Program
    {
        const float PixelsPerUnit = 4f;

        /// <summary>Steps discarded before timing, so JIT and plate-rect growth settle.</summary>
        const int WarmupSteps = 50;

        /// <summary>Trailing window used for the steady-state allocation figure.</summary>
        const int SteadyStateWindow = 200;

        /// <summary>Allocation budget scales with plate count: the residual is plate-rect
        /// resizing, so more growing plates legitimately means more of it.</summary>
        const int AllocBudgetBytesPerTwoPlates = 1024;

        /// <summary>Margin cells are dimmed so the reversibility buffer is visible.</summary>
        const float MarginDimFactor = 0.45f;

        /// <summary>Roughly how many steps plates need to stop expanding into the margin.
        /// Before this, allocation is dominated by rect resizing and says nothing about
        /// steady-state behaviour.</summary>
        const int StepsToEquilibrium = 2000;

        static int Main(string[] args)
        {
            if (args.Length > 0 && args[0] == "--help") { PrintUsage(); return 0; }
            if (args.Length > 0 && args[0] == "bench") { Bench(args); return 0; }
            if (args.Length > 0 && args[0] == "benchplates") { BenchPlates(args); return 0; }

            var opt = Options.Parse(args);

            var cfg = new TectonicConfig
            {
                CellsAcrossVisibleWidth = opt.Resolution,
                ParallelGeology = opt.Parallel,
                Seed = opt.Seed,
            };
            cfg.Rebuild();

            System.IO.Directory.CreateDirectory(opt.OutDir);

            var world = opt.Scenario == "radial"
                ? WorldBuilder.RadialScenario(cfg, opt.Satellites)
                : WorldBuilder.DefaultScenario(cfg);
            var host = new HeadlessHost(world);
            int steps = opt.Steps;
            string outDir = opt.OutDir;

            Console.WriteLine("Tectonic core — headless run");
            Console.WriteLine($"  scenario         {opt.Scenario}");
            Console.WriteLine($"  resolution       {cfg.CellsAcrossVisibleWidth} cells across visible width");
            Console.WriteLine($"  cell size        {cfg.CellSize:F4} world units (spacing {cfg.CellSpacing:F4})");
            Console.WriteLine($"  parallel geology {cfg.ParallelGeology}");
            if (cfg.ErosionClampedForStability)
                Console.WriteLine("  WARNING          erosion coefficient clamped for stability — reduce Timestep at this resolution");
            Console.WriteLine($"  visible region   {cfg.VisibleHalfExtent.x * 2} x {cfg.VisibleHalfExtent.y * 2} world units");
            Console.WriteLine($"  wall extent      {cfg.WallHalfExtent.x * 2} x {cfg.WallHalfExtent.y * 2}");
            Console.WriteLine($"  timestep         {cfg.Timestep} ({cfg.Timestep * TectonicConfig.ModelTimeToMillionYears} Myr/step)");
            int staticPlates = 0;
            foreach (var pl in world.Plates) if (pl.IsStatic) staticPlates++;
            Console.WriteLine($"  plates           {world.Plates.Count} ({staticPlates} anchored)");
            Console.WriteLine();
            Console.WriteLine("  step     Myr  cells  live/frozen   subd  volc  new   elev[min..max]  relief   thickness   alloc");
            Console.WriteLine("  ------------------------------------------------------------------------------------------");

            var sw = Stopwatch.StartNew();
            long simAllocTotal = 0;
            long steadyAlloc = 0; int steadySamples = 0;
            double stepMsTotal = 0;
            int sampled = 0;

            Render(world, System.IO.Path.Combine(outDir, "map-0000.png"));

            for (int s = 1; s <= steps; s++)
            {
                long before = GC.GetTotalAllocatedBytes(false);
                var t0 = Stopwatch.GetTimestamp();
                host.Step();
                var t1 = Stopwatch.GetTimestamp();
                long stepAlloc = GC.GetTotalAllocatedBytes(false) - before;

                if (s > WarmupSteps)
                {
                    stepMsTotal += (t1 - t0) * 1000.0 / Stopwatch.Frequency;
                    simAllocTotal += stepAlloc;
                    sampled++;
                    if (s > steps - SteadyStateWindow) { steadyAlloc += stepAlloc; steadySamples++; }
                }

                if (s % 150 == 0 || s == steps)
                {
                    var st = Stats.Collect(world);
                    long allocPerStep = sampled > 0 ? simAllocTotal / sampled : 0;
                    Console.WriteLine(
                        $"  {s,4}  {world.TimeInMillionYears,6:F1}  {st.Total,5}  {st.Live,5}/{st.Frozen,-5}  " +
                        $"{st.Subducting,4}  {st.Volcanic,4}  {world.CellsCreatedTotal,4}new  " +
                        $"[{st.MinElev,5:F2}..{st.MaxElev,5:F2}]  {st.MaxElev - st.MinElev,6:F2}  " +
                        $"{st.TotalThickness,9:F1}   {allocPerStep,6} B");
                }

                if (s == steps / 2) Render(world, System.IO.Path.Combine(outDir, $"map-{s:0000}.png"));
            }

            sw.Stop();
            Render(world, System.IO.Path.Combine(outDir, $"map-{steps:0000}.png"));

            double avgStepMs = sampled > 0 ? stepMsTotal / sampled : 0;
            var final = Stats.Collect(world);

            Console.WriteLine();
            Console.WriteLine($"  wall clock       {sw.ElapsedMilliseconds} ms for {steps} steps");
            Console.WriteLine($"  mean step        {avgStepMs:F3} ms  ->  {1000.0 / Math.Max(avgStepMs, 1e-6):F0} steps/sec headroom");
            Console.WriteLine($"  simulated time   {world.TimeInMillionYears:F1} million years");
            Console.WriteLine($"  settled          {world.IsSettled} (static for {world.StaticSteps} steps)");
            Console.WriteLine();

            // --- Invariants ---------------------------------------------------------
            int failures = 0;
            failures += Check("no NaN elevations", !final.SawNaN);
            failures += Check("crust column never empty on a live cell", !final.SawEmptyColumn);
            failures += Check("cells exist", final.Live > 0);
            failures += Check("subduction occurred", final.EverSubducted);
            failures += Check("relief developed (mountains + trench)", final.MaxElev - final.MinElev > 0.55f);
            failures += Check("plates still inside walls", final.InsideWalls);
            long allocPerSimStep = sampled > 0 ? simAllocTotal / sampled : 0;
            long steadyPerStep = steadySamples > 0 ? steadyAlloc / steadySamples : 0;
            Console.WriteLine($"  sim allocation   {allocPerSimStep} B/step overall, {steadyPerStep} B/step steady-state (last 200)");
            Console.WriteLine($"  cells created    {world.CellsCreatedTotal} over the run");
            Console.WriteLine();
            // Budget scales with plate count: the residual is plate-rect resizing as plates
            // expand into the margin, so more growing plates legitimately allocate more.
            long allocBudget = AllocBudgetBytesPerTwoPlates * math.max(1, world.Plates.Count / 2);
            Console.WriteLine($"  alloc budget     {allocBudget} B/step ({world.Plates.Count} plates)");

            // Only meaningful once plates have reached their equilibrium extent. Before that,
            // rect resizing is the dominant cost and the figure is large by design — asserting
            // on a short run would just make this check flaky rather than informative.
            if (steps >= StepsToEquilibrium)
                failures += Check("steady-state step allocation is negligible", steadyPerStep < allocBudget);
            else
                Console.WriteLine($"  [SKIP] steady-state allocation — needs >= {StepsToEquilibrium} steps " +
                                  "for plate extents to settle; this run was too short to judge");

            // Determinism: same seed, same result, independent of run.
            var worldA = WorldBuilder.DefaultScenario(new TectonicConfig());
            var hostA = new HeadlessHost(worldA);
            for (int i = 0; i < 120; i++) hostA.Step();
            var worldB = WorldBuilder.DefaultScenario(new TectonicConfig());
            var hostB = new HeadlessHost(worldB);
            for (int i = 0; i < 120; i++) hostB.Step();
            failures += Check("deterministic across runs", Stats.Fingerprint(worldA) == Stats.Fingerprint(worldB));

            Console.WriteLine();
            Console.WriteLine(failures == 0 ? "  ALL CHECKS PASSED" : $"  {failures} CHECK(S) FAILED");
            return failures == 0 ? 0 : 1;
        }

        /// <summary>Physical width of the exhibit display, for the hex-size column.</summary>
        const double DisplayWidthInches = 144.0;

        /// <summary>
        /// Resolution sweep. Cost scales with the square of the resolution knob, so this is
        /// the table to choose from rather than guessing from a display size.
        /// </summary>
        static void Bench(string[] args)
        {
            int steps = args.Length > 1 ? int.Parse(args[1], CultureInfo.InvariantCulture) : 200;
            int[] resolutions = { 120, 180, 240, 320, 420, 560 };

            Console.WriteLine($"Resolution sweep — {steps} steps each, {Environment.ProcessorCount} logical cores");
            Console.WriteLine();
            Console.WriteLine($"  across   cell    hex@{DisplayWidthInches / 12.0:F0}ft     cells    serial ms   parallel ms   speedup   mem MB   10Hz load");
            Console.WriteLine("  ------------------------------------------------------------------------------------");

            foreach (int res in resolutions)
            {
                double serial = RunOnce(res, steps, false, out int cells, out double mem);
                double par = RunOnce(res, steps, true, out _, out _);
                double load = par * 10.0 / 1000.0 * 100.0;
                var probe = new TectonicConfig { CellsAcrossVisibleWidth = res };
                probe.Rebuild();   // object initialisers run after the ctor, so re-derive
                float cellSize = probe.CellSize;
                double inchesPerCell = DisplayWidthInches / res;

                Console.WriteLine(
                    $"  {res,6}   {cellSize,5:F3}   {inchesPerCell,6:F2}\"   {cells,7}   {serial,9:F2}   {par,11:F2}   " +
                    $"{serial / Math.Max(par, 1e-9),6:F2}x   {mem,6:F1}   {load,7:F1}%");
            }

            Console.WriteLine();
            Console.WriteLine("  'parallel ms' uses Parallel.For over the geology gather — legal only because");
            Console.WriteLine("  the sweep writes no neighbours. '10Hz load' is the fraction of ONE core needed");
            Console.WriteLine("  to sustain 10 geology steps/sec.");
        }

        /// <summary>
        /// Plate-count sweep at fixed resolution. Collision detection is the phase that cares
        /// about plate count, so this is the table to consult when a scenario's plate count is
        /// still undecided.
        /// </summary>
        static void BenchPlates(string[] args)
        {
            int steps = args.Length > 1 ? int.Parse(args[1], CultureInfo.InvariantCulture) : 120;
            int resolution = args.Length > 2 ? int.Parse(args[2], CultureInfo.InvariantCulture) : 240;
            int[] counts = { 2, 4, 6, 8, 12, 16 };

            Console.WriteLine($"Plate-count sweep — {steps} steps each at {resolution} cells across");
            Console.WriteLine();
            Console.WriteLine("  satellites   plates    cells    ms/step   10Hz load   5Hz load");
            Console.WriteLine("  ---------------------------------------------------------------");

            foreach (int sats in counts)
            {
                var cfg = new TectonicConfig { CellsAcrossVisibleWidth = resolution };
                cfg.Rebuild();

                var world = WorldBuilder.RadialScenario(cfg, sats);
                var host = new HeadlessHost(world);
                for (int i = 0; i < 20; i++) host.Step();

                var sw = Stopwatch.StartNew();
                for (int i = 0; i < steps; i++) host.Step();
                sw.Stop();

                int cells = 0;
                foreach (var p in world.Plates) foreach (var c in p.Read) if (c.Alive) cells++;

                double ms = sw.Elapsed.TotalMilliseconds / steps;
                Console.WriteLine($"  {sats,10}   {world.Plates.Count,6}   {cells,6}   {ms,7:F2}   " +
                                  $"{ms * 10 / 10.0,8:F1}%   {ms * 5 / 10.0,7:F1}%");
            }
        }

        static double RunOnce(int resolution, int steps, bool parallel, out int cells, out double memMB)
        {
            var cfg = new TectonicConfig { CellsAcrossVisibleWidth = resolution, ParallelGeology = parallel };
            cfg.Rebuild();

            var world = WorldBuilder.DefaultScenario(cfg);
            var host = new HeadlessHost(world);

            for (int i = 0; i < 20; i++) host.Step();   // warm up JIT and let rects settle

            var sw = Stopwatch.StartNew();
            for (int i = 0; i < steps; i++) host.Step();
            sw.Stop();

            cells = 0;
            long bytes = 0;
            foreach (var p in world.Plates)
            {
                foreach (var c in p.Read) if (c.Alive) cells++;
                bytes += (long)p.Read.Length * 2 * System.Runtime.InteropServices.Marshal.SizeOf<Cell>();
                bytes += (long)p.Read.Length * 3 * sizeof(float);
            }
            memMB = bytes / (1024.0 * 1024.0);
            return sw.Elapsed.TotalMilliseconds / steps;
        }

        /// <summary>
        /// Named flags rather than positional arguments. The positional form silently dropped
        /// anything past the second argument, so documented commands ran with the wrong
        /// resolution — or threw, in the case of `bench`.
        /// </summary>
        sealed class Options
        {
            public int Steps = 1200;
            public string OutDir = "out";
            public int Resolution = 320;
            public bool Parallel = false;   // see TectonicConfig.ParallelGeology
            public string Scenario = "radial";
            public int Satellites = 6;
            public uint Seed = 0x7EC70127u;

            public static Options Parse(string[] args)
            {
                var o = new Options();
                for (int i = 0; i < args.Length; i++)
                {
                    string key = args[i];
                    string Next() => i + 1 < args.Length ? args[++i]
                        : throw new ArgumentException($"{key} needs a value");

                    switch (key)
                    {
                        case "--steps":      o.Steps = int.Parse(Next(), CultureInfo.InvariantCulture); break;
                        case "--out":        o.OutDir = Next(); break;
                        case "--res":        o.Resolution = int.Parse(Next(), CultureInfo.InvariantCulture); break;
                        case "--scenario":   o.Scenario = Next(); break;
                        case "--satellites": o.Satellites = int.Parse(Next(), CultureInfo.InvariantCulture); break;
                        case "--seed":       o.Seed = uint.Parse(Next(), CultureInfo.InvariantCulture); break;
                        case "--parallel":   o.Parallel = true; break;
                        case "--serial":     o.Parallel = false; break;
                        default: throw new ArgumentException($"unknown argument '{key}' (try --help)");
                    }
                }
                return o;
            }
        }

        static void PrintUsage()
        {
            Console.WriteLine(@"Tectonic geology core — headless runner

  run-headless.sh [flags]          simulate, render PNG maps, check invariants
  run-headless.sh bench [steps]    resolution sweep
  run-headless.sh benchplates [steps] [res]   plate-count sweep

Flags:
  --steps N        simulation steps           (default 1200)
  --out DIR        where to write PNG maps    (default out)
  --res N          cells across visible width (default 320)
  --scenario NAME  radial | margin            (default radial)
  --satellites N   mobile plates, radial only (default 6)
  --seed N         RNG seed                   (default 2126250279)
  --serial         run the geology on one thread (default)
  --parallel       run the geology across threads; worth it above ~320 cells across");
        }

        static int Check(string label, bool ok)
        {
            Console.WriteLine($"  [{(ok ? "PASS" : "FAIL")}] {label}");
            return ok ? 0 : 1;
        }

        // ---- Rendering --------------------------------------------------------------

        static void Render(World world, string path)
        {
            var cfg = world.Config;
            float2 half = cfg.WallHalfExtent;
            int w = (int)(half.x * 2 * PixelsPerUnit);
            int h = (int)(half.y * 2 * PixelsPerUnit);
            var rgb = new byte[w * h * 3];

            for (int py = 0; py < h; py++)
            {
                for (int px = 0; px < w; px++)
                {
                    float2 worldPos = new float2(
                        -half.x + (px + 0.5f) / PixelsPerUnit,
                         half.y - (py + 0.5f) / PixelsPerUnit);

                    int plate = world.TopPlateAt(worldPos, out int cellIndex);
                    byte r, g, b;

                    if (plate < 0)
                    {
                        // No crust at all — the abyss beyond any plate.
                        r = 6; g = 8; b = 16;
                    }
                    else
                    {
                        float3 col = SampleSmooth(world, world.Plates[plate], cellIndex, worldPos);
                        if (!world.IsInsideVisible(worldPos)) col *= MarginDimFactor;
                        col = math.saturate(col) * 255f;
                        r = (byte)col.x; g = (byte)col.y; b = (byte)col.z;
                    }

                    int o = (py * w + px) * 3;
                    rgb[o] = r; rgb[o + 1] = g; rgb[o + 2] = b;
                }
            }

            // Visible-region frame, so the margin (the reversibility buffer) is obvious.
            DrawFrame(rgb, w, h, cfg.VisibleHalfExtent, cfg.WallHalfExtent);

            Png.Write(path, w, h, rgb);
            Console.WriteLine($"  wrote {path} ({w}x{h})");
        }

        /// <summary>
        /// Barycentric interpolation across the hexagon's triangle fan, with each corner
        /// carrying the average of the three cells that meet there. This mirrors exactly what
        /// the Unity mesh does with vertex colours, so what gets verified here is what the
        /// projector shows.
        ///
        /// The point: on a 4K projector at 320 cells across, one hexagon is ~12 px and flat
        /// shading reads as a mosaic. Interpolating dissolves the grid without simulating a
        /// single extra cell — visual resolution decoupled from simulation resolution.
        /// </summary>
        static float3 SampleSmooth(World world, Plate plate, int index, float2 worldPos)
        {
            float cs = world.Config.CellSize;
            int2 axial = plate.AxialAt(index);

            float2 local = plate.WorldToLocal(worldPos);
            float2 v = local - HexLattice.AxialToLocal(axial, cs);

            float ang = math.atan2(v.y, v.x);
            if (ang < 0f) ang += 2f * math.PI;
            int k = (int)(ang / (math.PI / 3f)) % 6;

            float2 c1 = HexLattice.CornerOffset(k, cs);
            float2 c2 = HexLattice.CornerOffset((k + 1) % 6, cs);

            float det = c1.x * c2.y - c1.y * c2.x;
            float wc1 = 0f, wc2 = 0f;
            if (math.abs(det) > 1e-9f)
            {
                wc1 = (v.x * c2.y - v.y * c2.x) / det;
                wc2 = (c1.x * v.y - c1.y * v.x) / det;
            }
            wc1 = math.saturate(wc1);
            wc2 = math.saturate(wc2);
            float wCentre = math.max(0f, 1f - wc1 - wc2);

            float total = wCentre + wc1 + wc2;
            if (total < 1e-6f) return CellPalette.ColorOf(world, plate.Read[index]);

            float3 centreCol = CellPalette.ColorOf(world, plate.Read[index]);
            float3 col = centreCol * wCentre
                       + CellPalette.CornerColor(world, plate, axial, k, centreCol) * wc1
                       + CellPalette.CornerColor(world, plate, axial, (k + 1) % 6, centreCol) * wc2;
            return col / total;
        }

        static void DrawFrame(byte[] rgb, int w, int h, float2 visible, float2 wall)
        {
            int x0 = (int)((wall.x - visible.x) * PixelsPerUnit);
            int x1 = w - x0 - 1;
            int y0 = (int)((wall.y - visible.y) * PixelsPerUnit);
            int y1 = h - y0 - 1;

            void Px(int x, int y)
            {
                if (x < 0 || y < 0 || x >= w || y >= h) return;
                int o = (y * w + x) * 3;
                rgb[o] = 255; rgb[o + 1] = 255; rgb[o + 2] = 255;
            }

            for (int x = x0; x <= x1; x += 2) { Px(x, y0); Px(x, y1); }
            for (int y = y0; y <= y1; y += 2) { Px(x0, y); Px(x1, y); }
        }
    }
}
