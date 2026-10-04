using System.Threading.Tasks;
using Unity.Mathematics;

namespace Tectonic
{
    public sealed partial class World
    {
        // The geology sweep — the part that actually makes rock.
        //
        // THE CENTRAL DISCIPLINE: every operation here is a GATHER. A cell computes its own
        // new state by reading its own and its neighbours' OLD state, and writes only
        // itself. No cell ever writes into a neighbour.
        //
        // The JS model does the opposite: it is a Gauss-Seidel sweep where each cell pushes
        // sediment and metamorphism into its neighbours, so cells processed later in a step
        // see earlier cells' writes. That makes results depend on iteration order (the JS
        // source admits this causes visible divergence over long runs) and makes the loop
        // fundamentally serial.
        //
        // Reformulating scatter as gather is what buys us: order independence, reproducible
        // output, and a loop that can be handed to Parallel.For or Burst later without
        // changing a line of the physics.
        //
        // Mass conservation is preserved by making every transfer symmetric: outflow from a
        // cell and the matching inflow to its neighbour are both computed from the same old
        // state, so they agree exactly.
        //
        // Three passes, all parallel-safe:
        //   A1  elevation per cell
        //   A2  downhill gradient and how much material each cell sheds
        //   B   the gather, producing the new state
        /// <summary>Thickness ceiling for a column, which depends on whether it is oceanic.</summary>
        float ColumnCeiling(in Cell c, bool oceanic) =>
            oceanic ? math.min(c.MaxCrustThickness, Config.MaxOceanicThickness) : c.MaxCrustThickness;

        /// <summary>Below this a crust column no longer counts as crust.</summary>
        public const float MinViableThickness = 0.02f;

        void RunGeology(float dt)
        {
            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];
                ComputeElevations(plate);
                ComputeErosionBudget(plate, dt);
                GatherStep(plate, dt);
            }
        }

        // ---- Pass A1 ----------------------------------------------------------------

        void ComputeElevations(Plate plate)
        {
            var cells = plate.Read;
            var elev = plate.ScratchElevation;
            for (int i = 0; i < cells.Length; i++)
                elev[i] = cells[i].Alive ? Elevation(cells[i]) : float.NegativeInfinity;
        }

        // ---- Pass A2 ----------------------------------------------------------------

        void ComputeErosionBudget(Plate plate, float dt)
        {
            var cells = plate.Read;
            var elev = plate.ScratchElevation;
            var drop = plate.ScratchTotalDrop;
            var erodible = plate.ScratchErodible;

            for (int i = 0; i < cells.Length; i++)
            {
                drop[i] = 0f;
                erodible[i] = 0f;
                if (!cells[i].Alive || cells[i].Frozen) continue;

                int2 a = plate.AxialAt(i);
                float e = elev[i];
                float total = 0f;

                for (int d = 0; d < HexLattice.NeighborCount; d++)
                {
                    int ni = plate.Index(a + HexLattice.NeighborDirs[d]);
                    // Frozen neighbours are excluded deliberately. They are skipped by the
                    // gather, so they never accept a deposit — counting them as downhill sinks
                    // would route material to them and silently destroy it, breaking the
                    // conservation guarantee at the one-cell frozen rim.
                    if (ni < 0 || !cells[ni].Alive || cells[ni].Frozen) continue;
                    float diff = e - elev[ni];
                    if (diff > 0f) total += diff;
                }

                drop[i] = total;
                if (total <= 0f) continue;

                float loose = cells[i].Rock.SedimentThickness + Weathering(cells[i], e, dt);
                erodible[i] = math.min(loose, Config.ErosionCoefficient * dt * total);
            }
        }

        /// <summary>Bedrock converted to loose sediment in place, driven by height above sea level.</summary>
        float Weathering(in Cell c, float elevation, float dt) =>
            math.max(0f, elevation - Config.SeaLevel) * Config.WeatheringRate * dt;

        // ---- Pass B: the gather -----------------------------------------------------

        void GatherStep(Plate plate, float dt)
        {
            var read = plate.Read;
            var write = plate.Write;
            var elev = plate.ScratchElevation;
            var drop = plate.ScratchTotalDrop;
            var erodible = plate.ScratchErodible;

            float cs = Config.CellSize;
            float maxSubductionDist = Config.SubductionMaxDist;

            if (Config.ParallelGeology)
                Parallel.For(0, read.Length,
                    i => GatherCell(plate, read, write, elev, drop, erodible, i, cs, maxSubductionDist, dt));
            else
                for (int i = 0; i < read.Length; i++)
                    GatherCell(plate, read, write, elev, drop, erodible, i, cs, maxSubductionDist, dt);
        }

        /// <summary>
        /// One cell's update. Reads Read + scratch, writes only Write[i] — which is exactly
        /// why it is handed to Parallel.For today, and could go to an IJobParallelFor later,
        /// without touching a line of the physics.
        /// </summary>
        void GatherCell(Plate plate, Cell[] read, Cell[] write, float[] elev, float[] drop,
                        float[] erodible, int i, float cs, float maxSubductionDist, float dt)
        {
            // EVERY cell writes itself, so the Write buffer never needs seeding from Read.
            // Dead and frozen cells pass through unchanged.
            if (!read[i].Alive) { write[i] = Cell.Empty; return; }
            if (read[i].Frozen) { write[i] = read[i]; return; }

            Cell next = read[i];            // local copy; mutated then stored once

            // Killed by an earlier phase (e.g. ocean squeezed between closing continents).
            if ((read[i].Flags & CellFlags.Doomed) != 0) { write[i] = Cell.Empty; return; }
            int2 a = plate.AxialAt(i);
            float e = elev[i];
            var rng = CellRandom(plate.Id, i);

            // -- age is distance travelled, as in the JS model ---------------------
            float2 world = plate.CellWorld(a, cs);
            next.Age += math.max(math.length(plate.VelocityAtWorld(world)), Config.MinAgeRate) * dt;
            next.Age = math.min(next.Age, Config.FreshCrustMaxAge * Config.MatureCrustAge);

            // -- weathering: bedrock becomes loose sediment in place ---------------
            float weathered = Weathering(read[i], e, dt);
            if (weathered > 0f) WeatherBedrock(ref next, weathered);

            // -- erosion: shed material downhill ----------------------------------
            if (erodible[i] > 0f) RemoveSediment(ref next, erodible[i]);

            // -- deposition: gather what neighbours shed toward us -----------------
            float deposited = 0f;
            for (int d = 0; d < HexLattice.NeighborCount; d++)
            {
                int ni = plate.Index(a + HexLattice.NeighborDirs[d]);
                if (ni < 0 || !read[ni].Alive || drop[ni] <= 0f || erodible[ni] <= 0f) continue;
                float share = elev[ni] - e;
                if (share <= 0f) continue;
                deposited += erodible[ni] * (share / drop[ni]);
            }
            if (deposited > 0f) AddSediment(ref next, deposited, read[i].Rock.IsContinental);

            // -- subduction --------------------------------------------------------
            if ((read[i].Flags & CellFlags.Subduct) != 0)
            {
                float dist = math.max(read[i].SubductionDist, 0f) + read[i].ContactRelSpeed * dt;
                next.SubductionDist = dist;
                if (dist > maxSubductionDist) { write[i] = Cell.Empty; return; }

                // The slab sheds its sediment cover into the accretionary wedge as it dives.
                RemoveSediment(ref next, Config.SubductionSedimentRate * dt);
            }
            else
            {
                next.SubductionDist = -1f;
            }

            // -- trench bending, gathered from neighbours -------------------------
            float ownBend = SubductionProgress(read[i]) > 0f
                ? SubductionProgress(read[i]) + Config.TrenchMaxDepth
                : 0f;
            float bend = ownBend;
            for (int d = 0; d < HexLattice.NeighborCount; d++)
            {
                int ni = plate.Index(a + HexLattice.NeighborDirs[d]);
                if (ni < 0 || !read[ni].Alive) continue;
                bend = math.max(bend, read[ni].BendingProgress - Config.BendingFalloffPerCell);
            }
            next.BendingProgress = math.max(0f, bend);

            // -- young sea floor thickens as it ages ------------------------------
            if (read[i].Rock.IsOceanic)
            {
                float na = NormalizedAge(read[i]);
                if (na < 1f)
                {
                    float target = Config.BaseOceanicCrustThickness * (0.6f + 0.4f * na);
                    float current = next.Rock.Total;
                    if (current < target)
                    {
                        float add = math.min(target - current, Config.BasaltGabbroAccretionRate * dt);
                        next.Rock[RockProps.LayerBasalt] += add * 0.3f;
                        next.Rock[RockProps.LayerGabbro] += add * 0.7f;
                    }
                }
                else if (na > 0.7f)
                {
                    // Old sea floor accumulates a thin pelagic drape.
                    float sed = next.Rock[RockProps.LayerOceanicSediment];
                    if (sed < Config.MaxRegularSedimentThickness)
                        next.Rock[RockProps.LayerOceanicSediment] =
                            math.min(Config.MaxRegularSedimentThickness, sed + 0.02f * dt);
                }
            }

            // -- arc volcanism on the overriding plate ----------------------------
            if ((read[i].Flags & CellFlags.Volcanic) != 0 && next.UpliftCapacity > 0f)
            {
                float intensity = 0.35f + 0.65f * rng.NextFloat();
                next.VolcanicIntensity = intensity;
                float add = Config.VolcanicRockRate * intensity * dt;

                if (read[i].Rock.IsContinental)
                {
                    next.Rock[RockProps.LayerRhyolite] += add * 0.4f;
                    next.Rock[RockProps.LayerDiorite] += add * 0.6f;
                }
                else
                {
                    next.Rock[RockProps.LayerAndesite] += add * 0.5f;
                    next.Rock[RockProps.LayerDiorite] += add * 0.5f;
                }

                next.UpliftCapacity = math.max(0f, next.UpliftCapacity - dt);
                next.EruptionTime += dt;
            }
            else
            {
                next.VolcanicIntensity = math.max(0f, next.VolcanicIntensity - dt);
            }

            // -- mountain building where neither column will subduct --------------
            if ((read[i].Flags & CellFlags.Orogeny) != 0)
            {
                // Oceanic columns get a much lower ceiling than continental ones. Ocean
                // crust jammed between converging plates should stay ocean, not inflate
                // into a mountain range.
                float ceiling = ColumnCeiling(next, read[i].Rock.IsOceanic);

                float fold = Config.FoldRate * read[i].ContactRelSpeed * dt;
                if (fold > 0f && next.Rock.Total < ceiling)
                {
                    next.Rock.ScaleTo(math.min(ceiling, next.Rock.Total + fold));
                    next.Metamorphic = math.max(next.Metamorphic, Metamorphism.MediumGrade);
                }
            }

            // -- metamorphism spreads outward and never reverses ------------------
            float meta = next.Metamorphic;
            for (int d = 0; d < HexLattice.NeighborCount; d++)
            {
                int ni = plate.Index(a + HexLattice.NeighborDirs[d]);
                if (ni < 0 || !read[ni].Alive) continue;
                meta = math.max(meta, read[ni].Metamorphic * Config.MetamorphismDecayPerCell);
            }
            next.Metamorphic = math.saturate(meta);

            // -- earthquakes: short-lived markers near active subduction -----------
            if (next.EarthquakeLifespan > 0f)
            {
                next.EarthquakeLifespan -= dt;
            }
            else if ((read[i].Flags & CellFlags.Subduct) != 0 && rng.NextFloat() < Config.EarthquakeProbPerCell)
            {
                next.EarthquakeMagnitude = 1f + rng.NextFloat() * 8f;
                next.EarthquakeDepth = e - SubductionProgress(read[i]) * 0.6f;
                next.EarthquakeLifespan = 0.3f;
            }
            else
            {
                next.EarthquakeMagnitude = 0f;
            }

            // -- keep the column inside its per-cell ceiling ----------------------
            // Oceanic columns get a lower ceiling than continental ones, applied to the
            // total from ALL processes — folding, arc volcanics and deposition together.
            // Capping only the folding left arc volcanism free to lift whole ocean plates
            // above sea level.
            float ceilingFinal = ColumnCeiling(next, read[i].Rock.IsOceanic);
            float total = next.Rock.Total;
            if (total > ceilingFinal && total > 0f) next.Rock.ScaleTo(ceilingFinal);

            // A column worn away to nothing has been fully consumed. Retire the cell:
            // a live cell with no rock has no bottom layer, so it is neither oceanic nor
            // continental and every lithology test downstream becomes meaningless.
            if (total < MinViableThickness) { write[i] = Cell.Empty; return; }

            write[i] = next;
        }

        // ---- Column helpers ---------------------------------------------------------

        /// <summary>Convert the topmost bedrock layer into loose sediment, in place.</summary>
        static void WeatherBedrock(ref Cell c, float amount)
        {
            bool continental = c.Rock.IsContinental;
            for (int layer = 0; layer < RockProps.LayerCount && amount > 0f; layer++)
            {
                if (RockProps.IsSediment[layer]) continue;
                float take = math.min(c.Rock[layer], amount);
                if (take <= 0f) continue;
                c.Rock[layer] -= take;
                amount -= take;
                AddSediment(ref c, take, continental);
                break; // only the exposed layer weathers
            }
        }

        static readonly int[] SedimentOrder =
        {
            RockProps.LayerContinentalSediment,
            RockProps.LayerOceanicSediment,
        };

        static void RemoveSediment(ref Cell c, float amount)
        {
            var order = SedimentOrder;
            for (int k = 0; k < order.Length && amount > 0f; k++)
            {
                int layer = order[k];
                float take = math.min(c.Rock[layer], amount);
                c.Rock[layer] -= take;
                amount -= take;
            }
        }

        static void AddSediment(ref Cell c, float amount, bool continental)
        {
            int layer = continental ? RockProps.LayerContinentalSediment : RockProps.LayerOceanicSediment;
            c.Rock[layer] += amount;
        }

        // The old SweepDeadCells pass is gone. It walked every cell of every plate purely to
        // reset the ones that died this step — a whole extra pass that only existed because
        // the gather did not own the Write buffer end to end. Now that it does, a cell that
        // dies simply writes Cell.Empty and there is nothing left to sweep.
    }
}
