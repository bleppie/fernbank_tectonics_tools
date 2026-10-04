using System.Collections.Generic;
using Unity.Mathematics;

namespace Tectonic
{
    public sealed partial class World
    {
        struct NewCellRequest
        {
            public int PlateIndex;
            public int2 Axial;
            public Cell Source;
        }

        readonly List<NewCellRequest> _newCells = new List<NewCellRequest>(512);

        /// <summary>
        /// Crust creation at divergent boundaries.
        ///
        /// A plate edge facing genuinely open space grows into it. Two flavours, matching
        /// the JS model: if the neighbouring column is continental and thick enough to
        /// stretch, the new cell is thinned continental crust (rifting); otherwise it is
        /// fresh sea floor at age zero, which thickens over the following steps.
        ///
        /// This same path serves the walls. The walls stand in for spreading ridges just
        /// outside the region, so a plate retreating from a wall leaves new ocean behind
        /// rather than a void — no special case needed, and it gives the system a steady
        /// source of crust to pair with subduction as a sink.
        ///
        /// Creation is bounded by the walls, which is what keeps the cell count finite.
        /// </summary>
        void GenerateNewCrust(float dt)
        {
            float cs = Config.CellSize;
            _newCells.Clear();
            CellsCreatedLastStep = 0;

            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];
                var read = plate.Read;

                for (int i = 0; i < read.Length; i++)
                {
                    if (!read[i].Alive || !read[i].Boundary) continue;

                    int2 a = plate.AxialAt(i);

                    for (int d = 0; d < HexLattice.NeighborCount; d++)
                    {
                        int2 target = a + HexLattice.NeighborDirs[d];
                        int ti = plate.Index(target);
                        if (ti >= 0 && (read[ti].Alive || plate.Write[ti].Alive)) continue;

                        float2 world = plate.CellWorld(target, cs);
                        if (!IsInsideWalls(world)) continue;

                        // Only fill space no plate occupies — that is what makes this a
                        // divergent boundary rather than an overlap.
                        if (IsCoveredByAnyPlate(world, p)) continue;

                        // Require a couple of existing neighbours so plates grow as sheets
                        // rather than sprouting spindles into open water.
                        if (CountLiveNeighbors(plate, target) < Config.NewCrustMinNeighbors) continue;

                        _newCells.Add(new NewCellRequest { PlateIndex = p, Axial = target, Source = read[i] });
                    }
                }
            }

            // Applied after iteration: EnsureContains may reallocate the cell arrays.
            for (int r = 0; r < _newCells.Count; r++)
            {
                var req = _newCells[r];
                var plate = Plates[req.PlateIndex];
                plate.EnsureContains(req.Axial);

                int idx = plate.Index(req.Axial);
                if (idx < 0 || plate.Write[idx].Alive) continue;

                plate.Write[idx] = MakeNewCell(plate, req.Source, idx);
                CellsCreatedLastStep++;
                CellsCreatedTotal++;
            }
        }

        bool IsCoveredByAnyPlate(float2 world, int exceptPlate)
        {
            for (int q = 0; q < Plates.Count; q++)
            {
                if (q == exceptPlate) continue;
                if (!Plates[q].WorldBoundsContain(world)) continue;   // cheap reject
                if (Plates[q].TryCellAtWorld(world, Config.CellSize, out _, out _)) return true;
            }
            return false;
        }

        static int CountLiveNeighbors(Plate plate, int2 axial)
        {
            int count = 0;
            for (int d = 0; d < HexLattice.NeighborCount; d++)
            {
                int ni = plate.Index(axial + HexLattice.NeighborDirs[d]);
                if (ni >= 0 && plate.Read[ni].Alive) count++;
            }
            return count;
        }

        Cell MakeNewCell(Plate plate, in Cell src, int newIndex)
        {
            var rng = CellRandom(plate.Id, newIndex, 0x9E3779B9u);

            Cell c = Cell.Empty;
            c.Alive = true;
            c.MaxCrustThickness = Config.MaxCrustThicknessBase + rng.NextFloat();
            c.UpliftCapacity = 0.6f + 0.7f * rng.NextFloat();

            bool stretchable = src.Rock.IsContinental &&
                               src.Rock.Total - Config.ContinentalStretchPerCell > Config.MinContinentalCrustThickness;

            if (stretchable)
            {
                // Continental rifting: the new cell inherits a thinned version of its parent.
                c.Rock = src.Rock;
                c.Rock.ScaleTo(math.max(Config.MinContinentalCrustThickness,
                                        src.Rock.Total - Config.ContinentalStretchPerCell));
                c.Age = src.Age;
            }
            else
            {
                // Fresh sea floor, age zero, thin — it thickens over the next steps.
                float t = Config.BaseOceanicCrustThickness * 0.6f;
                c.Rock[RockProps.LayerBasalt] = t * 0.3f;
                c.Rock[RockProps.LayerGabbro] = t * 0.7f;
                c.Age = 0f;
            }

            return c;
        }
    }
}
