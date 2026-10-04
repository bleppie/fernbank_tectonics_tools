using System.Collections.Generic;
using Unity.Mathematics;

namespace Tectonic
{
    public sealed partial class World
    {
        readonly List<int> _bfsQueue = new List<int>(1 << 16);

        /// <summary>
        /// Contact detection. No colliders and no spatial index: a cell's world position is
        /// inverse-transformed into a candidate plate's frame and rounded to an axial
        /// coordinate. O(1) per probe.
        ///
        /// This deliberately does NOT use Unity's physics. Plates must interpenetrate —
        /// subduction is one plate sliding under another — and a contact solver exists to
        /// prevent exactly that. Unity's solver is used only for the walls, where
        /// non-penetration is the correct behaviour.
        ///
        /// Results are transient and written into Read, complete before geology begins.
        /// </summary>
        void DetectCollisions()
        {
            float cs = Config.CellSize;

            for (int p = 0; p < Plates.Count; p++)
            {
                var cells = Plates[p].Read;
                for (int i = 0; i < cells.Length; i++)
                    if (cells[i].Alive) cells[i].ResetContact();
            }

            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];
                var cells = plate.Read;

                for (int i = 0; i < cells.Length; i++)
                {
                    if (!cells[i].Alive || cells[i].Frozen) continue;

                    float2 world = plate.CellWorld(plate.AxialAt(i), cs);

                    // Walk outward through the density-sorted plate list; the nearest-density
                    // plate occupying this point is the one we interact with.
                    for (int step = 1; step < Plates.Count; step++)
                    {
                        int belowIdx = p + step;
                        if (belowIdx < Plates.Count &&
                            Plates[belowIdx].WorldBoundsContain(world) &&
                            Plates[belowIdx].TryCellAtWorld(world, cs, out int bi, out _))
                        {
                            Contact(belowIdx, bi, p, i, world);
                            break;
                        }

                        int aboveIdx = p - step;
                        if (aboveIdx >= 0 &&
                            Plates[aboveIdx].WorldBoundsContain(world) &&
                            Plates[aboveIdx].TryCellAtWorld(world, cs, out int ai, out _))
                        {
                            Contact(p, i, aboveIdx, ai, world);
                            break;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// The subduct-vs-mountain-building decision, carried over from the JS model's
        /// fields-collision. Two stages: density picks who is underneath, then lithology
        /// decides what actually happens.
        /// </summary>
        void Contact(int bottomPlate, int bottomIdx, int topPlate, int topIdx, float2 world)
        {
            var pb = Plates[bottomPlate];
            var pt = Plates[topPlate];

            ref Cell bottom = ref pb.Read[bottomIdx];
            ref Cell top = ref pt.Read[topIdx];

            float relSpeed = math.length(pt.VelocityAtWorld(world) - pb.VelocityAtWorld(world));

            bottom.CollidingPlate = topPlate;
            top.CollidingPlate = bottomPlate;
            bottom.ContactRelSpeed = relSpeed;
            top.ContactRelSpeed = relSpeed;

            bool bottomCanSubduct = bottom.Rock.CanSubduct(Config.SubductionBlockingThickness);

            if (bottomCanSubduct)
            {
                // Ocean dives under whatever is above it, and feeds an arc volcano up top.
                bottom.Flags |= CellFlags.Subduct;
                top.Flags |= CellFlags.Volcanic;

                if (bottom.ContinentBuffer && !top.Rock.IsContinental && !top.ContinentBuffer)
                    Drag(ref bottom, ref top, bottomPlate, topPlate);
            }
            else
            {
                bool topCanSubduct = top.Rock.CanSubduct(Config.SubductionBlockingThickness);

                if (!topCanSubduct)
                {
                    // Neither column will go down: crust piles up instead.
                    bottom.Flags |= CellFlags.Orogeny;
                    top.Flags |= CellFlags.Orogeny;
                    Drag(ref bottom, ref top, bottomPlate, topPlate);
                }
                else if (top.ContinentBuffer)
                {
                    // Ocean caught between two closing continents is consumed so they can meet.
                    // Flagged in Read, not written to Write: this phase runs before the geology
                    // gather, which rebuilds the whole cell from Read and would discard it.
                    top.Flags |= CellFlags.Doomed;
                }
                else
                {
                    // A continent will not dive under an ocean; both just brake.
                    Drag(ref bottom, ref top, bottomPlate, topPlate);
                }
            }
        }

        static void Drag(ref Cell bottom, ref Cell top, int bottomPlate, int topPlate)
        {
            bottom.DraggingPlate = topPlate;
            top.DraggingPlate = bottomPlate;
        }

        /// <summary>
        /// Multi-source BFS marking oceanic cells within a few cells of continental crust.
        /// Cheap, and it is what lets a closing ocean basin behave differently from open sea
        /// floor. Uses a reusable list with a head cursor rather than a queue that shifts.
        /// </summary>
        void UpdateContinentBuffers()
        {
            int maxDist = Config.ContinentBufferCells;

            for (int p = 0; p < Plates.Count; p++)
            {
                var plate = Plates[p];
                var cells = plate.Read;

                _bfsQueue.Clear();

                for (int i = 0; i < cells.Length; i++)
                {
                    if (!cells[i].Alive) { continue; }
                    bool continental = cells[i].Rock.IsContinental;
                    cells[i].ContinentBuffer = continental;
                    if (continental) _bfsQueue.Add(i);
                }

                // Distance is tracked by processing the frontier level by level.
                int head = 0;
                for (int d = 0; d < maxDist && head < _bfsQueue.Count; d++)
                {
                    int levelEnd = _bfsQueue.Count;
                    while (head < levelEnd)
                    {
                        int idx = _bfsQueue[head++];
                        int2 a = plate.AxialAt(idx);
                        for (int dir = 0; dir < HexLattice.NeighborCount; dir++)
                        {
                            int ni = plate.Index(a + HexLattice.NeighborDirs[dir]);
                            if (ni < 0 || !cells[ni].Alive || cells[ni].ContinentBuffer) continue;
                            cells[ni].ContinentBuffer = true;
                            _bfsQueue.Add(ni);
                        }
                    }
                }

                // No need to mirror into Write — the gather rebuilds each cell from Read.
            }
        }
    }
}
