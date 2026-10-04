using System;

namespace Tectonic
{
    [Flags]
    public enum CellFlags : ushort
    {
        None = 0,
        /// <summary>Cell exists. A cell that dies writes Cell.Empty during the geology gather.</summary>
        Alive = 1 << 0,
        /// <summary>Cell touches the edge of its plate. Drives new-crust generation at divergent margins.</summary>
        Boundary = 1 << 1,
        /// <summary>Within continent-buffer distance of continental crust.</summary>
        ContinentBuffer = 1 << 2,
        /// <summary>Outside the visible region. Retained and transported, but geology is skipped.</summary>
        Frozen = 1 << 3,
        /// <summary>Contact says this cell should be diving. Transient, set by collision detection.</summary>
        Subduct = 1 << 5,
        /// <summary>Contact between two columns that both refuse to subduct: mountain building. Transient.</summary>
        Orogeny = 1 << 6,
        /// <summary>Overriding cell above a subducting one: arc volcanism. Transient.</summary>
        Volcanic = 1 << 7,
        /// <summary>
        /// Marked for death by a phase that runs BEFORE the geology sweep. Such a phase cannot
        /// write the Write buffer directly — the gather overwrites the whole cell from Read
        /// afterwards — so it flags the cell in Read and the gather honours it.
        /// </summary>
        Doomed = 1 << 8,
    }

    /// <summary>
    /// A crust column: one thickness per rock type, ordered top to bottom.
    ///
    /// The JS model stored a variable-length list of {rock, thickness} sorted by a fixed
    /// order index — but because it only ever kept one layer per rock type, the list is
    /// really a fixed 11-slot vector. Collapsing it removes all the list splicing, all the
    /// allocation, and makes the whole cell blittable.
    /// </summary>
    public unsafe struct RockColumn
    {
        fixed float _t[RockProps.LayerCount];

        /// <summary>Thickness of a layer. Index 0 is the top of the column.</summary>
        public float this[int layer]
        {
            get => _t[layer];
            set => _t[layer] = value < 0f ? 0f : value;
        }

        public float Total
        {
            get
            {
                float sum = 0f;
                for (int i = 0; i < RockProps.LayerCount; i++) sum += _t[i];
                return sum;
            }
        }

        public float SedimentThickness
        {
            get
            {
                float sum = 0f;
                for (int i = 0; i < RockProps.LayerCount; i++)
                    if (RockProps.IsSediment[i]) sum += _t[i];
                return sum;
            }
        }

        /// <summary>
        /// Thickness showing above zero elevation. Isostasy is faked with a single ratio:
        /// only part of a layer stands proud, the rest is the mountain root. Oceanic
        /// sediment is exempt because it sits on top of the column rather than in it.
        /// </summary>
        public float ThicknessAboveZero(float ratio)
        {
            float sum = 0f;
            int oceanicSed = RockProps.LayerOceanicSediment;
            for (int i = 0; i < RockProps.LayerCount; i++)
                sum += _t[i] * (i == oceanicSed ? 1f : ratio);
            return sum;
        }

        /// <summary>Bottom-most non-empty layer. Identity is positional: Gabbro floor means
        /// oceanic crust, Granite floor means continental. Returns -1 for an empty column.</summary>
        public int BottomLayer
        {
            get
            {
                for (int i = RockProps.LayerCount - 1; i >= 0; i--)
                    if (_t[i] > 0f) return i;
                return -1;
            }
        }

        public bool IsOceanic => BottomLayer == RockProps.LayerGabbro;
        public bool IsContinental => BottomLayer == RockProps.LayerGranite;

        /// <summary>
        /// Whether this column can be carried down a subduction zone: no non-subductable
        /// layer thicker than the blocking threshold. The entire lithological rule.
        /// </summary>
        public bool CanSubduct(float blockingThickness)
        {
            for (int i = 0; i < RockProps.LayerCount; i++)
                if (!RockProps.CanSubduct[i] && _t[i] > blockingThickness) return false;
            return true;
        }

        public void Clear()
        {
            for (int i = 0; i < RockProps.LayerCount; i++) _t[i] = 0f;
        }

        /// <summary>Scale every layer so the column reaches a target total thickness.</summary>
        public void ScaleTo(float targetTotal)
        {
            float total = Total;
            if (total <= 1e-6f) return;
            float k = targetTotal / total;
            for (int i = 0; i < RockProps.LayerCount; i++) _t[i] *= k;
        }
    }

    /// <summary>
    /// Per-cell simulation state. Blittable and fixed-size, so a plate's Read/Write buffers
    /// swap by reference each step and can later be dropped into a NativeArray for Burst/Jobs
    /// without changing shape.
    /// </summary>
    public struct Cell
    {
        public RockColumn Rock;

        /// <summary>Distance travelled since formation, not elapsed time — as in the JS model.</summary>
        public float Age;

        public float Metamorphic;

        /// <summary>Per-cell randomised ceiling on crust thickness.</summary>
        public float MaxCrustThickness;

        /// <summary>Per-cell randomised budget consumed by uplift; once spent, growth stops.</summary>
        public float UpliftCapacity;

        /// <summary>Trench depression. Propagates outward from subducting cells.</summary>
        public float BendingProgress;

        /// <summary>Distance travelled down the subduction zone. Negative means not subducting.</summary>
        public float SubductionDist;

        public float VolcanicIntensity;
        public float EruptionTime;

        public float EarthquakeMagnitude;
        public float EarthquakeDepth;
        public float EarthquakeLifespan;

        public CellFlags Flags;

        /// <summary>Plate index this cell is in contact with, or -1. Transient, recomputed each step.</summary>
        public int CollidingPlate;

        /// <summary>Plate index dragging this cell, or -1. The only channel from collisions
        /// back into the rigid-body physics. Transient.</summary>
        public int DraggingPlate;

        /// <summary>Closing speed at the contact, world units per time unit. Transient.</summary>
        public float ContactRelSpeed;

        public bool Alive
        {
            get => (Flags & CellFlags.Alive) != 0;
            set => Flags = value ? Flags | CellFlags.Alive : Flags & ~CellFlags.Alive;
        }

        public bool Frozen
        {
            get => (Flags & CellFlags.Frozen) != 0;
            set => Flags = value ? Flags | CellFlags.Frozen : Flags & ~CellFlags.Frozen;
        }

        public bool Boundary
        {
            get => (Flags & CellFlags.Boundary) != 0;
            set => Flags = value ? Flags | CellFlags.Boundary : Flags & ~CellFlags.Boundary;
        }

        public bool ContinentBuffer
        {
            get => (Flags & CellFlags.ContinentBuffer) != 0;
            set => Flags = value ? Flags | CellFlags.ContinentBuffer : Flags & ~CellFlags.ContinentBuffer;
        }

        public bool Subducting => SubductionDist >= 0f;

        /// <summary>Clear the per-step contact state. Persistent geology is untouched.</summary>
        public void ResetContact()
        {
            CollidingPlate = -1;
            DraggingPlate = -1;
            ContactRelSpeed = 0f;
            Flags &= ~(CellFlags.Subduct | CellFlags.Orogeny | CellFlags.Volcanic | CellFlags.Doomed);
        }

        public static Cell Empty => new Cell
        {
            Flags = CellFlags.None,
            SubductionDist = -1f,
            CollidingPlate = -1,
            DraggingPlate = -1,
        };
    }
}
