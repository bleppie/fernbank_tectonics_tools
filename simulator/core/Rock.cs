namespace Tectonic
{
    /// <summary>
    /// Rock types. Numeric values are carried over from the JS model unchanged so saved
    /// scenarios stay readable across both implementations; the JS source warns
    /// "Do not change numeric values without strong reason".
    /// </summary>
    public enum Rock : byte
    {
        OceanicSediment = 0,
        Granite = 1,
        Basalt = 2,
        Gabbro = 3,
        Rhyolite = 4,
        Andesite = 5,
        Diorite = 6,
        ContinentalSediment = 7,
        Limestone = 8,
        Shale = 9,
        Sandstone = 10,
    }

    /// <summary>
    /// Static rock property tables.
    ///
    /// The key structural finding from the JS analysis: a crust column holds AT MOST ONE
    /// LAYER PER ROCK TYPE, kept sorted by a fixed order index. That makes the whole
    /// variable-length `IRockLayer[]` collapse into a fixed 11-float array indexed by
    /// layer position — no lists, no allocation, no sorting, and trivially double-bufferable.
    ///
    /// Everything here is indexed by LAYER INDEX (0 = top of column), not by the Rock enum.
    /// Use <see cref="LayerOf"/> to convert.
    /// </summary>
    public static class RockProps
    {
        public const int LayerCount = 11;

        /// <summary>Layer index (0 = top of the column) to rock type.</summary>
        public static readonly Rock[] ByLayer =
        {
            Rock.ContinentalSediment, // 0  top
            Rock.OceanicSediment,     // 1
            Rock.Rhyolite,            // 2
            Rock.Andesite,            // 3
            Rock.Diorite,             // 4
            Rock.Sandstone,           // 5
            Rock.Shale,               // 6
            Rock.Limestone,           // 7
            Rock.Granite,             // 8
            Rock.Basalt,              // 9
            Rock.Gabbro,              // 10 bottom
        };

        /// <summary>Rock type to layer index. Indexed by (byte)Rock.</summary>
        public static readonly int[] LayerOf = BuildLayerOf();

        /// <summary>
        /// Whether a layer can be carried down a subduction zone. Only the three oceanic
        /// rocks can. This single predicate is the entire lithological rule that decides
        /// subduction vs. mountain-building.
        /// </summary>
        public static readonly bool[] CanSubduct = BuildFlags(
            Rock.OceanicSediment, Rock.Basalt, Rock.Gabbro);

        /// <summary>Sediment layers erode, transfer and spread; everything else is bedrock.</summary>
        public static readonly bool[] IsSediment = BuildFlags(
            Rock.OceanicSediment, Rock.ContinentalSediment);

        // Convenience layer indices for the hot paths.
        public static readonly int LayerContinentalSediment = LayerOf[(int)Rock.ContinentalSediment];
        public static readonly int LayerOceanicSediment     = LayerOf[(int)Rock.OceanicSediment];
        public static readonly int LayerGranite             = LayerOf[(int)Rock.Granite];
        public static readonly int LayerBasalt              = LayerOf[(int)Rock.Basalt];
        public static readonly int LayerGabbro              = LayerOf[(int)Rock.Gabbro];
        public static readonly int LayerAndesite            = LayerOf[(int)Rock.Andesite];
        public static readonly int LayerDiorite             = LayerOf[(int)Rock.Diorite];
        public static readonly int LayerRhyolite            = LayerOf[(int)Rock.Rhyolite];
        public static readonly int LayerLimestone           = LayerOf[(int)Rock.Limestone];
        public static readonly int LayerShale               = LayerOf[(int)Rock.Shale];
        public static readonly int LayerSandstone           = LayerOf[(int)Rock.Sandstone];

        static int[] BuildLayerOf()
        {
            var map = new int[LayerCount];
            for (int layer = 0; layer < LayerCount; layer++) map[(int)ByLayer[layer]] = layer;
            return map;
        }

        static bool[] BuildFlags(params Rock[] rocks)
        {
            var flags = new bool[LayerCount];
            foreach (var r in rocks) flags[LayerOf[(int)r]] = true;
            return flags;
        }

    }

    /// <summary>Degree of metamorphism, as a continuous 0..1 factor.</summary>
    public static class Metamorphism
    {
        public const float MediumGrade = 0.5f;
    }
}
