using Unity.Mathematics;
using UnityEngine;

namespace Tectonic.UnityBridge
{
    /// <summary>
    /// Inspector-facing mirror of TectonicConfig.
    ///
    /// The core deliberately has no ScriptableObject, no attributes and no UnityEngine
    /// reference — it is a plain C# class so it can be built and tested headlessly. This
    /// asset is the adapter, and it is the only place the two vocabularies meet.
    /// </summary>
    [CreateAssetMenu(menuName = "Tectonic/Settings", fileName = "TectonicSettings")]
    public sealed class TectonicSettings : ScriptableObject
    {
        [Header("Resolution — startup only")]
        [Tooltip("Cells spanning the visible width. THE resolution knob. Cost scales with " +
                 "the SQUARE: doubling this quadruples cell count and step time. On a 12 ft " +
                 "display: 120 -> 1.2\" hexes, 320 -> 0.45\", 420 -> 0.34\". See the " +
                 "benchmark table in the README.")]
        [Range(60, 640)] public int CellsAcrossVisibleWidth = 320;

        [Tooltip("Run the per-cell geology across threads. Safe by construction — the sweep " +
                 "is a pure gather, so results are bit-identical either way (asserted by " +
                 "DeterminismTests). Worth ~1.3x above roughly 320 cells across; below that it " +
                 "buys nothing in a full step and costs ~14 KB/step of allocation.")]
        public bool ParallelGeology = false;

        [Header("Scenario")]
        [Tooltip("The exhibit topology: an anchored central continent ringed by mobile plates. " +
                 "Off runs the two-plate convergent-margin test scenario instead.")]
        public bool RadialScenario = true;

        [Range(2, 16)] public int Satellites = 6;

        [Header("Extent")]
        [Tooltip("128 x 72 world units is exactly 16:9, matching a 3840x2160 projector.")]
        public Vector2 VisibleHalfExtent = new Vector2(64f, 36f);

        [Tooltip("Margin between the visible edge and the walls. Crust out here is retained " +
                 "and transported but its geology is frozen — this is what makes a visitor's " +
                 "push reversible, and it bounds the memory cost.")]
        public Vector2 MarginHalfExtent = new Vector2(9.6f, 5.4f);

        [Header("Time")]
        [Tooltip("Model time units per geology step. 1 unit = 3 million years.")]
        public float Timestep = 0.04f;

        [Header("Forces")]
        public float DragCoefficient = 0.006f;

        [Tooltip("Resists motion relative to a plate in contact. The only feedback path from " +
                 "geology into rigid-body motion — this is what makes a continental collision " +
                 "actually arrest the plates instead of letting them interpenetrate forever.")]
        public float OrogenicDragCoefficient = 0.03f;

        [Tooltip("Inward push from the walls, which stand in for spreading ridges just " +
                 "outside the region. The ambient drive that keeps the map evolving between " +
                 "visitors.")]
        public float RidgePush = 0.0022f;
        public float RidgePushRange = 24f;

        [Header("Spatial extents — WORLD UNITS, never cells")]
        [Tooltip("Authoring these in world units is what keeps the geology looking like the " +
                 "same planet at any resolution. Per-cell values are derived at startup.")]
        public float SubductionWidth = 19f;
        public float ContinentBufferWidth = 12f;
        public float BendingFalloffDistance = 7f;
        public float MetamorphismHaloDistance = 3f;

        [Header("Geological rates")]
        [Tooltip("Surface-transport diffusivity in world units squared — resolution-independent.")]
        public float ErosionDiffusivity = 0.66f;
        public float WeatheringRate = 0.022f;
        public float FoldRate = 1.1f;   // MUST match TectonicConfig; ConfigMirrorTests asserts it
        public float VolcanicRockRate = 0.9f;
        public float BasaltGabbroAccretionRate = 0.6f;

        [Header("Exhibit")]
        [Tooltip("Below this relative plate motion the map counts as settled — the cue to " +
                 "begin the slow reset.")]
        public float StaticMotionThreshold = 0.004f;
        public int StaticStepsBeforeIdle = 600;

        [Header("Determinism")]
        public uint Seed = 0x7EC70127u;

        public TectonicConfig ToConfig()
        {
            var cfg = new TectonicConfig
            {
                CellsAcrossVisibleWidth = CellsAcrossVisibleWidth,
                ParallelGeology = ParallelGeology,
                VisibleHalfExtent = new float2(VisibleHalfExtent.x, VisibleHalfExtent.y),
                MarginHalfExtent = new float2(MarginHalfExtent.x, MarginHalfExtent.y),
                Timestep = Timestep,
                DragCoefficient = DragCoefficient,
                OrogenicDragCoefficient = OrogenicDragCoefficient,
                RidgePush = RidgePush,
                RidgePushRange = RidgePushRange,
                ErosionDiffusivity = ErosionDiffusivity,
                WeatheringRate = WeatheringRate,
                FoldRate = FoldRate,
                VolcanicRockRate = VolcanicRockRate,
                BasaltGabbroAccretionRate = BasaltGabbroAccretionRate,
                SubductionWidth = SubductionWidth,
                ContinentBufferWidth = ContinentBufferWidth,
                BendingFalloffDistance = BendingFalloffDistance,
                MetamorphismHaloDistance = MetamorphismHaloDistance,
                StaticMotionThreshold = StaticMotionThreshold,
                StaticStepsBeforeIdle = StaticStepsBeforeIdle,
                Seed = Seed,
            };
            cfg.Rebuild();   // derive per-cell values from the world-unit constants
            return cfg;
        }
    }
}
