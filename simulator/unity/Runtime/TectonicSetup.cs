using UnityEngine;

namespace Tectonic.UnityBridge
{
    /// <summary>
    /// Walls and the layer collision matrix.
    ///
    /// This is the crux of "let Unity do the physics" without letting Unity break the
    /// simulation. Two rules, set once:
    ///
    ///     Plates <-> Walls   ON    Unity solves these properly. Non-penetration against a
    ///                              static wall is exactly what a contact solver is for.
    ///     Plates <-> Plates  OFF   Plates MUST interpenetrate. Subduction is one plate
    ///                              sliding under another; a solver would spend every frame
    ///                              pushing them apart. Plate-vs-plate overlap is answered by
    ///                              the core in O(1) lattice arithmetic instead.
    ///
    /// Set these in Project Settings > Physics 2D > Layer Collision Matrix, or call
    /// ApplyLayerMatrix() at startup.
    /// </summary>
    public static class TectonicSetup
    {
        public const string PlateLayerName = "TectonicPlates";
        public const string WallLayerName = "TectonicWalls";

        public static int PlateLayer => LayerMask.NameToLayer(PlateLayerName);
        public static int WallLayer => LayerMask.NameToLayer(WallLayerName);

        public static bool ApplyLayerMatrix()
        {
            int plates = PlateLayer;
            int walls = WallLayer;
            if (plates < 0 || walls < 0)
            {
                Debug.LogError($"[Tectonic] Create layers '{PlateLayerName}' and '{WallLayerName}' " +
                               "in Project Settings > Tags and Layers. Without them plates would " +
                               "collide with each other and the contact solver would fight " +
                               "subduction every frame.");
                return false;
            }

            Physics2D.IgnoreLayerCollision(plates, plates, true);   // plates pass through each other
            Physics2D.IgnoreLayerCollision(plates, walls, false);   // plates stop at walls
            return true;
        }

        /// <summary>
        /// Four static edges just outside the visible region. The gap between the visible
        /// edge and the wall is the margin: crust out there is retained and transported but
        /// its geology is frozen, which is what makes a visitor's push reversible.
        /// </summary>
        public static void BuildWalls(Transform parent, TectonicConfig cfg)
        {
            if (!ApplyLayerMatrix()) return;

            var go = new GameObject("Walls");
            go.transform.SetParent(parent, false);
            if (WallLayer >= 0) go.layer = WallLayer;

            var rb = go.AddComponent<Rigidbody2D>();
            rb.bodyType = RigidbodyType2D.Static;

            float w = cfg.WallHalfExtent.x;
            float h = cfg.WallHalfExtent.y;
            const float thickness = 5f;

            AddWall(go, new Vector2(w + thickness * 0.5f, 0f), new Vector2(thickness, h * 2f + thickness * 2f));
            AddWall(go, new Vector2(-w - thickness * 0.5f, 0f), new Vector2(thickness, h * 2f + thickness * 2f));
            AddWall(go, new Vector2(0f, h + thickness * 0.5f), new Vector2(w * 2f + thickness * 2f, thickness));
            AddWall(go, new Vector2(0f, -h - thickness * 0.5f), new Vector2(w * 2f + thickness * 2f, thickness));
        }

        static void AddWall(GameObject parent, Vector2 offset, Vector2 size)
        {
            var box = parent.AddComponent<BoxCollider2D>();
            box.offset = offset;
            box.size = size;
        }
    }
}
