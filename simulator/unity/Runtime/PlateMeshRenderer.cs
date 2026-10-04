using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace Tectonic.UnityBridge
{
    /// <summary>
    /// Top-down view: one mesh per plate, one hexagon per live cell, colour carried in vertex
    /// colours. Rebuilt only when the core reports new data (every geology tick, not every
    /// frame), so at 10 Hz this is cheap even at ten thousand cells.
    ///
    /// Mesh-per-plate rather than one texture for the whole map, because plates move
    /// independently — this way the GPU transform does the work and the mesh only changes
    /// when the crust does.
    /// </summary>
    [RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
    public sealed class PlateMeshRenderer : MonoBehaviour
    {
        PlateBody _body;
        TectonicRunner _runner;
        Mesh _mesh;

        /// <summary>One material for every plate — colour is per-vertex, so nothing about it is
        /// per-plate. Also collapses N material states into one for the renderer.</summary>
        static Material _sharedMaterial;
        static Material SharedMaterial =>
            _sharedMaterial != null ? _sharedMaterial : _sharedMaterial = new Material(Shader.Find("Sprites/Default"));

        readonly List<Vector3> _verts = new List<Vector3>(1 << 16);
        readonly List<Color32> _colors = new List<Color32>(1 << 16);
        readonly List<int> _tris = new List<int>(1 << 17);

        public void Bind(PlateBody body, TectonicRunner runner)
        {
            _body = body;
            _runner = runner;

            _mesh = new Mesh { indexFormat = UnityEngine.Rendering.IndexFormat.UInt32 };
            _mesh.MarkDynamic();
            GetComponent<MeshFilter>().sharedMesh = _mesh;

            var mr = GetComponent<MeshRenderer>();
            mr.sharedMaterial = SharedMaterial;
            mr.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
            mr.receiveShadows = false;

            Rebuild();
        }

        /// <summary>Meshes are not scene objects, so destroying the GameObject does not free
        /// them. Plates die and rift continuously, so without this the exhibit leaks native
        /// memory for as long as it runs.</summary>
        void OnDestroy()
        {
            if (_mesh != null) Destroy(_mesh);
        }

        void LateUpdate()
        {
            if (_body == null || !_body.Dirty) return;
            Rebuild();
            _body.ClearDirty();
        }

        void Rebuild()
        {
            if (_body?.Plate == null || _runner?.World == null) return;

            var world = _runner.World;
            var plate = _body.Plate;
            float cellSize = world.Config.CellSize;

            _verts.Clear(); _colors.Clear(); _tris.Clear();

            var cells = plate.Read;
            for (int i = 0; i < cells.Length; i++)
            {
                if (!cells[i].Alive) continue;

                int2 axial = plate.AxialAt(i);
                float2 c = HexLattice.AxialToLocal(axial, cellSize);
                float3 centreCol = CellPalette.ColorOf(world, cells[i]);

                int baseIdx = _verts.Count;
                _verts.Add(new Vector3(c.x, c.y, 0f));
                _colors.Add(ToColor32(centreCol));

                // Each corner vertex gets the average of the three cells meeting there, so the
                // fragment shader interpolates across cell boundaries and the hex grid stops
                // being visible. At 320 cells across on a 4K projector a flat-shaded hex is
                // ~12 px and reads as a mosaic; this dissolves it without simulating more.
                for (int k = 0; k < 6; k++)
                {
                    float2 offset = HexLattice.CornerOffset(k, cellSize);
                    _verts.Add(new Vector3(c.x + offset.x, c.y + offset.y, 0f));
                    _colors.Add(ToColor32(CellPalette.CornerColor(world, plate, axial, k, centreCol)));
                }

                for (int k = 0; k < 6; k++)
                {
                    _tris.Add(baseIdx);
                    _tris.Add(baseIdx + 1 + k);
                    _tris.Add(baseIdx + 1 + (k + 1) % 6);
                }
            }

            _mesh.Clear();
            _mesh.SetVertices(_verts);
            _mesh.SetColors(_colors);
            _mesh.SetTriangles(_tris, 0);
            _mesh.RecalculateBounds();
        }

        /// <summary>
        /// Same palette as the headless renderer, so what you verify offline is what the
        /// exhibit shows. Process overlays (subduction, volcanism, quakes) are drawn on top
        /// of the elevation ramp because the exhibit is about the process, not the outcome.
        /// </summary>
        /// <summary>Quantises a linear 0..1 colour to an opaque <see cref="Color32"/>. The
        /// palette itself lives in <see cref="CellPalette"/> so the Unity mesh and the headless
        /// PNG renderer cannot drift apart.</summary>
        static Color32 ToColor32(float3 c)
        {
            c = math.saturate(c) * 255f;
            return new Color32((byte)c.x, (byte)c.y, (byte)c.z, 255);
        }
    }
}
