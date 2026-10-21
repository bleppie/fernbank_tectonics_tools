using System;
using System.Collections.Generic;
using System.IO;
using UnityEditor;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>Editor window: pick an SVG and an output folder, tweak the settings, press Generate.</summary>
    public class PlateMeshWindow : EditorWindow
    {
        [SerializeField] string svgPath = "";
        [SerializeField] string outputFolder = "Assets/PlateMeshes";
        [SerializeField] PlateMeshSettings settings = new PlateMeshSettings();

        SerializedObject serializedWindow;

        [MenuItem("Tools/Plate Meshes/SVG To Plate Meshes...")]
        static void Open() => GetWindow<PlateMeshWindow>("SVG Plate Meshes");

        void OnEnable() => serializedWindow = new SerializedObject(this);

        void OnGUI()
        {
            DrawSvgPicker();
            DrawFolderPicker();

            serializedWindow.Update();
            EditorGUILayout.PropertyField(serializedWindow.FindProperty("settings"), true);
            serializedWindow.ApplyModifiedProperties();

            EditorGUILayout.Space();
            using (new EditorGUI.DisabledScope(!File.Exists(svgPath)))
            {
                if (GUILayout.Button("Generate Meshes", GUILayout.Height(28)))
                    Generate();
            }
        }

        // ================================================================ GUI pieces

        void DrawSvgPicker()
        {
            using (new EditorGUILayout.HorizontalScope())
            {
                svgPath = EditorGUILayout.TextField("SVG file", svgPath);
                if (GUILayout.Button("Browse", GUILayout.Width(60)))
                {
                    string chosen = EditorUtility.OpenFilePanel("Choose SVG", Application.dataPath, "svg");
                    if (!string.IsNullOrEmpty(chosen)) svgPath = chosen;
                }
            }
        }

        void DrawFolderPicker()
        {
            using (new EditorGUILayout.HorizontalScope())
            {
                outputFolder = EditorGUILayout.TextField("Output folder", outputFolder);
                if (GUILayout.Button("Browse", GUILayout.Width(60)))
                {
                    string chosen = EditorUtility.OpenFolderPanel("Choose output folder", "Assets", "");
                    if (!string.IsNullOrEmpty(chosen)) outputFolder = ToProjectRelativePath(chosen);
                }
            }
        }

        // ================================================================ generate + save

        void Generate()
        {
            try
            {
                string svgText = File.ReadAllText(svgPath);
                List<Mesh> meshes = PlateMeshPipeline.Build(svgText, settings);

                foreach (Mesh mesh in meshes)
                    SaveMeshAsset(mesh, outputFolder);

                AssetDatabase.SaveAssets();
                AssetDatabase.Refresh();
                Debug.Log($"Saved {meshes.Count} plate meshes to {outputFolder}.");
            }
            catch (Exception exception)
            {
                Debug.LogException(exception);
            }
        }

        /// <summary>Creates the asset, or overwrites an existing one in place so scene and prefab references survive.</summary>
        static void SaveMeshAsset(Mesh mesh, string folder)
        {
            EnsureFolderExists(folder);
            string path = $"{folder}/{mesh.name}.asset";

            var existing = AssetDatabase.LoadAssetAtPath<Mesh>(path);
            if (existing == null)
            {
                AssetDatabase.CreateAsset(mesh, path);
            }
            else
            {
                CopyMeshData(mesh, existing);
                DestroyImmediate(mesh);
            }
        }

        /// <summary>
        /// Replaces the contents of an existing mesh asset (keeping its GUID) by clearing it and re-assigning every stream.
        /// UV channels 0 to 2 are copied with whatever number of components they have (plates: Vector2 and Vector4,
        /// ribbons: three Vector2 sets).
        /// </summary>
        static void CopyMeshData(Mesh source, Mesh target)
        {
            // Read everything from the source first: the target is cleared next.
            var uvChannels = new List<Array>();
            for (int channel = 0; channel < 3; channel++)
                uvChannels.Add(ReadUvChannel(source, channel));

            target.Clear();
            target.indexFormat = source.indexFormat;
            target.vertices = source.vertices;
            target.normals = source.normals;
            for (int channel = 0; channel < uvChannels.Count; channel++)
                WriteUvChannel(target, channel, uvChannels[channel]);
            target.SetTriangles(source.triangles, 0);
            target.RecalculateBounds();
            EditorUtility.SetDirty(target);
        }

        static Array ReadUvChannel(Mesh mesh, int channel)
        {
            switch (mesh.GetVertexAttributeDimension(UnityEngine.Rendering.VertexAttribute.TexCoord0 + channel))
            {
                case 2: { var list = new List<Vector2>(); mesh.GetUVs(channel, list); return list.ToArray(); }
                case 3: { var list = new List<Vector3>(); mesh.GetUVs(channel, list); return list.ToArray(); }
                case 4: { var list = new List<Vector4>(); mesh.GetUVs(channel, list); return list.ToArray(); }
                default: return null;   // channel not present
            }
        }

        static void WriteUvChannel(Mesh mesh, int channel, Array uvs)
        {
            switch (uvs)
            {
                case Vector2[] v2: mesh.SetUVs(channel, new List<Vector2>(v2)); break;
                case Vector3[] v3: mesh.SetUVs(channel, new List<Vector3>(v3)); break;
                case Vector4[] v4: mesh.SetUVs(channel, new List<Vector4>(v4)); break;
            }
        }

        static void EnsureFolderExists(string folder)
        {
            if (AssetDatabase.IsValidFolder(folder)) return;

            string parent = Path.GetDirectoryName(folder)?.Replace('\\', '/');
            if (string.IsNullOrEmpty(parent) || !folder.StartsWith("Assets"))
                throw new ArgumentException($"Output folder must be inside the project's Assets folder: '{folder}'");

            EnsureFolderExists(parent);
            AssetDatabase.CreateFolder(parent, Path.GetFileName(folder));
        }

        static string ToProjectRelativePath(string absolutePath)
        {
            absolutePath = absolutePath.Replace('\\', '/');
            string dataPath = Application.dataPath;   // ".../Project/Assets"
            if (absolutePath.StartsWith(dataPath))
                return "Assets" + absolutePath.Substring(dataPath.Length);

            Debug.LogWarning("Output folder must be inside the project's Assets folder.");
            return "Assets/PlateMeshes";
        }
    }
}
