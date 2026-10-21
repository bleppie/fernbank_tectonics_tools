using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using Unity.Mathematics;
using UnityEngine;

namespace PlateMeshes
{
    /// <summary>
    /// Reads closed outlines from SVG text: every &lt;polygon&gt; and every sub-path of every &lt;path&gt;.
    /// Curves are flattened into straight segments. Coordinates are returned exactly as written in the file
    /// (SVG units, y pointing down); element transforms are not applied.
    /// Supported path commands: M L H V C S Q T Z (absolute and relative). Arcs are not supported.
    /// </summary>
    public static class SvgPolygonReader
    {
        // Matches either a command letter or a number such as "-.45", "12.", "1e-5".
        // Numbers may be glued together in SVG ("-.45-78.1" is two numbers), which this pattern handles.
        static readonly Regex TokenPattern =
            new Regex(@"[A-Za-z]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?", RegexOptions.Compiled);

        public static List<float2[]> Read(string svgText, float curveTolerance)
        {
            var outlines = new List<float2[]>();
            foreach (XElement element in XDocument.Parse(svgText).Descendants())
            {
                string tag = element.Name.LocalName;
                if ((tag == "polygon" || tag == "path") && element.Attribute("transform") != null)
                    Debug.LogWarning($"SVG <{tag}> has a transform attribute, which is ignored.");

                if (tag == "polygon")
                    outlines.Add(ParsePolygonPoints(element.Attribute("points")?.Value ?? ""));
                else if (tag == "path")
                    outlines.AddRange(ParsePathData(element.Attribute("d")?.Value ?? "", curveTolerance));
            }
            return outlines;
        }

        // ---------------------------------------------------------------- <polygon>

        static float2[] ParsePolygonPoints(string text)
        {
            var reader = new TokenReader(text);
            var points = new List<float2>();
            while (!reader.AtEnd)
                points.Add(reader.ReadPoint());
            return points.ToArray();
        }

        // ---------------------------------------------------------------- <path>

        static List<float2[]> ParsePathData(string data, float tolerance)
        {
            var subpaths = new List<float2[]>();
            var reader = new TokenReader(data);
            var current = new List<float2>();

            float2 position = float2.zero;
            float2 subpathStart = float2.zero;
            float2 lastCubicControl = float2.zero;      // needed by S (smooth cubic)
            float2 lastQuadraticControl = float2.zero;  // needed by T (smooth quadratic)
            bool hasCubicControl = false;
            bool hasQuadraticControl = false;
            char command = ' ';

            while (!reader.AtEnd)
            {
                if (reader.NextIsCommand) command = reader.ReadCommand();
                else if (command == ' ') throw new FormatException("SVG path data must start with a command.");

                bool relative = char.IsLower(command);
                float2 origin = relative ? position : float2.zero;
                bool previousWasCubic = hasCubicControl;
                bool previousWasQuadratic = hasQuadraticControl;
                hasCubicControl = false;
                hasQuadraticControl = false;

                switch (char.ToUpperInvariant(command))
                {
                    case 'M':
                    {
                        FinishSubpath(subpaths, current);
                        position = origin + reader.ReadPoint();
                        subpathStart = position;
                        current = new List<float2> { position };
                        command = relative ? 'l' : 'L'; // further coordinate pairs after M are implicit line-tos
                        break;
                    }
                    case 'L':
                    {
                        position = origin + reader.ReadPoint();
                        current.Add(position);
                        break;
                    }
                    case 'H':
                    {
                        float x = reader.ReadNumber();
                        position.x = relative ? position.x + x : x;
                        current.Add(position);
                        break;
                    }
                    case 'V':
                    {
                        float y = reader.ReadNumber();
                        position.y = relative ? position.y + y : y;
                        current.Add(position);
                        break;
                    }
                    case 'C':
                    {
                        float2 control1 = origin + reader.ReadPoint();
                        float2 control2 = origin + reader.ReadPoint();
                        float2 end = origin + reader.ReadPoint();
                        AddCubic(current, position, control1, control2, end, tolerance);
                        position = end;
                        lastCubicControl = control2;
                        hasCubicControl = true;
                        break;
                    }
                    case 'S':
                    {
                        float2 control1 = previousWasCubic ? 2f * position - lastCubicControl : position;
                        float2 control2 = origin + reader.ReadPoint();
                        float2 end = origin + reader.ReadPoint();
                        AddCubic(current, position, control1, control2, end, tolerance);
                        position = end;
                        lastCubicControl = control2;
                        hasCubicControl = true;
                        break;
                    }
                    case 'Q':
                    {
                        float2 control = origin + reader.ReadPoint();
                        float2 end = origin + reader.ReadPoint();
                        AddQuadratic(current, position, control, end, tolerance);
                        position = end;
                        lastQuadraticControl = control;
                        hasQuadraticControl = true;
                        break;
                    }
                    case 'T':
                    {
                        float2 control = previousWasQuadratic ? 2f * position - lastQuadraticControl : position;
                        float2 end = origin + reader.ReadPoint();
                        AddQuadratic(current, position, control, end, tolerance);
                        position = end;
                        lastQuadraticControl = control;
                        hasQuadraticControl = true;
                        break;
                    }
                    case 'Z':
                    {
                        position = subpathStart;
                        FinishSubpath(subpaths, current);
                        current = new List<float2>();
                        command = ' ';
                        break;
                    }
                    case 'A':
                        throw new NotSupportedException("SVG arc commands (A/a) are not supported. Convert arcs to curves first.");
                    default:
                        throw new FormatException($"Unknown SVG path command '{command}'.");
                }
            }

            FinishSubpath(subpaths, current);
            return subpaths;
        }

        static void FinishSubpath(List<float2[]> subpaths, List<float2> points)
        {
            if (points.Count >= 3)
                subpaths.Add(points.ToArray());
        }

        // ---------------------------------------------------------------- curve flattening

        static void AddCubic(List<float2> points, float2 p0, float2 p1, float2 p2, float2 p3, float tolerance)
        {
            float secondDifference = math.max(math.length(p0 - 2f * p1 + p2), math.length(p1 - 2f * p2 + p3));
            int segments = SegmentCount(3, secondDifference, tolerance);
            for (int k = 1; k <= segments; k++)
            {
                float t = k / (float)segments;
                float s = 1f - t;
                points.Add(s * s * s * p0 + 3f * s * s * t * p1 + 3f * s * t * t * p2 + t * t * t * p3);
            }
        }

        static void AddQuadratic(List<float2> points, float2 p0, float2 p1, float2 p2, float tolerance)
        {
            float secondDifference = math.length(p0 - 2f * p1 + p2);
            int segments = SegmentCount(2, secondDifference, tolerance);
            for (int k = 1; k <= segments; k++)
            {
                float t = k / (float)segments;
                float s = 1f - t;
                points.Add(s * s * p0 + 2f * s * t * p1 + t * t * p2);
            }
        }

        /// <summary>Wang's formula: segments needed so a degree-d Bezier deviates less than the tolerance.</summary>
        static int SegmentCount(int degree, float maxSecondDifference, float tolerance)
        {
            float n = math.sqrt(degree * (degree - 1) / 8f * maxSecondDifference / tolerance);
            return math.max(1, (int)math.ceil(n));
        }

        // ---------------------------------------------------------------- tokenizer

        sealed class TokenReader
        {
            readonly MatchCollection tokens;
            int index;

            public TokenReader(string text) { tokens = TokenPattern.Matches(text); }

            public bool AtEnd => index >= tokens.Count;
            public bool NextIsCommand => char.IsLetter(tokens[index].Value[0]);
            public char ReadCommand() => tokens[index++].Value[0];

            public float ReadNumber()
            {
                if (AtEnd || NextIsCommand)
                    throw new FormatException("Unexpected end of coordinates in SVG data.");
                return float.Parse(tokens[index++].Value, NumberStyles.Float, CultureInfo.InvariantCulture);
            }

            public float2 ReadPoint()
            {
                float x = ReadNumber();
                float y = ReadNumber();
                return new float2(x, y);
            }
        }
    }
}
