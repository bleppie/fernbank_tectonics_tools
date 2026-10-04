using System;
using System.IO;

namespace Tectonic.Headless
{
    /// <summary>
    /// Minimal PNG writer (8-bit RGB, no filtering, stored deflate blocks). Just enough to
    /// look at the simulation without pulling in an image library.
    /// </summary>
    public static class Png
    {
        public static void Write(string path, int width, int height, byte[] rgb)
        {
            using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
            fs.Write(new byte[] { 0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A }, 0, 8);

            var ihdr = new byte[13];
            WriteBE(ihdr, 0, width);
            WriteBE(ihdr, 4, height);
            ihdr[8] = 8;   // bit depth
            ihdr[9] = 2;   // colour type: truecolour
            Chunk(fs, "IHDR", ihdr);

            // Raw scanlines, each prefixed with filter type 0.
            var raw = new byte[height * (1 + width * 3)];
            int o = 0;
            for (int y = 0; y < height; y++)
            {
                raw[o++] = 0;
                Buffer.BlockCopy(rgb, y * width * 3, raw, o, width * 3);
                o += width * 3;
            }

            Chunk(fs, "IDAT", ZlibStored(raw));
            Chunk(fs, "IEND", Array.Empty<byte>());
        }

        static byte[] ZlibStored(byte[] data)
        {
            using var ms = new MemoryStream();
            ms.WriteByte(0x78);
            ms.WriteByte(0x01);

            int pos = 0;
            while (pos < data.Length)
            {
                int len = Math.Min(65535, data.Length - pos);
                bool last = pos + len >= data.Length;
                ms.WriteByte((byte)(last ? 1 : 0));
                ms.WriteByte((byte)(len & 0xFF));
                ms.WriteByte((byte)(len >> 8));
                ms.WriteByte((byte)(~len & 0xFF));
                ms.WriteByte((byte)((~len >> 8) & 0xFF));
                ms.Write(data, pos, len);
                pos += len;
            }

            uint a = Adler32(data);
            ms.WriteByte((byte)(a >> 24)); ms.WriteByte((byte)(a >> 16));
            ms.WriteByte((byte)(a >> 8)); ms.WriteByte((byte)a);
            return ms.ToArray();
        }

        static void Chunk(Stream s, string type, byte[] data)
        {
            var len = new byte[4];
            WriteBE(len, 0, data.Length);
            s.Write(len, 0, 4);

            var full = new byte[4 + data.Length];
            for (int i = 0; i < 4; i++) full[i] = (byte)type[i];
            Buffer.BlockCopy(data, 0, full, 4, data.Length);
            s.Write(full, 0, full.Length);

            var crc = new byte[4];
            WriteBE(crc, 0, unchecked((int)Crc32(full)));
            s.Write(crc, 0, 4);
        }

        static void WriteBE(byte[] b, int o, int v)
        {
            b[o] = (byte)(v >> 24); b[o + 1] = (byte)(v >> 16);
            b[o + 2] = (byte)(v >> 8); b[o + 3] = (byte)v;
        }

        static uint Adler32(byte[] data)
        {
            uint a = 1, b = 0;
            foreach (byte x in data) { a = (a + x) % 65521; b = (b + a) % 65521; }
            return (b << 16) | a;
        }

        static readonly uint[] CrcTable = BuildCrc();

        static uint[] BuildCrc()
        {
            var t = new uint[256];
            for (uint i = 0; i < 256; i++)
            {
                uint c = i;
                for (int k = 0; k < 8; k++) c = (c & 1) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
                t[i] = c;
            }
            return t;
        }

        static uint Crc32(byte[] data)
        {
            uint c = 0xFFFFFFFFu;
            foreach (byte b in data) c = CrcTable[(c ^ b) & 0xFF] ^ (c >> 8);
            return c ^ 0xFFFFFFFFu;
        }
    }
}
