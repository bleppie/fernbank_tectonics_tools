#!/usr/bin/env bash
#
# Fetch build dependencies that are deliberately NOT committed to this repo.
#
# Unity.Mathematics is licensed under the Unity Companion License, not MIT. That licence
# covers using it in Unity-dependent projects, which this is — but redistributing its source
# inside a public repository is a separate act, so we fetch it instead of vendoring it.
#
# Inside Unity you don't need this at all: add the `com.unity.mathematics` package and the
# `Tectonic.Core` assembly definition picks it up. This script exists only so the headless
# build (which has no Unity install) can compile the same source.
#
set -euo pipefail
cd "$(dirname "$0")/.."

# Pinned so headless builds are reproducible. Bump deliberately, not incidentally.
UNITY_MATHEMATICS_REF="${UNITY_MATHEMATICS_REF:-1.3.2}"
DEST="thirdparty/unity-mathematics"

if [ -d "$DEST/src" ]; then
  echo "Unity.Mathematics already present at $DEST — delete it to re-fetch."
  exit 0
fi

echo "Fetching Unity.Mathematics @ $UNITY_MATHEMATICS_REF ..."
mkdir -p thirdparty
rm -rf "$DEST"
git clone --depth 1 --branch "$UNITY_MATHEMATICS_REF" \
    https://github.com/Unity-Technologies/Unity.Mathematics.git "$DEST" 2>/dev/null \
  || git clone --depth 1 https://github.com/Unity-Technologies/Unity.Mathematics.git "$DEST"

# Unity .meta files are noise outside a Unity project.
find "$DEST" -name '*.meta' -delete

echo "Done. Headless build can now run: ./run-headless.sh"
