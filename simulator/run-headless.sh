#!/usr/bin/env bash
#
# Run the geology core without Unity: simulate, render PNG maps, check invariants.
#
#   ./run-headless.sh                                  defaults
#   ./run-headless.sh --steps 2400 --res 320           a longer, finer run
#   ./run-headless.sh --scenario radial --satellites 8 the exhibit topology
#   ./run-headless.sh bench                            resolution sweep
#   ./run-headless.sh benchplates                      plate-count sweep
#
# Run with --help for the full flag list.
#
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d thirdparty/unity-mathematics/src ]; then
  echo "Dependencies missing. Running scripts/fetch-deps.sh ..." >&2
  ./scripts/fetch-deps.sh
fi

DOTNET="${DOTNET:-$(command -v dotnet || echo /tmp/dotnet/dotnet)}"
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1

# Forward ALL arguments. This previously forwarded only two, which silently dropped
# --res and made the documented `bench` command throw a FormatException.
exec "$DOTNET" run --project headless/Headless.csproj -c Release -- "$@"
