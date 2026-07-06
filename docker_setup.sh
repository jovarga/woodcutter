#!/bin/bash
# docker_setup.sh — all steps to build the MAF solver from source in a fresh Debian 13.5 container.

# This setup was tested inside a clean Debian 13.5 Docker environment created with:
#     docker run -it --rm debian:13.5
#
# CMake (CMakeLists.txt) is the single source of truth for configuration and track selection; this
# script installs the toolchain and drives CMake directly (no build.sh wrapper).  It builds ALL THREE
# competition tracks, leaving each static binary in its track build directory:
#     build/<track>/main_<track>     (exact | heuristic | lowerbound)
#
# HiGHS 1.15.0 is vendored in src/dependencies/, so NO network is needed for the build itself (only the
# apt-get toolchain install needs the Debian mirrors).  It is compiled ONCE into ./.highs and then merely
# LINKED by each track build, so the three builds do not recompile it.
#
# Usage:
#   bash docker_setup.sh                       # build all three tracks
#   bash docker_setup.sh heuristic lowerbound  # build only the listed track(s)
#
# To build a single track by hand (CMake directly, no this script):
#   cmake -S . -B build/<track> -DCMAKE_BUILD_TYPE=Release -DMAF_TRACK=<track>
#   cmake --build build/<track> -j --target maf          # -> build/<track>/main_<track>
#   (<track> is one of exact | heuristic | lowerbound; add -DHIGHS_INSTALL_DIR="$PWD/.highs" to link a
#    prebuilt HiGHS instead of rebuilding it in-tree.)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
HIGHS_DIR="$HERE/.highs"
HIGHS_SRC="$HERE/src/dependencies/HiGHS-1.15.0"

# Tracks to build (default: all three; or the ones named on the command line).
TRACKS=("$@"); [ ${#TRACKS[@]} -eq 0 ] && TRACKS=(exact heuristic lowerbound)
for t in "${TRACKS[@]}"; do
  case "$t" in exact|heuristic|lowerbound) ;;
    *) echo "error: unknown track '$t' (expected exact | heuristic | lowerbound)" >&2; exit 1 ;;
  esac
done

# 1. Toolchain: g++/gcc/make/libc-dev (build-essential) + cmake.  HiGHS is vendored and ZLIB is OFF, so
#    no zlib/unzip/wget are required.
if command -v apt-get >/dev/null 2>&1; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update
  apt-get install -y --no-install-recommends build-essential cmake
else
  echo "warning: apt-get not found; assuming build-essential + cmake are already installed" >&2
fi

command -v cmake >/dev/null 2>&1 || { echo "error: cmake is required" >&2; exit 1; }

# 2. Build HiGHS 1.15.0 ONCE into ./.highs (skip if already installed), so the per-track configures LINK
#    it instead of each rebuilding it.  These flags MUST match the in-tree fallback in CMakeLists.txt:
#    static, no zlib/BLAS, single-threaded.  Set FORCE_HIGHS=1 to rebuild it from scratch.
if [ "${FORCE_HIGHS:-0}" != "0" ]; then rm -rf "$HIGHS_DIR" "$HERE/.highs-build"; fi
if ! find "$HIGHS_DIR" -name 'highs-config.cmake' 2>/dev/null | grep -q .; then
  echo "==> building HiGHS 1.15.0 -> $HIGHS_DIR (one-time; reused by every track build)"
  cmake -S "$HIGHS_SRC" -B "$HERE/.highs-build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$HIGHS_DIR" \
    -DBUILD_SHARED_LIBS=OFF -DZLIB=OFF -DFAST_BUILD=ON -DHIPO=OFF \
    -DBUILD_OPENBLAS=OFF -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
    -DHIGHS_NO_DEFAULT_THREADS=ON
  cmake --build "$HERE/.highs-build" -j"$JOBS" --target install
else
  echo "==> reusing prebuilt HiGHS at $HIGHS_DIR (delete it or set FORCE_HIGHS=1 to rebuild)"
fi

# 3. Configure + build each track (linking the prebuilt HiGHS).
for TRACK in "${TRACKS[@]}"; do
  echo "==> building $TRACK track"
  cmake -S "$HERE" -B "$HERE/build/$TRACK" \
    -DCMAKE_BUILD_TYPE=Release -DMAF_TRACK="$TRACK" -DHIGHS_INSTALL_DIR="$HIGHS_DIR"
  cmake --build "$HERE/build/$TRACK" -j"$JOBS" --target maf
  echo "built: $HERE/build/$TRACK/main_${TRACK}"
done

echo "done. tracks: ${TRACKS[*]}   binaries: build/<track>/main_<track>"
