#!/bin/bash
# Clean all build artifacts, then (unless --ignore-highs-verification) verify that the vendored
# HiGHS 1.15.0 source matches the official upstream release.
#
# Steps:
#   1. Remove every build artifact — the MAF build (cpp/build) AND the HiGHS build/install trees
#      under the vendored source (so HiGHS-1.15.0/ is left as pristine source, comparable to a
#      fresh extraction of the official release).
#   2. Ensure ~/HiGHS-1.15.0.zip is present; if not, download it from the official GitHub release
#      tag (wget -O ~/HiGHS-1.15.0.zip <v1.15.0.zip>).
#   3. Extract it to /tmp (-> /tmp/HiGHS-1.15.0).
#   4. Diff it file-by-file against src/dependencies/HiGHS-1.15.0.
#
# Options:
#   --ignore-highs-verification   do step 1 only (remove build artifacts); skip verification.
#
# Exit 0 iff the vendored source is byte-identical to the upstream release (or --ignore- given).
set -euo pipefail

VERIFY=1
for arg in "$@"; do
  case "$arg" in
    --ignore-highs-verification) VERIFY=0 ;;
    *) echo "error: unknown option '$arg' (expected --ignore-highs-verification)" >&2; exit 1 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENDORED="$ROOT/src/dependencies/HiGHS-1.15.0"
ZIP="$HOME/HiGHS-1.15.0.zip"
URL="https://github.com/ERGO-Code/HiGHS/archive/refs/tags/v1.15.0.zip"
TMP="/tmp/HiGHS-1.15.0"

# --- 1. remove all build artifacts (MAF builds + the prebuilt-HiGHS cache) ----------------------
echo "[1/4] removing build artifacts ..."
# All MAF build trees (build, build_*), the persistent prebuilt-HiGHS cache (.highs / .highs-build,
# compiled once by docker_setup.sh), and any legacy HiGHS build/install trees under the vendored source — so
# HiGHS-1.15.0/ is left as pristine source.  Only these directories (never source files such as
# build_webdemo.sh) are removed.  NB: this forces the NEXT build to recompile HiGHS from scratch.
rm -rf "$ROOT"/build "$ROOT"/build_* "$ROOT/.highs" "$ROOT/.highs-build"
for d in build build-static install install-static; do
  [ -d "$VENDORED/$d" ] && rm -rf "$VENDORED/$d"
done

if [ "$VERIFY" -eq 0 ]; then
  echo "cleanup done (--ignore-highs-verification: skipped HiGHS source verification)."
  exit 0
fi

# --- 2. ensure the reference archive ~/HiGHS-1.15.0.zip ----------------------------------------
if [ -f "$ZIP" ]; then
  echo "[2/4] using cached $ZIP"
else
  command -v wget >/dev/null 2>&1 || { echo "error: wget is required to fetch $URL" >&2; exit 1; }
  echo "[2/4] $ZIP not found; downloading ..."
  wget -O "$ZIP" "$URL"
fi

# --- 3. extract the reference to /tmp ----------------------------------------------------------
echo "[3/4] extracting reference to $TMP ..."
rm -rf "$TMP"
unzip -q "$ZIP" -d /tmp
[ -d "$TMP" ] || { echo "error: expected $TMP after extraction" >&2; exit 1; }

# --- 4. compare file-by-file against the vendored source ---------------------------------------
echo "[4/4] comparing vendored $VENDORED against upstream v1.15.0 ..."
[ -d "$VENDORED" ] || { echo "error: vendored source missing: $VENDORED" >&2; exit 1; }
if diff -rq "$TMP" "$VENDORED"; then
  echo "OK: vendored HiGHS 1.15.0 source is IDENTICAL to upstream v1.15.0"
  rc=0
else
  echo "MISMATCH: the differences listed above exist between upstream and the vendored source" >&2
  rc=1
fi

rm -rf "$TMP"
exit $rc
