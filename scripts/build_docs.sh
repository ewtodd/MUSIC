#!/usr/bin/env bash
# Build the MUSIC API reference.
#
# tooling/Makefile is the internal per-dataset build recipe and is never run
# directly, so the docs get their own entry point rather than a target there.
# Substitutes the @DOC_*@ placeholders in Doxyfile.in and runs doxygen; the
# generated Doxyfile is written into the output directory, not the source tree.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$REPO/docs}"

command -v doxygen >/dev/null 2>&1 || {
  echo "doxygen not found - enter the dev shell first with 'nix develop'" >&2
  exit 1
}
: "${DOXYGEN_AWESOME_CSS:?set DOXYGEN_AWESOME_CSS, or run inside 'nix develop'}"
test -f "$DOXYGEN_AWESOME_CSS/doxygen-awesome.css" || {
  echo "doxygen-awesome-css not found at $DOXYGEN_AWESOME_CSS" >&2
  exit 1
}

# The flake is the version source; fall back to the git description so a
# working-tree build still labels itself usefully.
VERSION="$(git -C "$REPO" describe --tags --always --dirty 2>/dev/null || echo unknown)"
HAVE_DOT="$(command -v dot >/dev/null 2>&1 && echo YES || echo NO)"

mkdir -p "$OUT"
sed -e "s|@DOC_VERSION@|${VERSION}|g" \
    -e "s|@DOC_SOURCE_DIR@|${REPO}|g" \
    -e "s|@DOC_OUTPUT_DIR@|${OUT}|g" \
    -e "s|@DOC_HAVE_DOT@|${HAVE_DOT}|g" \
    -e "s|@DOC_AWESOME_DIR@|${DOXYGEN_AWESOME_CSS}|g" \
    "$REPO/Doxyfile.in" > "$OUT/Doxyfile"

doxygen "$OUT/Doxyfile"
echo "API documentation: $OUT/html/index.html"
echo "Doxygen warnings:  $OUT/doxygen-warnings.log"
