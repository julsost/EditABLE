#!/usr/bin/env bash
# export.sh — build & export EditABLE (Shinylive) and patch index.html for GitHub Pages
# Usage:
#   ./export.sh [--src editable_app] [--export app_for_export] [--docs docs] [--base /RepoName]
#
# Notes:
# - Automatically computes relPath at runtime so the same build works on localhost and GitHub Pages.
# - If you know your repo path (e.g. /EditABLE), you can force it with --base; otherwise it auto-detects.

set -euo pipefail

#######################################
# Defaults (override via flags)
#######################################
SRC_DIR="editable_app"
EXPORT_DIR="app_for_export"
DOCS_DIR="docs"
FORCED_BASE=""

#######################################
# Parse flags
#######################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    --src)    SRC_DIR="${2:?}"; shift 2 ;;
    --export) EXPORT_DIR="${2:?}"; shift 2 ;;
    --docs)   DOCS_DIR="${2:?}"; shift 2 ;;
    --base)   FORCED_BASE="${2:?}"; shift 2 ;;
    -h|--help)
      sed -n '1,25p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

#######################################
# Helpers
#######################################
need() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command: $1" >&2
    exit 1
  }
}

copy_if_exists() {
  # Usage: copy_if_exists <src> <dest_dir>
  local src="$1"
  local dest="$2"
  if [[ -e "$src" ]]; then
    cp -R "$src" "$dest"/
  else
    echo "(!) Skipping missing: $src"
  fi
}

preserve_docs_specials() {
  # Preserve CNAME and .nojekyll if they exist
  local dir="$1"
  local tmp="$(mktemp -d)"
  for f in CNAME .nojekyll; do
    if [[ -f "$dir/$f" ]]; then
      mv "$dir/$f" "$tmp"/
    fi
  done
  rm -rf "$dir"/*
  mkdir -p "$dir"
  for f in "$tmp"/*; do
    [[ -e "$f" ]] && mv "$f" "$dir"/
  done
  rmdir "$tmp"
}

#######################################
# Checks
#######################################
need shinylive
need perl

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Source directory not found: $SRC_DIR" >&2
  exit 1
fi

#######################################
# Clean/prepare
#######################################
echo "=== Cleaning old export ==="
rm -rf "$EXPORT_DIR"
mkdir -p "$EXPORT_DIR"
mkdir -p "$DOCS_DIR"
preserve_docs_specials "$DOCS_DIR"

#######################################
# Copy app sources
#######################################
echo "=== Copying Python modules ==="
for f in \
  "$SRC_DIR/app.py" \
  "$SRC_DIR/utils.py" \
  "$SRC_DIR/calculate_scores.py" \
  "$SRC_DIR/make_json.py" \
  "$SRC_DIR/pegfinder.py" \
  "$SRC_DIR/__init__.py"
do
  copy_if_exists "$f" "$EXPORT_DIR"
done

echo "=== Copying package directories ==="
copy_if_exists "$SRC_DIR/jsonscript" "$EXPORT_DIR"

echo "=== Copying image assets ==="
for f in \
  "$SRC_DIR/EditABLE-logos_transparent.png" \
  "$SRC_DIR/integrase_diagram.png" \
  "$SRC_DIR/prime_editing_diagram.png" \
  "$SRC_DIR/SOM_Web_vert_LG.png"
do
  copy_if_exists "$f" "$EXPORT_DIR"
done

#######################################
# Export with Shinylive
#######################################
echo "=== Exporting with Shinylive ==="
shinylive export "$EXPORT_DIR" "$DOCS_DIR"

INDEX="$DOCS_DIR/index.html"
if [[ ! -f "$INDEX" ]]; then
  echo "Export failed: $INDEX not found" >&2
  exit 1
fi

#######################################
# Patch index.html
# - Insert a small script that computes window.__SL_SCOPE__
# - Rewrite runExportedApp({...}) to set relPath: window.__SL_SCOPE__ || ""
#######################################
echo "=== Patching index.html to auto-detect relPath ==="

# 1) Insert scope helper once, right before the FIRST <script type="module"> tag.
#    We match any <script ... type="module" ...> with attrs in any order.
perl -0777 -i -pe '
  my $forced = $ENV{FORCED_BASE} // "";
  my $marker = "BEGIN_REL_PATH_AUTODETECT";
  if (index($_, $marker) == -1) {
    s~(<script\b[^>]*\btype=(["'\''])module\2[^>]*>)~<!-- BEGIN_REL_PATH_AUTODETECT -->
<script>
  // Computes scope for both localhost and GitHub Pages project sites.
  // If FORCED_BASE is provided (e.g., "/EditABLE"), it will be used.
  (function () {
    var forced = '."'$FORCED_BASE'".q{};
    if (forced && typeof forced === "string") {
      window.__SL_SCOPE__ = forced;
      return;
    }
    // directory containing index.html
    var p = location.pathname.replace(/\/index\.html?$/, "");
    // strip trailing slash (but keep leading)
    p = p.replace(/\/$/, "");
    // Special case: root path becomes ""
    window.__SL_SCOPE__ = (p === "" || p === "/") ? "" : p;
  })();
</script>
<!-- END_REL_PATH_AUTODETECT -->
$1~is;
  }
' "$INDEX"

# 2) Replace relPath in runExportedApp options object (idempotent).
perl -0777 -i -pe '
  s/(runExportedApp\s*\(\s*\{\s*[^}]*?)relPath\s*:\s*["'\''][^"'\'']*["'\'']/$1relPath: window.__SL_SCOPE__ || ""/s;
' "$INDEX"

#######################################
# Done
#######################################
echo "=== Done ==="
echo "Serve locally with:"
echo "  python3 -m http.server --directory \"$DOCS_DIR\" --bind localhost 8008"