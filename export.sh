#!/usr/bin/env bash
set -e  # stop on error

# Define folders
SRC_DIR="editable_app"
EXPORT_DIR="app_for_export"
DOCS_DIR="docs"

echo "=== Cleaning old export ==="
rm -rf "$EXPORT_DIR" "$DOCS_DIR"/*
mkdir -p "$EXPORT_DIR"

echo "=== Copying Python modules ==="
cp "$SRC_DIR"/app.py "$EXPORT_DIR"/
cp "$SRC_DIR"/utils.py "$EXPORT_DIR"/
cp "$SRC_DIR"/calculate_scores.py "$EXPORT_DIR"/
cp "$SRC_DIR"/make_json.py "$EXPORT_DIR"/
cp "$SRC_DIR"/pegfinder.py "$EXPORT_DIR"/
cp "$SRC_DIR"/__init__.py "$EXPORT_DIR"/

echo "=== Copying package directories ==="
cp -r "$SRC_DIR"/jsonscript "$EXPORT_DIR"/

echo "=== Copying image assets ==="
cp "$SRC_DIR"/EditABLE-logos_transparent.png "$EXPORT_DIR"/
cp "$SRC_DIR"/integrase_diagram.png "$EXPORT_DIR"/
cp "$SRC_DIR"/prime_editing_diagram.png "$EXPORT_DIR"/
cp "$SRC_DIR"/SOM_Web_vert_LG.png "$EXPORT_DIR"/

echo "=== Exporting to docs/ ==="
shinylive export "$EXPORT_DIR" "$DOCS_DIR"

echo "=== Export complete ==="
echo "Now you can test locally with:"
echo "  python3 -m http.server --directory $DOCS_DIR --bind localhost 8008"