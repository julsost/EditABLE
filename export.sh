#!/usr/bin/env bash
set -euo pipefail  # stop on error, undefined vars, and fail on pipe errors

SRC_DIR="editable_app"
EXPORT_DIR="app_for_export"
DOCS_DIR="docs"

echo "=== Cleaning old export ==="
rm -rf "$EXPORT_DIR" "$DOCS_DIR"/*
mkdir -p "$EXPORT_DIR"

echo "=== Copying Python modules ==="
cp "$SRC_DIR"/{app.py,utils.py,calculate_scores.py,make_json.py,pegfinder.py,__init__.py} "$EXPORT_DIR"/

echo "=== Copying package directories ==="
cp -r "$SRC_DIR"/jsonscript "$EXPORT_DIR"/

echo "=== Copying image assets ==="
cp "$SRC_DIR"/{EditABLE-logos_transparent.png,integrase_diagram.png,prime_editing_diagram.png,SOM_Web_vert_LG.png} "$EXPORT_DIR"/

echo "=== Exporting with Shinylive ==="
shinylive export "$EXPORT_DIR" "$DOCS_DIR"

echo "=== Patching index.html for dynamic relPath ==="
perl -0777 -i -pe '
s#<script type="module">\s*import\s+runExportedApp\s+from\s+"\.\/shinylive\/shinylive\.js";\s*runExportedApp\(\{.*?}\);\s*</script>#<script type="module">
  import runExportedApp from "./shinylive/shinylive.js";
  const relPath=((p)=>{const i=p.lastIndexOf("/");return i<=0?"":p.slice(1,i);})(window.location.pathname);
  runExportedApp({ id: "root", appEngine: "python", relPath: relPath || "." });
</script>#s' "$DOCS_DIR/index.html"




echo "=== Export complete ==="
echo "Test locally with:"
echo "  python3 -m http.server --directory \"$DOCS_DIR\" --bind localhost 8008"