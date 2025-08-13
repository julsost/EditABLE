import json
import os

# Directory containing your files
source_dir = "."  # change if your files are in a subfolder

# Extensions to include
include_exts = {".py", ".txt", ".md"}

output = []

for filename in os.listdir(source_dir):
    filepath = os.path.join(source_dir, filename)

    if os.path.isfile(filepath) and os.path.splitext(filename)[1] in include_exts:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()
        output.append({
            "name": filename,
            "content": content
        })

# Write to JSON
with open("files.json", "w", encoding="utf-8") as out_f:
    json.dump(output, out_f, ensure_ascii=False, indent=2)

print(f"Saved {len(output)} files to files.json")
