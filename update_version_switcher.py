# .github/scripts/update_version_switcher.py
import json
from pathlib import Path
from packaging.version import Version, InvalidVersion

ROOT = Path(__file__).resolve().parents[2]  # repo root
DOCS_DIR = ROOT / "docs"
OUT_FILE = ROOT / "version_switcher.json"

BASE_URL = "https://martinpdes.github.io/PyMieSim/docs"

# Collect all directories in docs/
versions = []
for p in DOCS_DIR.iterdir():
    if not p.is_dir():
        continue
    name = p.name
    if name == "latest":
        continue
    # Optional: only treat things that look like versions
    try:
        # Strip leading "v" if present: v3.8.6 -> 3.8.6
        ver_str = name[1:] if name.startswith("v") else name
        _ = Version(ver_str)
    except InvalidVersion:
        continue
    versions.append(name)

# Sort versions (highest first) using packaging.version
def sort_key(name: str):
    ver_str = name[1:] if name.startswith("v") else name
    return Version(ver_str)

versions = sorted(versions, key=sort_key, reverse=True)

entries = []

# Always add "latest" first
entries.append({
    "name": "latest",
    "version": "latest",
    "url": f"{BASE_URL}/latest",
})

# Then add all other versions
for v in versions:
    entries.append({
        "name": v,
        "version": v,
        "url": f"{BASE_URL}/{v}",
    })

OUT_FILE.write_text(json.dumps(entries, indent=2) + "\n")
print(f"Wrote {OUT_FILE} with {len(entries)} entries")
