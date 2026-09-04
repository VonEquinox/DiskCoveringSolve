#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-quick}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
SOURCE="$ROOT/proof_bundle"
PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  if [[ -x "$ROOT/.venv/bin/python" ]]; then
    PYTHON_BIN="$ROOT/.venv/bin/python"
  elif [[ -x "$ROOT/../.venv/bin/python" ]]; then
    PYTHON_BIN="$ROOT/../.venv/bin/python"
  else
    PYTHON_BIN=python3
  fi
fi

# Resolve both a command name (for example, python3) and a user-supplied path
# to the executable that will still be valid after we change directories.
PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  echo "Python interpreter not found or not executable. Set PYTHON_BIN=/absolute/path/to/python." >&2
  exit 127
fi

if [[ "$MODE" != "quick" && "$MODE" != "full" ]]; then
  echo "Usage: $0 [quick|full]" >&2
  exit 2
fi

"$PYTHON_BIN" - <<'PY'
import importlib
for name in ("numpy", "scipy", "mpmath"):
    importlib.import_module(name)
print("Python dependencies available")
PY

TMP="$(mktemp -d "${TMPDIR:-/tmp}/cover11-replay.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
cp -a "$SOURCE" "$TMP/data"
DATA="$TMP/data"

DATA="$DATA" PYTHON_BIN_FOR_PATCH="$PYTHON_BIN" "$PYTHON_BIN" - <<'PY'
import os
from pathlib import Path
root = Path(os.environ["DATA"])
for path in root.glob("*.py"):
    text = path.read_text()
    path.write_text(text.replace("/mnt/data", str(root)))
script = root / "verify_all.sh"
text = script.read_text()
start = text.index('if [[ "$HERE" != "/mnt/data" ]]; then')
end = text.index('cd /mnt/data', start) + len('cd /mnt/data')
text = text[:start] + 'cd "$HERE"' + text[end:]
text = text.replace("run python ", f"run {os.environ.get('PYTHON_BIN_FOR_PATCH', 'python3')} ")
text = text.replace(
    "printf '\\nALL REQUESTED COVER11 VERIFICATIONS PASSED\\n'",
    "if [[ \"$MODE\" == \"full\" ]]; then\n"
    "  printf '\\nCOVER11 FULL LEAF-BY-LEAF REPLAY PASSED\\n'\n"
    "else\n"
    "  printf '\\nCOVER11 QUICK AUDIT PASSED (STORED LEAF SUMMARIES ONLY)\\n'\n"
    "fi",
)
script.write_text(text)
PY

# Patch the interpreter with a shell-safe absolute path.
"$PYTHON_BIN" - "$DATA/verify_all.sh" "$PYTHON_BIN" <<'PY'
from pathlib import Path
import shlex, sys
path = Path(sys.argv[1])
python_bin = shlex.quote(str(Path(sys.argv[2]).resolve()))
text = path.read_text().replace("run python3 ", f"run {python_bin} ")
path.write_text(text)
PY
chmod +x "$DATA/verify_all.sh"

cd "$DATA"
exec ./verify_all.sh "$MODE"
