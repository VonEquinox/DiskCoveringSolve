#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUNDLE="$ROOT/proof_bundle"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  if [[ -x "$ROOT/../.venv/bin/python" ]]; then
    PYTHON_BIN="$ROOT/../.venv/bin/python"
  else
    PYTHON_BIN=python3
  fi
fi

PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  echo "Python interpreter not found. Set PYTHON_BIN=/absolute/path/to/python." >&2
  exit 127
fi

cd "$BUNDLE"

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c SHA256SUMS
else
  shasum -a 256 -c SHA256SUMS
fi

"$PYTHON_BIN" cover12_master_verify.py

echo "COVER12 EXACT MASTER REPLAY PASSED"
