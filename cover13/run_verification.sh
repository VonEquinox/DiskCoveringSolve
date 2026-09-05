#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUNDLE="$ROOT/proof_bundle"

PYTHON_BIN="${PYTHON_BIN:-python3}"
PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  echo "Python 3.10+ not found. Set PYTHON_BIN=/absolute/path/to/python." >&2
  exit 127
fi

# Homebrew installs Boost as keg-only on some macOS versions.  Expose the
# headers automatically without changing the certificate verifier itself.
if [[ "$(uname -s)" == "Darwin" ]] && command -v brew >/dev/null 2>&1; then
  BOOST_PREFIX="$(brew --prefix boost 2>/dev/null || true)"
  if [[ -n "$BOOST_PREFIX" && -f "$BOOST_PREFIX/include/boost/multiprecision/cpp_int.hpp" ]]; then
    export CPLUS_INCLUDE_PATH="$BOOST_PREFIX/include${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
  fi
fi

"$PYTHON_BIN" "$ROOT/verify_manifest.py"
cd "$BUNDLE"
"$PYTHON_BIN" -B verify_all.py

echo "COVER13 EXACT MASTER REPLAY PASSED"
