#!/usr/bin/env python3
"""Verify the fixed SHA-256 manifest for the cover14 bundle."""
from __future__ import annotations

import hashlib
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parent / "proof_bundle"
manifest = ROOT / "SHA256SUMS"
count = 0

for raw in manifest.read_text(encoding="ascii").splitlines():
    if not raw:
        continue
    expected, rel = raw.split("  ", 1)
    logical = PurePosixPath(rel)
    if logical.is_absolute() or ".." in logical.parts:
        raise SystemExit(f"INVALID PATH: {rel}")
    path = ROOT / logical
    if not path.is_file():
        raise SystemExit(f"MISSING: {rel}")
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != expected:
        raise SystemExit(
            f"HASH MISMATCH: {rel}\nexpected {expected}\nactual   {got}"
        )
    count += 1

print(f"COVER14 FIXED MANIFEST PASSED ({count} files)")
