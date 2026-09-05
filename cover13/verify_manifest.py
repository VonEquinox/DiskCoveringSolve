#!/usr/bin/env python3
"""Verify the immutable SHA-256 input manifest for the cover13 bundle."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "proof_bundle"
manifest = json.loads((ROOT / "MANIFEST_SHA256.json").read_text(encoding="utf-8"))
for rel, expected in manifest.items():
    path = ROOT / rel
    if not path.is_file():
        raise SystemExit(f"MISSING: {rel}")
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != expected:
        raise SystemExit(f"HASH MISMATCH: {rel}\nexpected {expected}\nactual   {got}")
print(f"COVER13 INPUT MANIFEST PASSED ({len(manifest)} files)")
