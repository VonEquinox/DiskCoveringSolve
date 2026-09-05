#!/usr/bin/env python3
"""Verify fixed archives and replay Cover15 in disposable working directories."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
BUNDLES = {"primary": "proof_bundle", "alternative": "alternative_bundle"}


def verify_manifest():
    manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
    if len(manifest) != 146:
        raise RuntimeError("Expected 146 archived files")
    for rel, expected in manifest.items():
        name = Path(rel)
        if name.is_absolute() or ".." in name.parts:
            raise RuntimeError(f"Invalid manifest path: {rel}")
        digest = hashlib.sha256()
        with (ROOT / name).open("rb") as source:
            for block in iter(lambda: source.read(1 << 20), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise RuntimeError(f"SHA-256 mismatch: {rel}")
    if (ROOT / "docs/r15_proof_zh.md").read_bytes() != (ROOT / "alternative_bundle/PROOF_zh.md").read_bytes():
        raise RuntimeError("Standalone manuscript differs from alternative package")
    print("COVER15 ARCHIVE MANIFEST PASSED (146 files)", flush=True)


def environment():
    env = os.environ.copy()
    if sys.platform == "darwin" and shutil.which("brew"):
        result = subprocess.run(["brew", "--prefix", "boost"], capture_output=True, text=True)
        if result.returncode == 0:
            include = Path(result.stdout.strip()) / "include"
            if (include / "boost/multiprecision/cpp_int.hpp").is_file():
                env["CPLUS_INCLUDE_PATH"] = str(include) + (
                    os.pathsep + env["CPLUS_INCLUDE_PATH"] if env.get("CPLUS_INCLUDE_PATH") else ""
                )
    return env


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=[*BUNDLES, "both"], default="both")
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=2, choices=range(1, 17), help="Primary C++ concurrency; alternative uses two workers")
    args = parser.parse_args()
    verify_manifest()
    if args.integrity_only:
        return
    env = environment()
    variants = list(BUNDLES) if args.variant == "both" else [args.variant]
    for variant in variants:
        with tempfile.TemporaryDirectory(prefix=f"cover15_{variant}_") as tmp:
            work = Path(tmp) / "bundle"
            shutil.copytree(ROOT / BUNDLES[variant], work, ignore=shutil.ignore_patterns("__pycache__", "replay_workspace", "_runtime"))
            command = [sys.executable, "-S", "-B", "verify_all.py"]
            if variant == "primary":
                command += ["--regenerate", "--jobs", str(args.jobs)]
            subprocess.run(command, cwd=work, env=env, check=True)
        print(f"COVER15 {variant.upper()} EXACT MASTER REPLAY PASSED", flush=True)
    subprocess.run([sys.executable, "-S", "-B", str(ROOT / "verify_bundle_linkage.py")], check=True)
    verify_manifest()
    print("COVER15 REQUESTED REPLAYS AND ROOT LINKAGE PASSED", flush=True)


if __name__ == "__main__":
    main()
