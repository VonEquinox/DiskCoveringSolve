#!/usr/bin/env python3
"""Restore immutable split inputs and replay every original Cover20 stage."""
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
LARGE = "proof_bundle/enumeration/noncandidate20_trees.jsonl.gz"
STAGES = {"verify_root20.py", "verify_upper20.py", "anchor_isolation20.py",
          "verify_interfaces20.py", "compile_peel_dfs20", "compile_peel_bfs20",
          "compile_compare_maps20", "verify_peeling20.py", "verify_recursion_tests20.py",
          "verify_forest20.py", "test_fail_closed20.py"}
STAGES.update(f"{mode}_B{b}" for mode in ("dfs", "bfs") for b in range(12, 20))


def safe_path(rel, folder):
    path = Path(rel)
    if path.is_absolute() or ".." in path.parts or len(path.parts) < 2 or path.parts[0] != folder:
        raise ValueError(f"Unsafe input path: {rel}")
    return ROOT / path


def consume(path, total, destination=None):
    own = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            own.update(block)
            total.update(block)
            size += len(block)
            if destination is not None:
                destination.write(block)
    return size, own.hexdigest()


def verify_manifest():
    inputs = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
    storage = json.loads((ROOT / "STORAGE_MANIFEST.json").read_text())
    if len(inputs) != 51 or set(storage) != {LARGE} or LARGE not in inputs:
        raise RuntimeError("Expected 51 logical files and one split certificate")
    physical = set()
    for rel, expected in inputs.items():
        path = safe_path(rel, "proof_bundle")
        total = hashlib.sha256()
        if rel in storage:
            entry = storage[rel]
            if set(entry) != {"bytes", "parts"} or len(entry["parts"]) != 3:
                raise RuntimeError("Invalid split manifest")
            size = 0
            for part in entry["parts"]:
                if set(part) != {"path", "bytes", "sha256"} or part["path"] in physical:
                    raise RuntimeError("Invalid or repeated part")
                part_path = safe_path(part["path"], "proof_parts")
                count, digest = consume(part_path, total)
                if count != part["bytes"] or digest != part["sha256"]:
                    raise RuntimeError(f"Part hash or length mismatch: {part_path}")
                size += count
                physical.add(part["path"])
            if size != entry["bytes"]:
                raise RuntimeError("Reconstructed length mismatch")
        else:
            consume(path, total)
            physical.add(rel)
        if total.hexdigest() != expected:
            raise RuntimeError(f"Original input SHA-256 mismatch: {rel}")
    actual = {str(p.relative_to(ROOT)) for folder in ("proof_bundle", "proof_parts")
              for p in (ROOT / folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    if actual != physical:
        raise RuntimeError("Unexpected physical input file set")
    if (ROOT / "docs/r20_proof_zh.md").read_bytes() != (ROOT / "proof_bundle/PROOF_zh.md").read_bytes():
        raise RuntimeError("Published proof differs from the ZIP proof")
    print("COVER20 ORIGINAL INPUT MANIFEST PASSED (51 logical files; 3 verified parts)", flush=True)
    return inputs, storage


def materialize(work, inputs, storage):
    shutil.copytree(ROOT / "proof_bundle", work, ignore=shutil.ignore_patterns("__pycache__"))
    for rel, entry in storage.items():
        target = work / Path(rel).relative_to("proof_bundle")
        digest = hashlib.sha256()
        size = 0
        with target.open("xb") as output:
            for part in entry["parts"]:
                count, part_hash = consume(safe_path(part["path"], "proof_parts"), digest, output)
                if count != part["bytes"] or part_hash != part["sha256"]:
                    raise RuntimeError("Part changed during reconstruction")
                size += count
        if digest.hexdigest() != inputs[rel] or size != entry["bytes"]:
            raise RuntimeError("Reconstruction differs from original ZIP input")
    for rel, expected in inputs.items():
        digest = hashlib.sha256()
        consume(work / Path(rel).relative_to("proof_bundle"), digest)
        if digest.hexdigest() != expected:
            raise RuntimeError(f"Materialized input mismatch: {rel}")


def report_target(value):
    if value is None:
        return None
    target = Path(value).resolve()
    if target.exists() or any(target.is_relative_to((ROOT / folder).resolve()) for folder in ("proof_bundle", "proof_parts")):
        raise ValueError("Reports must use a new directory outside immutable inputs")
    return target


def validate_result(result):
    if result.get("verified") is not True or result.get("unresolved_leaves") != 0:
        raise RuntimeError("Missing complete successful replay")
    stages = [s.get("step") for s in result.get("steps", [])]
    if len(stages) != len(STAGES) or set(stages) != STAGES:
        raise RuntimeError("Every original stage must execute exactly once")
    expected = {"kkt_dimension": 198, "surviving_orbits": 562169,
                "rooted_survivors": 14747883, "noncandidate_cases": 562168,
                "candidate_cases": 1, "unresolved_combinatorial_branches": 0,
                "recurrence_regression_cases": 33, "anchor_isolation_radius": "11/100"}
    if any(result.get(key) != value for key, value in expected.items()):
        raise RuntimeError("Unexpected or incomplete Cover20 result")
    if result.get("negative_tests", {}).get("negative_test_count") != 43:
        raise RuntimeError("Missing negative tests")
    candidate = result.get("candidate_tree", {})
    noncandidate = result.get("noncandidate_forest", {})
    if (candidate.get("nodes"), candidate.get("counts", {}).get("LOCAL"),
            noncandidate.get("nodes")) != (6171, 102, 657498):
        raise RuntimeError("Incorrect candidate or noncandidate tree totals")


def replay(workers, output=None):
    inputs, storage = verify_manifest()
    output = report_target(output)
    with tempfile.TemporaryDirectory(prefix="cover20_replay_") as td:
        work = Path(td) / "bundle"
        materialize(work, inputs, storage)
        # The original archive contains no success reports. Always start fresh.
        if (work / "reports").exists():
            shutil.rmtree(work / "reports")
        for old in (work / "core").glob("*_verified.json"):
            old.unlink()
        env = dict(os.environ)
        env.pop("PYTHONOPTIMIZE", None)
        env.pop("PYTHONPATH", None)
        env["COVER20_WORKERS"] = str(workers)
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        subprocess.run([sys.executable, "-S", "-B", "verify_all.py"], cwd=work, env=env, check=True)
        result = json.loads((work / "reports/MASTER_VERIFIED.json").read_text())
        validate_result(result)
        from verify_interfaces_extra import verify
        extra = verify(work)
        verify_manifest()
        if output is not None:
            shutil.copytree(work / "reports", output)
            (output / "independent_interfaces.json").write_text(json.dumps(extra, indent=2) + "\n")
    print("COVER20 EXACT MASTER REPLAY AND INPUT RECONSTRUCTION PASSED", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 17))
    parser.add_argument("--report-dir", type=Path)
    args = parser.parse_args()
    if args.integrity_only:
        verify_manifest()
        return
    if sys.version_info < (3, 10):
        raise RuntimeError("Python 3.10+ required")
    replay(args.workers, args.report_dir)


if __name__ == "__main__":
    main()
