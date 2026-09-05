"""Synthetic runner tests, separate from the mathematical certificate replay."""
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_verification as runner


class RunnerTests(unittest.TestCase):
    def setUp(self):
        output = patch("sys.stdout", new_callable=io.StringIO)
        output.start()
        self.addCleanup(output.stop)
        temp = tempfile.TemporaryDirectory(prefix="cover19_runner_test_")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        (self.root / "docs").mkdir()
        for variant, folder in runner.BUNDLES.items():
            bundle = self.root / folder
            (bundle / "core").mkdir(parents=True)
            names = ["PROOF_zh.md"]
            if variant == "alternative":
                names += ["reports/MASTER_VERIFIED.json", "reports/legacy.log", "core/fixture_verified.json"]
            names += [f"input{i}" for i in range(runner.COUNTS[variant]-len(names))]
            for name in names:
                path = bundle / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("fixture")
            proof = "r19_proof_zh.md" if variant == "primary" else "r19_alternative_proof_zh.md"
            (self.root / "docs" / proof).write_text("fixture")
        manifest = {str(p.relative_to(self.root)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for folder in runner.BUNDLES.values() for p in (self.root / folder).rglob("*") if p.is_file()}
        (self.root / "INPUT_MANIFEST.json").write_text(json.dumps(manifest))
        root_patch = patch.object(runner, "ROOT", self.root)
        root_patch.start()
        self.addCleanup(root_patch.stop)
        self.variant = "primary"
        self.corrupt = None

    def fake_run(self, command, *, cwd, check, env):
        self.assertEqual(command[-3:], ["-S", "-B", "verify_all.py"])
        self.assertTrue(check)
        self.assertNotIn("PYTHONOPTIMIZE", env)
        self.assertNotIn("PYTHONPATH", env)
        self.work = Path(cwd)
        self.assertEqual(list((self.work / "reports").iterdir()), [])
        self.assertFalse(list((self.work / "core").glob("*_verified.json")))
        key = "step" if self.variant == "primary" else "stage"
        orbit_key = "surviving_topology_orbits" if self.variant == "primary" else "surviving_orbits"
        result = {"verified": True, "unresolved_leaves": 0,
                  "exact_squared_radius": "1/13", "exact_radius_squared": "1/13",
                  orbit_key: 24127, "noncandidate_forest": {"cases": 24126, "noncandidate_cases": 24126},
                  "steps": [{key: s} for s in runner.STEPS[self.variant]],
                  "rejection_tests": {"rejection_tests": 35}, "negative_tests": {"rejection_tests": 36},
                  "candidate_global_lower": {"positive_rods": 42}, "candidate_cases": 1,
                  "positive_stress_graph": {"geometric_dimension": 70},
                  "unpruned_recurrence_tests": 26 if self.variant == "primary" else {"cases": 29}}
        if self.corrupt:
            self.corrupt(result)
        (self.work / "reports/MASTER_VERIFIED.json").write_text(json.dumps(result))

    def test_manifest_and_proofs(self):
        runner.verify_manifest()

    def test_corrupted_input(self):
        (self.root / "proof_bundle/input0").write_text("modified")
        with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
            runner.verify_manifest()

    def test_missing_input(self):
        (self.root / "alternative_bundle/input0").unlink()
        with self.assertRaises(FileNotFoundError):
            runner.verify_manifest()

    def test_extra_input(self):
        (self.root / "proof_bundle/extra.py").write_text("modified")
        with self.assertRaisesRegex(RuntimeError, "Unexpected archive file set"):
            runner.verify_manifest()

    def test_unsafe_path(self):
        path = self.root / "INPUT_MANIFEST.json"
        manifest = json.loads(path.read_text())
        manifest["../outside"] = manifest.pop("proof_bundle/input0")
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(RuntimeError, "Invalid manifest path"):
            runner.verify_manifest()

    def test_proof_mismatch(self):
        (self.root / "docs/r19_alternative_proof_zh.md").write_text("modified")
        with self.assertRaisesRegex(RuntimeError, "Standalone proof"):
            runner.verify_manifest()

    def test_export_and_cleanup_both_variants(self):
        for variant in runner.BUNDLES:
            self.variant = variant
            output = self.root / f"fresh_{variant}"
            with patch.object(runner.subprocess, "run", side_effect=self.fake_run):
                runner.replay(variant, output)
            self.assertFalse(self.work.exists())
            self.assertTrue((output / "MASTER_VERIFIED.json").exists())
            self.assertTrue((output / "core").is_dir())
            self.assertFalse((output / "legacy.log").exists())
            runner.verify_manifest()

    def test_stale_success(self):
        with patch.object(runner.subprocess, "run", return_value=None):
            with self.assertRaises(FileNotFoundError):
                runner.replay("alternative")

    def test_missing_stage(self):
        self.corrupt = lambda result: result["steps"].pop()
        with patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            with self.assertRaisesRegex(RuntimeError, "Every accepting stage"):
                runner.replay("primary")

    def test_unresolved_result(self):
        self.corrupt = lambda result: result.update(unresolved_leaves=1)
        with patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            with self.assertRaisesRegex(RuntimeError, "Missing successful"):
                runner.replay("primary")

    def test_missing_tests(self):
        self.corrupt = lambda result: result.update(unpruned_recurrence_tests=0)
        with patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            with self.assertRaisesRegex(RuntimeError, "Missing regression"):
                runner.replay("primary")

    def test_protected_output(self):
        for target in [self.root / "docs", self.root / "proof_bundle/new", self.root / "alternative_bundle/new"]:
            with self.assertRaises(ValueError):
                runner.replay("primary", target)


if __name__ == "__main__":
    unittest.main()
