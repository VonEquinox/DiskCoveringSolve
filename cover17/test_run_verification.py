"""Synthetic runner regressions, separate from the mathematical full replay."""
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
        temp = tempfile.TemporaryDirectory(prefix="cover17_runner_test_")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.bundle = self.root / "proof_bundle"
        for folder in ["proof_bundle/core", "proof_bundle/reports", "docs"]:
            (self.root / folder).mkdir(parents=True)
        for rel in ["PROOF_zh.md", "reports/MASTER_VERIFIED.json", "reports/legacy.log", "core/fixture_verified.json"]:
            (self.bundle / rel).write_text("fixture")
        for i in range(42):
            (self.bundle / f"input{i}").write_text(str(i))
        (self.root / "docs/r17_proof_zh.md").write_text("fixture")
        manifest = {p.relative_to(self.root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in self.bundle.rglob("*") if p.is_file()}
        (self.root / "INPUT_MANIFEST.json").write_text(json.dumps(manifest))
        root_patch = patch.object(runner, "ROOT", self.root)
        root_patch.start()
        self.addCleanup(root_patch.stop)

    def test_manifest_corruption(self):
        runner.verify_manifest()
        (self.bundle / "input0").write_text("modified")
        with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
            runner.verify_manifest()

    def test_missing_input(self):
        (self.bundle / "input0").unlink()
        with self.assertRaises(FileNotFoundError):
            runner.verify_manifest()

    def test_proof_mismatch(self):
        (self.root / "docs/r17_proof_zh.md").write_text("modified")
        with self.assertRaisesRegex(RuntimeError, "Standalone proof"):
            runner.verify_manifest()

    def test_unsafe_manifest(self):
        file = self.root / "INPUT_MANIFEST.json"
        manifest = json.loads(file.read_text())
        manifest["../outside"] = manifest.pop("proof_bundle/input0")
        file.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(RuntimeError, "Invalid manifest path"):
            runner.verify_manifest()

    def fake_run(self, command, *, cwd, check):
        self.assertEqual(command[-3:], ["-S", "-B", "verify_all.py"])
        self.assertTrue(check)
        self.work = Path(cwd)
        self.assertEqual(list((self.work / "reports").iterdir()), [])
        self.assertFalse(list((self.work / "core").glob("*_verified.json")))
        result = {"verified": True, "unresolved_leaves": 0, "kkt_dimension": 157,
                  "surviving_topology_orbits": 1344,
                  "steps": [{"step": s} for s in self.steps]}
        (self.work / "reports/MASTER_VERIFIED.json").write_text(json.dumps(result))

    def test_temporary_replay_export_and_cleanup(self):
        self.steps = runner.STEPS
        output = self.root / "fresh"
        with patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            runner.replay(output)
        self.assertFalse(self.work.exists())
        self.assertTrue((output / "MASTER_VERIFIED.json").exists())
        self.assertTrue((output / "core").is_dir())
        self.assertFalse((output / "legacy.log").exists())
        runner.verify_manifest()

    def test_skipped_accepting_stage(self):
        self.steps = runner.STEPS[:-1]
        with patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            with self.assertRaisesRegex(RuntimeError, "seven accepting stages"):
                runner.replay()
        self.assertFalse(self.work.exists())

    def test_stale_success_not_accepted(self):
        with patch.object(runner.subprocess, "run", return_value=None):
            with self.assertRaises(FileNotFoundError):
                runner.replay()

    def test_protected_output_paths(self):
        for target in [self.bundle / "new_reports", self.root / "docs"]:
            with self.assertRaises(ValueError):
                runner.replay(target)


if __name__ == "__main__":
    unittest.main()
