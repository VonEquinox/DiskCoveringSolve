"""Small split-input and fail-closed wrapper regressions, not proof replay."""
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_verification as runner
import verify_interfaces_extra as extra


class RunnerTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="cover20_runner_test_")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        root_patch = patch.object(runner, "ROOT", self.root)
        root_patch.start()
        self.addCleanup(root_patch.stop)
        silence = patch("sys.stdout", new_callable=io.StringIO)
        silence.start()
        self.addCleanup(silence.stop)
        (self.root / "proof_bundle/core").mkdir(parents=True)
        (self.root / "proof_bundle/enumeration").mkdir()
        (self.root / "proof_parts").mkdir()
        (self.root / "docs").mkdir()
        self.inputs = {}
        for name in ["PROOF_zh.md"] + [f"core/input{i}" for i in range(49)]:
            rel = "proof_bundle/" + name
            (self.root / rel).write_bytes(b"fixture")
            self.inputs[rel] = hashlib.sha256(b"fixture").hexdigest()
        (self.root / "docs/r20_proof_zh.md").write_bytes(b"fixture")
        parts = []
        for i, block in enumerate((b"abc", b"def", b"ghi")):
            rel = f"proof_parts/part{i}"
            (self.root / rel).write_bytes(block)
            parts.append({"path": rel, "bytes": len(block), "sha256": hashlib.sha256(block).hexdigest()})
        self.inputs[runner.LARGE] = hashlib.sha256(b"abcdefghi").hexdigest()
        self.storage = {runner.LARGE: {"bytes": 9, "parts": parts}}
        self.save_manifests()

    def save_manifests(self):
        (self.root / "INPUT_MANIFEST.json").write_text(json.dumps(self.inputs))
        (self.root / "STORAGE_MANIFEST.json").write_text(json.dumps(self.storage))

    def result(self):
        return {"verified": True, "unresolved_leaves": 0,
                "steps": [{"step": s} for s in sorted(runner.STAGES)],
                "kkt_dimension": 198, "surviving_orbits": 562169, "rooted_survivors": 14747883,
                "noncandidate_cases": 562168, "candidate_cases": 1,
                "unresolved_combinatorial_branches": 0, "recurrence_regression_cases": 33,
                "anchor_isolation_radius": "11/100", "negative_tests": {"negative_test_count": 43},
                "candidate_tree": {"nodes": 6171, "counts": {"LOCAL": 102}},
                "noncandidate_forest": {"nodes": 657498}}

    def test_reconstruct_original_bytes(self):
        inputs, storage = runner.verify_manifest()
        runner.materialize(self.root / "fresh", inputs, storage)
        self.assertEqual((self.root / "fresh/enumeration/noncandidate20_trees.jsonl.gz").read_bytes(), b"abcdefghi")
        runner.verify_manifest()

    def test_changed_part(self):
        (self.root / "proof_parts/part0").write_bytes(b"bad")
        with self.assertRaisesRegex(RuntimeError, "Part hash"):
            runner.verify_manifest()

    def test_missing_part(self):
        (self.root / "proof_parts/part0").unlink()
        with self.assertRaises(FileNotFoundError):
            runner.verify_manifest()

    def test_reordered_parts(self):
        self.storage[runner.LARGE]["parts"].reverse()
        self.save_manifests()
        with self.assertRaisesRegex(RuntimeError, "Original input SHA"):
            runner.verify_manifest()

    def test_duplicate_part(self):
        self.storage[runner.LARGE]["parts"][1] = self.storage[runner.LARGE]["parts"][0]
        self.save_manifests()
        with self.assertRaisesRegex(RuntimeError, "repeated part"):
            runner.verify_manifest()

    def test_wrong_length(self):
        self.storage[runner.LARGE]["bytes"] = 10
        self.save_manifests()
        with self.assertRaisesRegex(RuntimeError, "Reconstructed length"):
            runner.verify_manifest()

    def test_unsafe_path(self):
        self.storage[runner.LARGE]["parts"][0]["path"] = "../escape"
        self.save_manifests()
        with self.assertRaises(ValueError):
            runner.verify_manifest()

    def test_unexpected_file(self):
        (self.root / "proof_parts/extra").write_text("extra")
        with self.assertRaisesRegex(RuntimeError, "physical input"):
            runner.verify_manifest()

    def test_changed_unsplit_input(self):
        (self.root / "proof_bundle/core/input0").write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "Original input SHA"):
            runner.verify_manifest()

    def test_published_proof_mismatch(self):
        (self.root / "docs/r20_proof_zh.md").write_text("wrong version")
        with self.assertRaisesRegex(RuntimeError, "Published proof"):
            runner.verify_manifest()

    def test_changed_during_reconstruction(self):
        inputs, storage = runner.verify_manifest()
        (self.root / "proof_parts/part0").write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "Part changed"):
            runner.materialize(self.root / "fresh", inputs, storage)

    def test_stage_completeness(self):
        result = self.result()
        runner.validate_result(result)
        result["steps"][-1] = result["steps"][0]
        with self.assertRaisesRegex(RuntimeError, "exactly once"):
            runner.validate_result(result)

    def test_unresolved_result(self):
        result = self.result()
        result["unresolved_leaves"] = 1
        with self.assertRaisesRegex(RuntimeError, "complete successful"):
            runner.validate_result(result)

    def test_missing_regressions(self):
        result = self.result()
        result["negative_tests"] = {}
        with self.assertRaisesRegex(RuntimeError, "Missing negative"):
            runner.validate_result(result)

    def test_protected_output(self):
        for path in (self.root, self.root / "proof_parts/new", self.root / "proof_bundle/new"):
            with self.assertRaises(ValueError):
                runner.report_target(path)

    def test_success_requires_new_report(self):
        with patch.object(runner.subprocess, "run", return_value=None):
            with self.assertRaises(FileNotFoundError):
                runner.replay(1)

    def test_replay_export_and_cleanup(self):
        paths = []
        def fake_run(args, *, cwd, env, check):
            self.assertEqual(args[-3:], ["-S", "-B", "verify_all.py"])
            self.assertNotIn("PYTHONOPTIMIZE", env)
            self.assertNotIn("PYTHONPATH", env)
            self.assertEqual(env["COVER20_WORKERS"], "2")
            self.assertTrue(check)
            paths.append(cwd)
            self.assertFalse((cwd / "reports").exists())
            (cwd / "reports").mkdir()
            (cwd / "reports/MASTER_VERIFIED.json").write_text(json.dumps(self.result()))
        with patch.object(runner.subprocess, "run", side_effect=fake_run), patch.object(extra, "verify", return_value={"verified": True}):
            runner.replay(2, self.root / "export")
        self.assertFalse(paths[0].exists())
        self.assertTrue((self.root / "export/independent_interfaces.json").exists())
        runner.verify_manifest()


if __name__ == "__main__":
    unittest.main()
