"""Fast runner tests; these do not replace the full mathematical replay."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import run_verification as runner


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.output = patch("sys.stdout", new_callable=io.StringIO)
        self.output.start()
        self.addCleanup(self.output.stop)
        self.temp = tempfile.TemporaryDirectory(prefix="cover16_runner_test_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle = self.root / "proof_bundle"
        (self.bundle / "reports").mkdir(parents=True)
        (self.root / "docs").mkdir()
        (self.bundle / "PROOF_zh.md").write_text("Fixture proof\n")
        (self.root / "docs/r16_proof_zh.md").write_text("Fixture proof\n")
        (self.bundle / "reports/MASTER_VERIFIED.json").write_text('{"verified": true}')
        for i in range(88):
            (self.bundle / f"input{i}.txt").write_text(str(i))
        manifest = {
            p.relative_to(self.root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in self.bundle.rglob("*") if p.is_file()
        }
        (self.root / "INPUT_MANIFEST.json").write_text(json.dumps(manifest))
        self.patch_root = patch.object(runner, "ROOT", self.root)
        self.patch_root.start()
        self.addCleanup(self.patch_root.stop)

    def test_manifest_accepts_original_and_rejects_modified_input(self):
        runner.verify_manifest()
        (self.bundle / "input0.txt").write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
            runner.verify_manifest()

    def test_missing_input_rejected(self):
        (self.bundle / "input0.txt").unlink()
        with self.assertRaises(FileNotFoundError):
            runner.verify_manifest()

    def test_manuscript_mismatch_rejected(self):
        (self.root / "docs/r16_proof_zh.md").write_text("different")
        with self.assertRaisesRegex(RuntimeError, "Standalone manuscript"):
            runner.verify_manifest()

    def test_unsafe_manifest_path_rejected(self):
        p = self.root / "INPUT_MANIFEST.json"
        manifest = json.loads(p.read_text())
        manifest["../outside.txt"] = manifest.pop("proof_bundle/input0.txt")
        p.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(RuntimeError, "Invalid manifest path"):
            runner.verify_manifest()

    def fake_run(self, command, *, cwd, env, check):
        self.assertEqual(command[-2:], ["--jobs", "2"])
        self.assertNotIn("--reuse-enumeration", command)
        self.assertTrue(check)
        self.work = Path(cwd)
        self.assertNotEqual(self.work, self.bundle)
        master = self.work / "reports/MASTER_VERIFIED.json"
        self.assertFalse(master.exists(), "Submitted success report must be removed")
        master.write_text(json.dumps({"verified": True, "unresolved_leaves": 0,
                                     "enumeration_mode": self.mode}))
        return subprocess.CompletedProcess(command, 0)

    def test_temporary_copy_report_export_and_cleanup(self):
        self.mode = "regenerated_and_independently_audited"
        output = self.root / "fresh_reports"
        with patch.object(runner, "environment", return_value={}), \
                patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            runner.replay(2, output)
        self.assertFalse(self.work.exists())
        self.assertTrue((output / "MASTER_VERIFIED.json").exists())
        runner.verify_manifest()

    def test_reused_census_rejected_and_cleaned(self):
        self.mode = "reused_untrusted_tables_independently_audited"
        with patch.object(runner, "environment", return_value={}), \
                patch.object(runner.subprocess, "run", side_effect=self.fake_run):
            with self.assertRaisesRegex(RuntimeError, "regeneration"):
                runner.replay(2)
        self.assertFalse(self.work.exists())

    def test_missing_fresh_success_cannot_use_old_report(self):
        with patch.object(runner, "environment", return_value={}), \
                patch.object(runner.subprocess, "run", return_value=None):
            with self.assertRaises(FileNotFoundError):
                runner.replay(1)

    def test_archive_or_existing_export_target_rejected(self):
        for target in [self.bundle / "new_reports", self.root / "docs"]:
            with self.assertRaises(ValueError):
                runner.replay(1, target)


if __name__ == "__main__":
    unittest.main()
