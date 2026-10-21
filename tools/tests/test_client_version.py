"""Revision changes refresh C++ metadata without rebuilding unchanged stamps."""

import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location(
    "client_version", Path(__file__).resolve().parents[2] / "client/get_version.py")
version = importlib.util.module_from_spec(spec)
spec.loader.exec_module(version)

TEMPLATE = '''const std::string GIT_HASH = "@GIT_HASH@";
const std::string GIT_HASH_FULL = "@GIT_HASH_FULL@";
const std::string BUILD_DATE = "@BUILD_DATE@";
const std::string FEATURES_STRING = "MPI enabled";
'''


class VersionHeaderTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.head = self.root / "HEAD"
        self.head.write_text("a" * 40)
        self.template = self.root / "template.h.in"
        self.template.write_text(TEMPLATE)
        self.output = self.root / "version.h"
        self.git = patch.object(version.subprocess, "check_output", side_effect=self.read_head).start()
        self.addCleanup(patch.stopall)

    def read_head(self, command, **kwargs):
        self.assertEqual(command, ["git", "rev-parse", "HEAD"])
        self.assertEqual(kwargs["cwd"], self.root)
        return self.head.read_bytes() + b"\n"

    def generate(self):
        return version.write_header(self.template, self.output, self.root)

    def test_commit_change_refreshes_both_revision_fields(self):
        self.assertTrue(self.generate())
        self.assertIn('GIT_HASH = "aaaaaaaa"', self.output.read_text())
        self.head.write_text("b" * 40)
        self.assertTrue(self.generate())
        text = self.output.read_text()
        self.assertIn('GIT_HASH = "bbbbbbbb"', text)
        self.assertIn('GIT_HASH_FULL = "' + "b" * 40 + '"', text)
        self.assertNotIn("@GIT_HASH", text)

    def test_unchanged_commit_preserves_contents_and_mtime(self):
        self.generate()
        os.utime(self.output, ns=(1234567890000000000, 1234567890000000000))
        expected = self.output.read_bytes()
        with patch.object(version, "get_build_date", side_effect=AssertionError("unchanged header must retain date")):
            self.assertFalse(self.generate())
        self.assertEqual(self.output.read_bytes(), expected)
        self.assertEqual(self.output.stat().st_mtime_ns, 1234567890000000000)

    def test_feature_change_updates_header_with_same_commit(self):
        self.generate()
        self.template.write_text(TEMPLATE.replace("MPI enabled", "MPI disabled"))
        self.assertTrue(self.generate())
        self.assertIn("MPI disabled", self.output.read_text())

    def test_archive_without_git_has_explicit_unknown_revision(self):
        self.git.side_effect = subprocess.CalledProcessError(128, ["git"])
        self.assertTrue(self.generate())
        self.assertIn('GIT_HASH = "unknown"', self.output.read_text())
        self.assertIn('GIT_HASH_FULL = "unknown"', self.output.read_text())

    def test_invalid_revision_is_not_embedded(self):
        self.head.write_text('invalid"revision')
        self.generate()
        self.assertIn('GIT_HASH_FULL = "unknown"', self.output.read_text())


if __name__ == "__main__":
    unittest.main()
