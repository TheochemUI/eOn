"""Exercise Meson's real interpreter expansion without compiling C++ sources."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class MesonVersionHeaderTest(unittest.TestCase):
    def test_real_custom_target_generates_and_preserves_revision_header(self):
        source = os.environ.get("EON_VERSION_SOURCE")
        if source is None:
            source = Path(__file__).resolve().parents[2]
        source = Path(source)
        helper = (source / "client/get_version.py").resolve()
        revision = subprocess.check_output(
            ["git", "-C", str(helper.parent), "rev-parse", "HEAD"], text=True).strip()
        git_dir = subprocess.check_output(
            ["git", "-C", str(helper.parent), "rev-parse", "--absolute-git-dir"], text=True).strip()
        environment = dict(os.environ, GIT_DIR=git_dir)
        meson = (source / "client/meson.build").read_text()
        start = meson.index("version_header = custom_target(")
        stop = meson.index("version_dependency =", start)
        target = meson[start:stop]
        with tempfile.TemporaryDirectory(prefix="eon-version-meson-") as directory:
            project = Path(directory)
            # A non-executable .py program requires Meson's interpreter prefix.
            # The isolated copy reads the source repository's Git metadata.
            shutil.copyfile(helper, project / "get_version.py")
            (project / "get_version.py").chmod(0o644)
            template = (source / "client/version.h.in").read_text()
            for key, value in {
                "VERSION": "3.6.0", "VERSION_MAJOR": "3", "VERSION_MINOR": "6",
                "VERSION_PATCH": "0", "OS_INFO": "test", "ARCH": "test",
                "FEATURES": "  MPI: enabled\\n",
            }.items():
                template = template.replace("@" + key + "@", value)
            (project / "template.h.in").write_text(template)
            (project / "meson.build").write_text(
                "project('version-header-fixture', version: '1')\n"
                "py = import('python').find_installation()\n"
                "version_script = find_program('get_version.py')\n"
                "version_template = files('template.h.in')\n" + target)
            subprocess.run(["meson", "setup", str(project / "build"), str(project)], env=environment, check=True)
            subprocess.run(["meson", "compile", "-C", str(project / "build")], env=environment, check=True)
            header = project / "build/version.h"
            text = header.read_text()
            self.assertIn('GIT_HASH_FULL = "' + revision + '"', text)
            self.assertIn('GIT_HASH = "' + revision[:8] + '"', text)
            self.assertIn("MPI: enabled", text)
            self.assertNotIn("@BUILD_DATE@", text)
            os.utime(header, ns=(1234567890000000000, 1234567890000000000))
            subprocess.run(["meson", "compile", "-C", str(project / "build")], env=environment, check=True)
            self.assertEqual(header.read_text(), text)
            self.assertEqual(header.stat().st_mtime_ns, 1234567890000000000)


if __name__ == "__main__":
    unittest.main()
