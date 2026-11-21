"""Run a native test with writable copies of its fixture inputs."""

import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


def main():
    executable, *args = sys.argv[1:]
    split = args.index("--")
    systems = Path(os.environ["EON_TEST_SYSTEMS_DIR"])
    with TemporaryDirectory() as directory:
        root = Path(directory)
        for name in args[:split]:
            dest = root / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(systems / name, dest)
        workdir = root / Path.cwd().relative_to(systems)
        env = dict(os.environ, TMPDIR=directory, TMP=directory, TEMP=directory)
        return subprocess.run([executable, *args[split + 1 :]], cwd=workdir, env=env).returncode


if __name__ == "__main__":
    sys.exit(main())
