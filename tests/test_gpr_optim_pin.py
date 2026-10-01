"""The GP dimer pin is one commit, named in the wrap and the workflow, and
the private gpr_optim tree never enters the eOn repository."""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REV = "b55c89e2115388f901839aba2a5808bfcef06f68"


def _git(*args):
    return subprocess.run(
        ["git", "-C", str(ROOT), *args], capture_output=True, text=True
    )


def test_wrap_and_workflow_pin_gpr_optim_develop():
    wrap = (ROOT / "subprojects" / "gpr_optim.wrap").read_text()
    workflow = (ROOT / ".github" / "workflows" / "ci_build_gprd.yml").read_text()
    found = re.search(r"^revision = ([0-9a-f]{40})$", wrap, re.M)
    assert found is not None
    assert found.group(1) == REV
    assert REV in workflow


def test_private_gpr_optim_checkout_stays_untracked():
    # A local checkout of the private repository is ignored, and nothing
    # under subprojects/gpr_optim is tracked by eOn.
    ignored = _git("check-ignore", "-q", "subprojects/gpr_optim/meson.build")
    assert ignored.returncode == 0
    tracked = _git("ls-files", "--", "subprojects/gpr_optim")
    assert tracked.returncode == 0
    assert tracked.stdout.strip() == ""


def test_gp_dimer_workflow_does_not_use_a_private_key():
    workflow = (ROOT / ".github" / "workflows" / "ci_build_gprd.yml").read_text()
    assert "SUBMODULE_PRIVATE" not in workflow
    assert "ssh-private-key" not in workflow
    assert "AtomicGPDimer" in workflow
