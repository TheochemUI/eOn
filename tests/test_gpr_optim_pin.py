"""The GP dimer checkout stays untracked, and eOn does not carry its wrap."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _git(*args):
    return subprocess.run(
        ["git", "-C", str(ROOT), *args], capture_output=True, text=True
    )


def test_gpr_optim_wrap_is_not_in_the_tree():
    assert not (ROOT / "subprojects" / "gpr_optim.wrap").exists()
    tracked = _git(
        "ls-files",
        "--",
        "subprojects/gpr_optim.wrap",
        "subprojects/rgmin.wrap",
        "subprojects/anneal.wrap",
        "subprojects/gpr_optim",
    )
    assert tracked.returncode == 0
    assert tracked.stdout.strip() == ""
    for path in (
        "subprojects/gpr_optim.wrap",
        "subprojects/rgmin.wrap",
        "subprojects/anneal.wrap",
        "subprojects/gpr_optim/meson.build",
    ):
        ignored = _git("check-ignore", "-q", path)
        assert ignored.returncode == 0, path
    workflow = (ROOT / ".github" / "workflows" / "ci_build_gprd.yml").read_text()
    assert "gpr_optim.wrap" not in workflow.replace(
        "test ! -e subprojects/gpr_optim.wrap", ""
    )


def test_gp_dimer_workflow_does_not_use_a_private_key():
    workflow = (ROOT / ".github" / "workflows" / "ci_build_gprd.yml").read_text()
    assert "SUBMODULE_PRIVATE" not in workflow
    assert "ssh-private-key" not in workflow
    assert "AtomicGPDimer" in workflow
