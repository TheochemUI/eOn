"""The GP dimer pin is one commit, named in the wrap and the workflow."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REV = "b55c89e2115388f901839aba2a5808bfcef06f68"


def test_wrap_workflow_and_tree_pin_gpr_optim_develop():
    wrap = (ROOT / "subprojects" / "gpr_optim.wrap").read_text()
    workflow = (ROOT / ".github" / "workflows" / "ci_build_gprd.yml").read_text()
    found = re.search(r"^revision = ([0-9a-f]{40})$", wrap, re.M)
    assert found is not None
    assert found.group(1) == REV
    assert REV in workflow
    assert (ROOT / "subprojects" / "gpr_optim" / "REVISION").read_text().strip() == REV
    assert (ROOT / "subprojects" / "gpr_optim" / "meson.build").is_file()


def test_gp_dimer_workflow_does_not_use_a_private_key():
    workflow = (ROOT / ".github" / "workflows" / "ci_build_gprd.yml").read_text()
    assert "SUBMODULE_PRIVATE" not in workflow
    assert "ssh-private-key" not in workflow
    assert "AtomicGPDimer" in workflow
