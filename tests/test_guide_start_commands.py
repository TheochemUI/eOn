"""Guide commands, client paths, and sample keys that start a run."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UTILS = ROOT / "eon" / "_utils.py"
MESON = ROOT / "client" / "meson.build"
TROUBLE = (
    ROOT / "docs" / "source" / "tutorials" / "CECAM_LTS_MAP_2024" / "troubleshooting.md"
)
DISPLACE = ROOT / "docs" / "source" / "tutorials" / "displacement_scripts.md"
PARREP = ROOT / "docs" / "source" / "tutorials" / "parrep.md"
VIZ = ROOT / "docs" / "source" / "tutorials" / "visualization.md"
BASIN = ROOT / "examples" / "basin-hopping" / "config.ini"
REPLICA = ROOT / "examples" / "parallel-replica" / ".config.ini.tur"


def test_displacement_guide_matches_the_server_argv():
    source = UTILS.read_text(encoding="utf-8")
    assert "[sys.executable, str(sconf.script_path), tmpf.name]" in source
    page = DISPLACE.read_text(encoding="utf-8")
    assert "only argument" in page
    assert "uvx ptmdisp.py" not in page
    assert "uvx adsorbate_region.py" not in page
    assert "python ptmdisp.py pos.con" in page


def test_install_and_replica_commands_start():
    trouble = TROUBLE.read_text(encoding="utf-8")
    assert "meson install -C bbdir" in trouble
    assert not re.search(r"(?m)^meson install bbdir\s*$", trouble)
    parrep = PARREP.read_text(encoding="utf-8")
    assert "for i in {0..2}; do python -m eon.server; done" in parrep
    assert "states/1/reactant.con/" not in parrep
    assert "`states/1/reactant.con`" in parrep


def test_samples_use_the_client_binary_and_live_keys():
    meson = MESON.read_text(encoding="utf-8")
    assert "executable(\n    'eonclient'," in meson
    for path in (BASIN, REPLICA):
        text = path.read_text(encoding="utf-8")
        assert "eonclient" in text
        assert "../../client/" not in text
    viz = VIZ.read_text(encoding="utf-8")
    assert "ci_mmf_penalty_strength" not in viz
    assert "ci_mmf_penalty_base" not in viz
