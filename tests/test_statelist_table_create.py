"""A states directory without state_table still reports zero states."""

from types import SimpleNamespace

from eon.statelist import StateList


def test_resume_creates_a_missing_state_table(tmp_path):
    root = tmp_path / "states dir"
    (root / "0").mkdir(parents=True)
    cfg = SimpleNamespace(
        path_states=str(root),
        comp_eps_e=1.0e-3,
        comp_eps_r=0.1,
        comp_use_identical=False,
    )
    listing = StateList(object, config=cfg)
    assert listing.get_num_states() == 0
    assert (root / "state_table").read_text() == ""


def test_resume_leaves_an_existing_state_table(tmp_path):
    root = tmp_path / "states dir"
    (root / "0").mkdir(parents=True)
    table = root / "state_table"
    original = "      0        -1.25000\n"
    table.write_text(original)
    cfg = SimpleNamespace(
        path_states=str(root),
        comp_eps_e=1.0e-3,
        comp_eps_r=0.1,
        comp_use_identical=False,
    )
    listing = StateList(object, config=cfg)
    assert listing.get_num_states() == 1
    assert table.read_text() == original
