"""pyeonclient wheels include a win-64 metatomic build."""

from pathlib import Path

WORKFLOW = (
    Path(__file__).resolve().parents[1]
    / ".github"
    / "workflows"
    / "pyeonclient-wheels.yml"
)


def test_win64_metatomic_wheel_is_in_the_matrix():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "name: abi3-win\n" in text
    assert "name: abi3-win-metatomic\n" in text
    assert "os: windows-2022" in text
    assert 'cibw_build: "cp312-win_amd64"' in text
    assert "name: abi3-manylinux-metatomic" in text
    assert "CIBW_BEFORE_BUILD_LINUX" in text
    assert "CIBW_BEFORE_BUILD_WINDOWS" in text
    win = text.split("name: abi3-win-metatomic", 1)[1].split("steps:", 1)[0]
    assert "variant: metatomic" in win
    assert "windows-2022" in win
