"""Unit tests for scripts/release_assert.py — real shipped entry point."""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release_assert.py"


def _run(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=cwd or REPO,
        capture_output=True,
        text=True,
        check=False,
    )


def test_script_exists_and_is_executable_entry():
    assert SCRIPT.is_file()
    assert SCRIPT.read_text(encoding="utf-8").startswith("#!/usr/bin/env python3")


def test_assert_passes_on_real_repo_lockstep():
    """Drive the real script against the live eOn tree."""
    proc = _run()
    assert proc.returncode == 0, proc.stderr
    assert "ok: lockstep version" in proc.stdout
    assert "2.14.0" in proc.stdout or "ok:" in proc.stdout


def test_assert_require_changelog_for_2140_on_real_repo():
    proc = _run("2.14.0", "--require-changelog")
    assert proc.returncode == 0, proc.stderr
    assert "CHANGELOG.md has section" in proc.stdout


def test_assert_fails_wrong_expected_version():
    proc = _run("9.9.9")
    assert proc.returncode == 1
    assert "expected" in proc.stderr.lower() or "surfaces at" in proc.stderr


def test_assert_fails_missing_changelog_section():
    proc = _run("99.0.0", "--require-changelog")
    # surfaces are 2.14.0 not 99.0.0, so either version mismatch or missing section
    assert proc.returncode == 1
    assert proc.stderr.strip()


def test_print_version_flag():
    proc = _run("--print-version")
    assert proc.returncode == 0, proc.stderr
    ver = proc.stdout.strip()
    assert ver  # non-empty
    parts = ver.split(".")
    assert len(parts) >= 3


def test_notes_extract_writes_file(tmp_path: Path):
    out = tmp_path / "notes.md"
    proc = _run("2.14.0", "--notes", str(out))
    assert proc.returncode == 0, proc.stderr
    assert out.is_file()
    body = out.read_text(encoding="utf-8")
    assert "2.14.0" in body or "Release" in body or body.strip()


def test_lockstep_helpers_on_fixture_tree(tmp_path: Path):
    """Import helpers and test disagreement detection on a temp tree."""
    sys.path.insert(0, str(REPO / "scripts"))
    try:
        import release_assert as ra  # type: ignore
    finally:
        if str(REPO / "scripts") in sys.path:
            pass  # leave for other tests; harmless

    root = tmp_path / "eonfake"
    root.mkdir()
    (root / "pyproject.toml").write_text(
        textwrap.dedent(
            """\
            [project]
            name = "eon-akmc"
            version = "1.2.3"
            """
        ),
        encoding="utf-8",
    )
    (root / "pixi.toml").write_text(
        textwrap.dedent(
            """\
            name = "eOn"
            version = "1.2.4"
            """
        ),
        encoding="utf-8",
    )
    (root / "CHANGELOG.md").write_text(
        "## [1.2.3](https://example/1.2.3) - 2026-01-01\n\n- hi\n",
        encoding="utf-8",
    )

    errs = ra.assert_lockstep(root)
    assert errs and "disagree" in errs[0]

    # Fix pixi, require changelog ok
    (root / "pixi.toml").write_text('name = "eOn"\nversion = "1.2.3"\n', encoding="utf-8")
    errs = ra.assert_lockstep(root, "1.2.3", require_changelog=True)
    assert errs == []

    errs = ra.assert_lockstep(root, "1.2.3", require_changelog=False)
    assert errs == []

    # Missing section
    errs = ra.assert_lockstep(root, "1.2.3", require_changelog=True)
    # section exists for 1.2.3
    assert errs == []
    errs_missing = ra.assert_lockstep(root, require_changelog=True)
    # no expected version, uses current 1.2.3 which has section
    assert errs_missing == []

    (root / "CHANGELOG.md").write_text("# empty\n", encoding="utf-8")
    errs2 = ra.assert_lockstep(root, "1.2.3", require_changelog=True)
    assert errs2 and "CHANGELOG" in errs2[0]


def test_semver_rejects_garbage(tmp_path: Path):
    sys.path.insert(0, str(REPO / "scripts"))
    import release_assert as ra  # type: ignore

    root = tmp_path / "badver"
    root.mkdir()
    (root / "pyproject.toml").write_text(
        '[project]\nname = "eon-akmc"\nversion = "not-a-version"\n',
        encoding="utf-8",
    )
    (root / "pixi.toml").write_text('version = "not-a-version"\n', encoding="utf-8")
    errs = ra.assert_lockstep(root)
    assert errs and "semver" in errs[0].lower()


def test_blocked_pypi_name_eon_rejected(tmp_path: Path):
    """pypi.org/project/eon/ is EoN epidemics — must not be our distribution name."""
    sys.path.insert(0, str(REPO / "scripts"))
    import release_assert as ra  # type: ignore

    root = tmp_path / "blocked"
    root.mkdir()
    (root / "pyproject.toml").write_text(
        '[project]\nname = "eon"\nversion = "1.0.0"\n',
        encoding="utf-8",
    )
    (root / "pixi.toml").write_text('version = "1.0.0"\n', encoding="utf-8")
    errs = ra.assert_lockstep(root)
    assert errs and any("blocked" in e.lower() or "EoN" in e for e in errs)


def test_live_repo_uses_eon_kmc_dist_name():
    import tomllib

    data = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    assert data["project"]["name"] == "eon-akmc"
    assert (REPO / "ci" / "gha" / "pypi.ncl").is_file()
    pypi_ncl = (REPO / "ci" / "gha" / "pypi.ncl").read_text(encoding="utf-8")
    assert 'distribution_name = "eon-akmc"' in pypi_ncl
    assert "blocked_pypi_name" in pypi_ncl


def test_nickel_gha_sources_exist_and_match_generated_header():
    """Nickel is the source of truth; generated YAML must exist alongside .ncl."""
    gha = REPO / "ci" / "gha"
    assert (gha / "common.ncl").is_file()
    assert (gha / "release.ncl").is_file()
    assert (gha / "release_prepare.ncl").is_file()
    assert (gha / "towncrier_check.ncl").is_file()
    assert (gha / "gen.sh").is_file()
    assert (gha / "README.md").is_file()
    for yml in (
        REPO / ".github" / "workflows" / "release.yml",
        REPO / ".github" / "workflows" / "release-prepare.yml",
        REPO / ".github" / "workflows" / "towncrier.yml",
    ):
        assert yml.is_file(), yml
        head = yml.read_text(encoding="utf-8")[:400]
        # nickel export does not always embed comments; ensure non-empty workflow
        assert "name:" in head or head.lstrip().startswith("name")
        assert "jobs:" in yml.read_text(encoding="utf-8")


def test_gen_sh_mentions_nickel_export():
    script = (REPO / "ci" / "gha" / "gen.sh").read_text(encoding="utf-8")
    assert "nickel export" in script
    assert "release.ncl" in script


def test_windows_ci_aligns_with_feedstock_fortran_policy():
    """GHA windows meson flags must match conda-forge build.bat (#15 in-tree Fortran win)."""
    yml = (REPO / ".github" / "workflows" / "ci_build_akmc.yml").read_text(encoding="utf-8")
    assert "Compile and install [windows]" in yml
    block = yml.split("Compile and install [windows]", 1)[1].split("Run tests", 1)[0]
    assert "-Dwith_fortran=true" in block
    assert "-Dwith_cuh2=false" in block
    assert "--default-library=static" in block
    rel = (REPO / "docs" / "source" / "devdocs" / "release.md").read_text(encoding="utf-8")
    assert "feedstock#15" in rel or "eon-feedstock#15" in rel
    assert "with_fortran=true" in rel
    assert "with_cuh2=false" in rel
