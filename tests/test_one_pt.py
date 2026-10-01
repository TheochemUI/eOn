import shutil
import sys

import pytest
import sh
from pathlib import Path

p = Path(str(sh.pwd())) # Hacky way to get project root
# eonclient = sh.Command(str(p).strip()+"/client/build/eonclient")

# Same saddle on this cell: a different Pt site is about an angstrom away,
# and both configs set the force tolerance at 0.001.
_SAME_SADDLE_ENERGY_EV = 1e-3
_SAME_SADDLE_POSITION_A = 1e-2


def read_results(path):
    """results.dat holds one `<value> <key>` record per line.

    Keys are read by name because their order is a property of the job that
    wrote them: SaddleSearchJob emits `job_type` between the termination text
    and the seed, so a positional index means a different field per job.
    """
    records = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                records[fields[-1]] = " ".join(fields[:-1])
    return records

def test_one_pt_morse_dimer(datadir, shared_datadir, eonclient, monkeypatch):
    ddir=f"{shared_datadir}/client/one_Pt_on_frozenSurface"
    sh.cp(f"{datadir}/morse_dimer.ini",f"{ddir}/config.ini")
    monkeypatch.chdir(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'config.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff
    eonclient() # Runs eon
    results = read_results(f"{ddir}/results.dat")
    assert results["termination_reason"] == "0"
    # magic_enum writes the PotType identifier exactly as declared.
    assert results["potential_type"] == "MORSE_PT"

def test_one_pt_morse_gprdimer(datadir, shared_datadir, eonclient, monkeypatch):
    ddir=f"{shared_datadir}/client/one_Pt_on_frozenSurface"
    sh.cp(f"{datadir}/morse_gprdimer.ini",f"{ddir}/config.ini")
    monkeypatch.chdir(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'config.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff
    try:
        eonclient()
    except sh.ErrorReturnCode as exc:
        # A build with -Dwith_gprd=disabled must refuse gprdimer rather
        # than silently run ImprovedDimer (status 12).
        text = str(exc)
        if Path("client.log").is_file():
            text += Path("client.log").read_text()
        assert "with_gprd" in text
        return
    results = read_results(f"{ddir}/results.dat")
    assert results["termination_reason"] == "0"
    # magic_enum writes the PotType identifier exactly as declared.
    assert results["potential_type"] == "MORSE_PT"


def _con_positions(path):
    """Coordinate rows from a classic con or from a symbol/x/y/z frame."""
    lines = Path(path).read_text(errors="replace").splitlines()
    rows = []
    if any(line.startswith("Coordinates of component") for line in lines):
        started = False
        for line in lines:
            if line.startswith("Coordinates of component"):
                started = True
                continue
            if not started:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append(tuple(float(x) for x in parts[:3]))
            except ValueError:
                continue
        return rows
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            try:
                rows.append(tuple(float(x) for x in parts[1:4]))
                continue
            except ValueError:
                pass
        if len(parts) >= 3:
            try:
                rows.append(tuple(float(x) for x in parts[:3]))
            except ValueError:
                continue
    return rows


def _stage_morse_pt(shared_datadir, datadir, root, name, ini_name):
    src = Path(shared_datadir) / "client" / "one_Pt_on_frozenSurface"
    dst = Path(root) / name
    shutil.copytree(src, dst)
    shutil.copy(Path(datadir) / ini_name, dst / "config.ini")
    return dst


def _run_eonclient(eonclient, workdir):
    try:
        eonclient(_cwd=str(workdir))
    except sh.ErrorReturnCode as exc:
        text = str(exc)
        log = Path(workdir) / "client.log"
        if log.is_file():
            text += log.read_text(errors="replace")
        return None, text
    return read_results(Path(workdir) / "results.dat"), ""


def test_morse_pt_gprdimer_matches_dimer_fewer_calls(
    datadir, shared_datadir, eonclient, tmp_path
):
    """gprdimer on the Morse Pt cell matches the dimer saddle in fewer calls."""
    dimer_dir = _stage_morse_pt(
        shared_datadir, datadir, tmp_path, "dimer", "morse_dimer.ini"
    )
    gpr_dir = _stage_morse_pt(
        shared_datadir, datadir, tmp_path, "gprdimer", "morse_gprdimer.ini"
    )
    dimer, dimer_err = _run_eonclient(eonclient, dimer_dir)
    gpr, gpr_err = _run_eonclient(eonclient, gpr_dir)
    assert dimer is not None, dimer_err[-2000:]
    if gpr is None and "with_gprd" in gpr_err:
        # gpr_optim is private: the GP dimer is built only where a checkout
        # sits in subprojects/gpr_optim.
        checkout = (
            Path(__file__).resolve().parents[1]
            / "subprojects"
            / "gpr_optim"
            / "meson.build"
        )
        if sys.platform != "linux" or not checkout.is_file():
            pytest.skip(
                "GP dimer needs a Linux build with the private gpr_optim "
                "checkout in subprojects/gpr_optim (-Dwith_gprd=auto)"
            )
        pytest.fail(
            "linux eonclient refused min_mode_method=gprdimer although "
            "subprojects/gpr_optim is checked out. -Dwith_gprd=auto should "
            "link it. " + gpr_err[-2000:]
        )
    assert gpr is not None, gpr_err[-2000:]
    assert dimer["termination_reason"] == "0", dimer
    assert gpr["termination_reason"] == "0", gpr
    assert dimer["potential_type"] == "MORSE_PT"
    assert gpr["potential_type"] == "MORSE_PT"
    dimer_energy = float(dimer["potential_energy_saddle"])
    gpr_energy = float(gpr["potential_energy_saddle"])
    assert abs(dimer_energy - gpr_energy) < _SAME_SADDLE_ENERGY_EV, (
        dimer_energy,
        gpr_energy,
    )
    dimer_pos = _con_positions(dimer_dir / "saddle.con")
    gpr_pos = _con_positions(gpr_dir / "saddle.con")
    assert len(dimer_pos) == len(gpr_pos) > 0
    max_dr = max(
        max(abs(a - b) for a, b in zip(left, right))
        for left, right in zip(dimer_pos, gpr_pos)
    )
    assert max_dr < _SAME_SADDLE_POSITION_A, max_dr
    dimer_calls = int(float(dimer["total_force_calls"]))
    gpr_calls = int(float(gpr["total_force_calls"]))
    assert gpr_calls > 0
    assert gpr_calls < dimer_calls, (gpr_calls, dimer_calls)
