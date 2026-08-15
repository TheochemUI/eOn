import sh
from pathlib import Path

p = Path(str(sh.pwd())) # Hacky way to get project root
# eonclient = sh.Command(str(p).strip()+"/client/build/eonclient")


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
        # Basic-eon is built without -Dwith_gprd; the request must refuse
        # rather than silently run ImprovedDimer (status 12).
        text = str(exc)
        if Path("client.log").is_file():
            text += Path("client.log").read_text()
        assert "with_gprd" in text
        return
    results = read_results(f"{ddir}/results.dat")
    assert results["termination_reason"] == "0"
    # magic_enum writes the PotType identifier exactly as declared.
    assert results["potential_type"] == "MORSE_PT"

# Broken AMS
# def test_one_pt_ams_dimer(datadir, shared_datadir, eonclient):
#     ddir=f"{shared_datadir}/one_Pt_on_frozenSurface"
#     sh.cp(f"{datadir}/ams_io_dimer.ini",f"{ddir}/config.ini")
#     sh.cd(ddir)
#     files = sh.ls()
#     diff = set(files.split()) ^ {'config.ini', 'displacement.con', 'direction.dat', 'pos.con'}
#     assert not diff
#     eonclient() # Runs eon
#     with open(f"{ddir}/results.dat", 'r') as res:
#         resText = res.readlines()
#         assert resText[0] == "0 termination_reason\n"
#         assert resText[3].split()[0] == "ams"
