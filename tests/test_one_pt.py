import sh
import pytest
from pathlib import Path

p = Path(str(sh.pwd())) # Hacky way to get project root
# eonclient = sh.Command(str(p).strip()+"/client/build/eonclient")
eonclient = sh.eonclient

def test_one_pt_morse_dimer(datadir,shared_datadir):
    ddir=f"{shared_datadir}/client/one_Pt_on_frozenSurface"
    sh.cp(f"{datadir}/morse_dimer.ini",f"{ddir}/config.ini")
    sh.cd(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'config.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff
    eonclient() # Runs eon
    with open(f"{ddir}/results.dat", 'r') as res:
        resText = res.readlines()
        assert resText[0] == "0 termination_reason\n"
        assert resText[3].split()[0] == "morse_pt"

def test_one_pt_morse_gprdimer(datadir,shared_datadir):
    ddir=f"{shared_datadir}/client/one_Pt_on_frozenSurface"
    sh.cp(f"{datadir}/morse_gprdimer.ini",f"{ddir}/config.ini")
    sh.cd(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'config.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff
    eonclient() # Runs eon
    with open(f"{ddir}/results.dat", 'r') as res:
        resText = res.readlines()
        assert resText[0] == "0 termination_reason\n"
        assert resText[3].split()[0] == "morse_pt"

# Broken AMS
# def test_one_pt_ams_dimer(datadir,shared_datadir):
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
