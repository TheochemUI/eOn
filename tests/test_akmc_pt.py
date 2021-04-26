import sh
import os
import pytest
from pathlib import Path

p = Path(str(sh.pwd())) # Hacky way to get project root
eonclient = sh.Command(str(p).strip()+"/client/build/eonclient")
eon = sh.Command(str(p).strip()+"/bin/eon")
os.environ['PATH'] += os.pathsep + str(p).strip()+"/bin"

def test_akmc_pt_morse_dimer(datadir,shared_datadir):
    ddir=f"{shared_datadir}/server/one_Pt_oneLayer"
    sh.cp(f"{datadir}/morse_dimer.ini",f"{ddir}/config.ini")
    sh.cd(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'config.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff
    eon() # Runs eon
    eon() # Runs eon
    eon() # Runs eon
