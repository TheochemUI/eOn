import sh
import os
import pytest
from pathlib import Path

p = Path(str(sh.pwd())) # Hacky way to get project root
new_env = os.environ.copy()
new_env["PATH"] += os.pathsep + str(p).strip()+"/bin"
new_env["PYTHONPATH"] += os.pathsep + str(p).strip()+"/eon"
eonclient = sh.Command(str(p).strip()+"/client/build/eonclient")
eon = sh.Command(str(p).strip()+"/bin/eon")

def test_akmc_pt_morse_dimer(datadir,shared_datadir):
    ddir=f"{shared_datadir}/server/Pt_Heptamer_oneLayer"
    sh.cp(f"{datadir}/morse_dimer.ini",f"{ddir}/config.ini")
    sh.cd(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'config.ini', 'pos.con'}
    assert not diff
    eon(_env=new_env) # Runs eon
    eon(_env=new_env) # Runs eon
    eon(_env=new_env) # Runs eon
    # TODO: Actually parse an output
