import pytest
import sh

@pytest.fixture(autouse=True)
def one_pt_morse_dimer(shared_datadir):
    ddir=f"{shared_datadir}/client/one_Pt_on_frozenSurface"
    data_dir=f"{shared_datadir}/client/jobs/dimer"
    sh.cp(f"{data_dir}/morse_dimer.ini",f"{ddir}")
    sh.cd(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'morse_dimer.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff
