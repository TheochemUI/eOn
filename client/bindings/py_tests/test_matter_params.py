import pytest
import pyeonclient as ec
import sh

@pytest.fixture
def one_pt_morse_dimer(shared_datadir):
    ddir=f"{shared_datadir}/client/one_Pt_on_frozenSurface"
    data_dir=f"{shared_datadir}/client/jobs/dimer"
    sh.cp(f"{data_dir}/morse_dimer.ini",f"{ddir}")
    sh.cd(ddir)
    files = sh.ls()
    diff = set(files.split()) ^ {'morse_dimer.ini', 'displacement.con', 'direction.dat', 'pos.con'}
    assert not diff

def test_parameters():
    params = ec.Parameters()
    assert params.potential == "lj"
    params.potential = "morse"
    assert params.potential == "morse"

class TestMatter:

    def test_load_con(self, one_pt_morse_dimer):
        params = ec.Parameters()
        m1 = ec.Matter(params)
        m1.con2matter("pos.con")

    def test_get_potential(self, one_pt_morse_dimer):
        params = ec.Parameters()
        m1 = ec.Matter(params)
        assert m1.pot_energy == 0.
        m1.con2matter("pos.con")
        assert m1.pot_energy == -21.720527098776707
