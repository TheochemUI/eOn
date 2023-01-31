import pytest
import pyeonclient as ec
import sh

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
