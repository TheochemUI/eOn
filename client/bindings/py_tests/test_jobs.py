import io
import contextlib

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


class TestSaddleSearchJob:
    def test_saddlesearchjob(shared_datadir):
        params = ec.Parameters()
        ec.log_init(params, "blah.log")
        params.load("morse_dimer.ini")
        ssj = ec.SaddleSearchJob(params)
        assert ssj.run() == ['results.dat', 'mode.dat', 'saddle.con']
        # Until the logging library is updated, no results can be parsed
