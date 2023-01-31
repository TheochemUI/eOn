import io
import contextlib

import pytest
import pyeonclient as ec
import sh

class TestSaddleSearchJob:
    def test_saddlesearchjob(shared_datadir):
        params = ec.Parameters()
        ec.log_init(params, "blah.log")
        params.load("morse_dimer.ini")
        ssj = ec.SaddleSearchJob(params)
        assert ssj.run() == ['results.dat', 'mode.dat', 'saddle.con']
        # Until the logging library is updated, no results can be parsed
