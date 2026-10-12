"""Loud skip for the Si3N4 point when the kit repo cannot be read.

The OpenCPMD workflow runs this file only after the token check fails.
The skip reason is the record that the energy comparison did not run.
"""

import pytest


def test_si3n4_point_skipped():
    pytest.skip(
        "TheochemUI/sige_repro is not readable by this token, so the Si3N4 point is skipped."
    )
