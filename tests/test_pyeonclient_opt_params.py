"""pyeonclient can set optimizer / refine / NEB methods (TheochemUI/eOn#406)."""

from __future__ import annotations

import pytest

pyec = pytest.importorskip("pyeonclient")


def test_opt_method_and_refine_roundtrip():
    params = pyec.Parameters()
    assert params.opt_method == pyec.OptType.CG
    params.opt_method = pyec.OptType.FIRE
    assert params.opt_method == pyec.OptType.FIRE
    params.refine_opt_method = pyec.OptType.LBFGS
    params.refine_threshold = 0.5
    assert params.refine_opt_method == pyec.OptType.LBFGS
    assert params.refine_threshold == pytest.approx(0.5)
    params.neb_opt_method = pyec.OptType.FIRE
    assert params.neb_opt_method == pyec.OptType.FIRE
    params.refine_opt_method = pyec.OptType.None_
    assert params.refine_opt_method == pyec.OptType.None_
    assert params.dimer_opt_method == pyec.OptType.CG
    params.dimer_opt_method = pyec.OptType.LBFGS
    assert params.dimer_opt_method == pyec.OptType.LBFGS


def test_dftd_pot_types_and_params():
    assert pyec.PotType.DFTD3 != pyec.PotType.LJ
    assert pyec.PotType.DFTD4 != pyec.PotType.DFTD3
    params = pyec.Parameters()
    params.potential = pyec.PotType.DFTD3
    params.dftd_functional = "pbe"
    params.dftd_atm = True
    params.d3_damping = "bj"
    assert params.potential == pyec.PotType.DFTD3
    assert params.dftd_functional == "pbe"
    assert params.d3_damping == "bj"


def test_expr_pot_lj_plus_morse():
    params = pyec.Parameters()
    params.potential = pyec.PotType.EXPR
    params.expr_expression = "0.5*lj + morse"
    params.expr_terms = "lj,morse"
    pot = pyec.make_potential(params.potential, params)
    assert pot is not None


def test_mopac_pot_type_and_params():
    params = pyec.Parameters()
    params.potential = pyec.PotType.MOPAC
    params.mopac_model = 4
    params.mopac_charge = 0
    assert params.potential == pyec.PotType.MOPAC
    assert params.mopac_model == 4
