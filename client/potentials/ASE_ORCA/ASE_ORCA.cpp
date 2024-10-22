/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/

#include "ASE_ORCA.h"
#include <fenv.h>

namespace eonc {

// XXX: This always assumes that charge is 0, mult is 1
// ASE default ----------------------------^ ---------^
// See also: https://gitlab.com/ase/ase/-/issues/1357
ASEOrcaPot::ASEOrcaPot(const ASEOrcaPot::Params &a_p) {
  // Fix for gh-184, see
  // https://github.com/numpy/numpy/issues/20504#issuecomment-985542508
  fenv_t orig_feenv;
  feholdexcept(&orig_feenv);
  _ase = py::module_::import("ase");
  fesetenv(&orig_feenv);
  py::module_ ase_orca = py::module_::import("ase.calculators.orca");
  py::module_ psutil = py::module_::import("psutil");

  // Set up ORCA profile and calculator
  py::object OrcaProfile = ase_orca.attr("OrcaProfile");
  py::object ORCA = ase_orca.attr("ORCA");
  size_t nproc = (a_p.orca_nproc == "auto")
                     ? py::cast<int>(psutil.attr("cpu_count")(false))
                     : std::stoi(a_p.orca_nproc);

  this->_calc =
      ORCA("profile"_a = OrcaProfile(py::str(a_p.orca_path)),
           "orcasimpleinput"_a = a_p.simpleinput,
           "orcablocks"_a = py::str(fmt::format("%pal nprocs {} end", nproc)),
           "directory"_a = ".");
}

void ASEOrcaPot::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  // TODO(rg) Test after type maps
  MatrixType positions = MatrixType::Map(fip.pos, fip.nAtoms, 3);
  MatrixType boxx = MatrixType::Map(fip.box, 3, 3);
  Vector<size_t> atmnmrs = Vector<size_t>::Map(fip.atmnrs, fip.nAtoms);
  py::object atoms = this->_ase.attr("Atoms")(
      "numbers"_a = atmnmrs, "positions"_a = positions, "cell"_a = boxx);
  atoms.attr("set_calculator")(this->_calc);
  atoms.attr("set_pbc")(std::tuple<bool, bool, bool>(true, true, true));
  double py_e = py::cast<double>(atoms.attr("get_potential_energy")());
  MatrixType py_force = py::cast<MatrixType>(atoms.attr("get_forces")());

  // Populate the output parameters
  efvd->energy = py_e;
  for (size_t idx = 0; idx < fip.nAtoms; ++idx) {
    efvd->F[3 * idx] = py_force(idx, 0);
    efvd->F[3 * idx + 1] = py_force(idx, 1);
    efvd->F[3 * idx + 2] = py_force(idx, 2);
  }
  return;
}

} // namespace eonc
