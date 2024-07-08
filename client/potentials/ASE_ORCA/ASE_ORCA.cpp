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
#include "../../EnvHelpers.hpp"

ASEOrcaPot::ASEOrcaPot(std::shared_ptr<Parameters> a_params)
    : Potential(PotType::ASE_ORCA, a_params) {
  counter = 1;
  py::module_ sys = py::module_::import("sys");
  ase = py::module_::import("ase");
  py::module_ ase_orca = py::module_::import("ase.calculators.orca");
  py::module_ psutil = py::module_::import("psutil");
  std::string orcpth = helper_functions::get_value_from_env_or_param(
      "ORCA_COMMAND", a_params->aseorca.orca_path, "", "", true);
  std::string orca_simpleinput = helper_functions::get_value_from_env_or_param(
      "ORCA_SIMPLEINPUT", a_params->aseorca.simpleinput, "ENGRAD HF-3c",
      "Using ENGRAD HF-3c as a default input, set simpleinput or the "
      "environment variable ORCA_SIMPLEINPUT.\n");

  // Set up ORCA profile and calculator
  py::object OrcaProfile = ase_orca.attr("OrcaProfile");
  py::object ORCA = ase_orca.attr("ORCA");
  size_t nproc{0};

  if (a_params->aseorca.orca_nproc == "auto") {
    nproc = py::cast<int>(psutil.attr("cpu_count")(false));
  } else {
    nproc = std::stoi(a_params->aseorca.orca_nproc);
  }

  this->calc =
      ORCA("profile"_a = OrcaProfile(py::str(orcpth)),
           "orcasimpleinput"_a = orca_simpleinput,
           "orcablocks"_a = py::str(fmt::format("%pal nprocs {} end", nproc)),
           "directory"_a = ".");
};

void ASEOrcaPot::force(long nAtoms, const double *R, const int *atomicNrs,
                       double *F, double *U, double *variance,
                       const double *box) {
  variance = nullptr;
  // TODO(rg) Test after type maps
  MatrixType positions = MatrixType::Map(R, nAtoms, 3);
  MatrixType boxx = MatrixType::Map(box, 3, 3);
  Vector<int> atmnmrs = Vector<int>::Map(atomicNrs, nAtoms);
  py::object atoms = this->ase.attr("Atoms")(
      "symbols"_a = atmnmrs, "positions"_a = positions, "cell"_a = boxx);
  atoms.attr("set_calculator")(this->calc);
  atoms.attr("set_pbc")(std::tuple<bool, bool, bool>(true, true, true));
  double py_e = py::cast<double>(atoms.attr("get_potential_energy")());
  MatrixType py_force = py::cast<MatrixType>(atoms.attr("get_forces")());

  // Populate the output parameters
  *U = py_e;
  for (long i = 0; i < nAtoms; ++i) {
    F[3 * i] = py_force(i, 0);
    F[3 * i + 1] = py_force(i, 1);
    F[3 * i + 2] = py_force(i, 2);
  }
  counter++;
  return;
}
