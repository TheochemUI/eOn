//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ASE_ORCA.h"
#include "Eigen/src/Core/Matrix.h"

ASEOrcaPot::ASEOrcaPot(shared_ptr<Parameters> a_params)
    : Potential(PotType::ASE_ORCA, a_params) {
  counter = 1;
  py::module_ sys = py::module_::import("sys");
  ase = py::module_::import("ase");
  py::module_ ase_orca = py::module_::import("ase.calculators.orca");
  py::module_ psutil = py::module_::import("psutil");

  // Set up ORCA profile and calculator
  py::object OrcaProfile = ase_orca.attr("OrcaProfile");
  py::object ORCA = ase_orca.attr("ORCA");

  std::string orca_simpleinput(fmt::format("ENGRAD {} {} {} {} UHF NOSOSCF", a_params->orca_pot, a_params->orca_basis, a_params->orca_grid, a_params->orca_extra_sline));
  this->calc = ORCA("profile"_a = OrcaProfile(py::str(a_params->orca_path)),
                    "orcasimpleinput"_a = orca_simpleinput,
                    "orcablocks"_a = py::str(fmt::format(
                        "%pal nprocs {} end",
                        py::cast<int>(psutil.attr("cpu_count")(false)))),
                    "directory"_a = ".");
};

void ASEOrcaPot::force(long nAtoms, const double *R, const int *atomicNrs,
                       double *F, double *U, double *variance,
                       const double *box) {
  variance = nullptr;
  Eigen::MatrixXd positions =
      Eigen::Map<Eigen::MatrixXd>(const_cast<double *>(R), nAtoms, 3);
  Eigen::MatrixXd boxx =
      Eigen::Map<Eigen::MatrixXd>(const_cast<double *>(box), 3, 3);
  Eigen::VectorXi atmnmrs =
      Eigen::Map<Eigen::VectorXi>(const_cast<int *>(atomicNrs), nAtoms);
  py::object atoms = this->ase.attr("Atoms")(
      "symbols"_a = atmnmrs, "positions"_a = positions, "cell"_a = boxx);
  atoms.attr("set_calculator")(this->calc);
  atoms.attr("set_pbc")(std::tuple<bool, bool, bool>(true, true, true));
  double py_e = py::cast<double>(atoms.attr("get_potential_energy")());
  Eigen::MatrixXd py_force =
      py::cast<Eigen::MatrixXd>(atoms.attr("get_forces")());

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
