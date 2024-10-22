//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ASE_NWCHEM.h"
#include "../../EnvHelpers.hpp"
#include <fenv.h>

// TODO(rg): Clean this up.
ASENwchemPot::ASENwchemPot(std::shared_ptr<Parameters> a_params)
    : Potential(PotType::ASE_NWCHEM, a_params) {
  counter = 0;
  py::module_ sys = py::module_::import("sys");
  // Fix for gh-184, see
  // https://github.com/numpy/numpy/issues/20504#issuecomment-985542508
  fenv_t orig_feenv;
  feholdexcept(&orig_feenv);
  ase = py::module_::import("ase");
  fesetenv(&orig_feenv);
  py::module_ ase_nwchem = py::module_::import("ase.calculators.nwchem");
  py::module_ psutil = py::module_::import("psutil");
  std::string nwchempth = helper_functions::get_value_from_env_or_param(
      "NWCHEM_COMMAND", a_params->nwchem_path, "", "", true);

  // Set up NWCHEM arguments
  // TODO(rg): Stop hardcoding these
  py::object NWCHEM = ase_nwchem.attr("NWChem");
  size_t nproc{0};
  size_t mult = 1; // 1 for singlet, 2 for doublet

  // TODO(rg): Use
  if (a_params->nwchem_nproc == "auto") {
    nproc = py::cast<int>(psutil.attr("cpu_count")(false));
  } else {
    nproc = std::stoi(a_params->nwchem_nproc);
  }

  // TODO(rg): Directory should be set by the user, and created here
  // dont_verify so we always get an energy and gradient
  this->calc =
      NWCHEM("label"_a = "_eonpot_engrad",
             "set"_a = py::dict("geom:dont_verify"_a = true),
             "command"_a = py::str(fmt::format(
                 "mpirun -n {} {} PREFIX.nwi > PREFIX.nwo", nproc, nwchempth)),
             "memory"_a = py::str("2 gb"),
             "scf"_a = py::dict("nopen"_a = mult - 1, "thresh"_a = 1e-6,
                                "maxiter"_a = 200),
             "basis"_a = py::str("3-21G"), "task"_a = py::str("gradient"),
             "directory"_a = ".");
};

void ASENwchemPot::force(long nAtoms, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  variance = nullptr;
  Eigen::MatrixXd positions =
      Eigen::Map<Eigen::MatrixXd>(const_cast<double *>(R), nAtoms, 3);
  Eigen::MatrixXd boxx =
      Eigen::Map<Eigen::MatrixXd>(const_cast<double *>(box), 3, 3);
  Eigen::VectorXi atmnmrs =
      Eigen::Map<Eigen::VectorXi>(const_cast<int *>(atomicNrs), nAtoms);
  // XXX: NWChem refuses to perform SCF for anything but a molecule, so no box
  // or pbc can be passed
  py::object atoms =
      this->ase.attr("Atoms")("symbols"_a = atmnmrs, "positions"_a = positions);
  atoms.attr("set_calculator")(this->calc);
  // atoms.attr("center")();
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
