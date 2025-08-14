#include "ASE_NWCHEM.h"
#include "../../EnvHelpers.hpp"
#include "../../fpe_handler.h"
#include <string>

// TODO(rg): Clean this up.
ASENwchemPot::ASENwchemPot(std::shared_ptr<Parameters> a_params)
    : Potential(PotType::ASE_NWCHEM, a_params) {
  counter = 0;
  py::module_ sys = py::module_::import("sys");
  // Fix for gh-184, see
  // https://github.com/numpy/numpy/issues/20504#issuecomment-985542508
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();
  ase = py::module_::import("ase");
  fpeh.restore_fpe();
  py::module_ ase_nwchem = py::module_::import("ase.calculators.nwchem");
  py::module_ psutil = py::module_::import("psutil");
  std::string nwchempth = helper_functions::get_value_from_env_or_param(
      "NWCHEM_COMMAND", a_params->nwchem_path, "", "", true);
  std::string nwc_mult = helper_functions::get_value_from_env_or_param(
      "NWCHEM_MULTIPLICITY", a_params->nwchem_multiplicity, "1",
      "Using 1 as a default multiplicity, i.e. an RHF calculation suitable for "
      "closed shell molecules, set multiplicity or the "
      "environment variable NWCHEM_MULTIPLICITY.\n");

  // Set up NWCHEM arguments
  // TODO(rg): Stop hardcoding these
  py::object NWCHEM = ase_nwchem.attr("NWChem");
  size_t nproc{0};
  auto mult = std::stoi(nwc_mult); // 1 for singlet, 2 for doublet

  // TODO(rg): Use
  if (a_params->nwchem_nproc == "auto") {
    nproc = py::cast<int>(psutil.attr("cpu_count")(false));
  } else {
    nproc = std::stoi(a_params->nwchem_nproc);
  }

  // TODO(rg): Directory should be set by the user, and created here
  // dont_verify so we always get an energy and gradient
  // Common NWCHEM parameters
  py::dict nwchem_params = py::dict(
      "label"_a = "_eonpot_engrad",
      "set"_a = py::dict("geom:dont_verify"_a = true),
      "command"_a = py::str(fmt::format(
          "mpirun -n {} {} PREFIX.nwi > PREFIX.nwo", nproc, nwchempth)),
      "memory"_a = py::str("2 gb"),
      "scf"_a = py::dict("nopen"_a = mult - 1,
                         "thresh"_a = a_params->nwchem_scf_thresh,
                         "maxiter"_a = a_params->nwchem_scf_maxiter),
      "basis"_a = py::str("3-21G"), "task"_a = py::str("gradient"),
      "directory"_a = ".");

  // Set flag for doublet (mult == 2)
  if (mult == 2) {
    nwchem_params["scf"]["uhf"] = py::none();
  }

  // Check for invalid spin multiplicity
  if (mult != 1 && mult != 2) {
    throw std::runtime_error("Unknown spin multiplicity, we support 1 for "
                             "singlet and 2 for doublet ONLY.");
  }

  this->calc = NWCHEM(**nwchem_params);
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
  atoms.attr("calc") = this->calc;
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
