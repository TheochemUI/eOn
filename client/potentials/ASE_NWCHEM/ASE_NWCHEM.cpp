#include "eon/potentials/ASE_NWCHEM/ASE_NWCHEM.h"
#include "eon/EnvHelpers.hpp"
#include "eon/EonLogger.h"
#include "eon/PyGuard.h"
#include "eon/fpe_handler.h"

#include <atomic>
#include <format>
#include <stdexcept>
#include <string>
#include <system_error>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

namespace {

std::filesystem::path makeAseWorkDir(const char *prefix) {
  static std::atomic<long> instanceCount{0};
#ifdef _WIN32
  const long pid = static_cast<long>(_getpid());
#else
  const long pid = static_cast<long>(getpid());
#endif
  std::error_code ec;
  auto dir = std::filesystem::temp_directory_path(ec);
  if (ec) {
    throw std::runtime_error(std::format(
        "ASE-NWChem: cannot resolve temp directory: {}", ec.message()));
  }
  dir /= std::format("{}_{}_{}", prefix, pid, instanceCount.fetch_add(1));
  std::filesystem::create_directories(dir, ec);
  if (ec) {
    throw std::runtime_error(
        std::format("ASE-NWChem: cannot create work directory {}: {}",
                    dir.string(), ec.message()));
  }
  return dir;
}

} // namespace

// TODO(rg): Clean this up.
ASENwchemPot::ASENwchemPot(const Parameters &a_params)
    : Potential(PotType::ASE_NWCHEM, a_params) {
  eonc::ensure_interpreter();
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
  std::string nwchempth = eonc::helpers::get_value_from_env_or_param(
      "NWCHEM_COMMAND", a_params.ase_nwchem_options.path, "", "", true);
  std::string nwc_mult = eonc::helpers::get_value_from_env_or_param(
      "NWCHEM_MULTIPLICITY", a_params.ase_nwchem_options.multiplicity, "1",
      "Using 1 as a default multiplicity, i.e. an RHF calculation suitable for "
      "closed shell molecules, set multiplicity or the "
      "environment variable NWCHEM_MULTIPLICITY.\n");

  // Set up NWCHEM arguments
  // TODO(rg): Stop hardcoding these
  py::object NWCHEM = ase_nwchem.attr("NWChem");
  size_t nproc{0};
  auto mult = std::stoi(nwc_mult); // 1 for singlet, 2 for doublet

  // TODO(rg): Use
  if (a_params.ase_nwchem_options.nproc == "auto") {
    nproc = py::cast<int>(psutil.attr("cpu_count")(false));
  } else {
    nproc = std::stoi(a_params.ase_nwchem_options.nproc);
  }

  // dont_verify so we always get an energy and gradient
  // mpi_launcher: mpirun (default) or srun on Slurm nodes (issue #193)
  const std::string &launcher = a_params.ase_nwchem_options.mpi_launcher;
  std::string mpi_cmd;
  if (launcher == "srun") {
    // srun uses -n for tasks; avoid OpenMPI-specific flags.
    mpi_cmd =
        std::format("srun -n {} {} PREFIX.nwi > PREFIX.nwo", nproc, nwchempth);
  } else {
    // Default and any other launcher treated as OpenMPI-style mpirun -n.
    mpi_cmd = std::format("{} -n {} {} PREFIX.nwi > PREFIX.nwo", launcher,
                          nproc, nwchempth);
  }

  if (mult != 1 && mult != 2) {
    throw std::runtime_error("Unknown spin multiplicity, we support 1 for "
                             "singlet and 2 for doublet ONLY.");
  }

  // One directory per calculator instance so two LocalInProcess jobs in
  // the same cwd do not write PREFIX.nwi / PREFIX.nwo over each other.
  workDir = makeAseWorkDir("eon_ase_nwchem");

  // Common NWCHEM parameters
  py::dict nwchem_params = py::dict(
      "label"_a = "_eonpot_engrad",
      "set"_a = py::dict("geom:dont_verify"_a = true),
      "command"_a = py::str(mpi_cmd), "memory"_a = py::str("2 gb"),
      "scf"_a = py::dict("nopen"_a = mult - 1,
                         "thresh"_a = a_params.ase_nwchem_options.scf_thresh,
                         "maxiter"_a = a_params.ase_nwchem_options.scf_maxiter),
      "basis"_a = py::str("3-21G"), "task"_a = py::str("gradient"),
      "directory"_a = workDir.string());

  // Set flag for doublet (mult == 2)
  if (mult == 2) {
    nwchem_params["scf"]["uhf"] = py::none();
  }

  try {
    this->calc = NWCHEM(**nwchem_params);
  } catch (...) {
    std::error_code ec;
    std::filesystem::remove_all(workDir, ec);
    workDir.clear();
    throw;
  }
};

ASENwchemPot::~ASENwchemPot() {
  QUILL_LOG_INFO(eonc::log::get(), "[ASENwchem] called potential {} times",
                 counter);
  if (!workDir.empty()) {
    std::error_code ec;
    std::filesystem::remove_all(workDir, ec);
  }
}

void ASENwchemPot::force(long nAtoms, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  variance = nullptr;
  try {
    AtomMatrix positions = AtomMatrix::Map(const_cast<double *>(R), nAtoms, 3);
    RotationMatrix boxx = RotationMatrix::Map(const_cast<double *>(box), 3, 3);
    Eigen::VectorXi atmnmrs =
        Eigen::Map<Eigen::VectorXi>(const_cast<int *>(atomicNrs), nAtoms);
    // XXX: NWChem refuses to perform SCF for anything but a molecule, so no box
    // or pbc can be passed
    py::object atoms = this->ase.attr("Atoms")("symbols"_a = atmnmrs,
                                               "positions"_a = positions);
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
  } catch (py::error_already_set &e) {
    throw std::runtime_error(std::string("ASE-NWChem Python error: ") +
                             e.what());
  } catch (const std::exception &e) {
    throw std::runtime_error(std::string("ASE-NWChem C++ exception: ") +
                             e.what());
  }

  counter++;
  return;
}
