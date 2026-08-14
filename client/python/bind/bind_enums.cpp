#include "eon/BaseStructures.h"
#include "eon/ConFileIO.h"
#include "eon/Matter.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <magic_enum/magic_enum.hpp>
#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_enums(nb::module_ &m) {
  nb::enum_<eonc::io::IoStatus>(m, "IoStatus")
      .value("Ok", eonc::io::IoStatus::Ok)
      .value("ReadError", eonc::io::IoStatus::ReadError)
      .value("WriteError", eonc::io::IoStatus::WriteError)
      .value("AppendError", eonc::io::IoStatus::AppendError)
      .value("OpenError", eonc::io::IoStatus::OpenError)
      .value("InvalidArgument", eonc::io::IoStatus::InvalidArgument);

  m.def(
      "io_ok", [](eonc::io::IoStatus s) { return eonc::io::io_ok(s); },
      nb::arg("status"));
  m.def(
      "io_status_name",
      [](eonc::io::IoStatus s) {
        return std::string(eonc::io::io_status_name(s));
      },
      nb::arg("status"));

  nb::enum_<eonc::PbcConvention>(m, "PbcConvention")
      .value("Legacy", eonc::PbcConvention::Legacy)
      .value("MinimumImage", eonc::PbcConvention::MinimumImage);

  // Full PotType surface (matches config names / magic_enum)
  nb::enum_<eonc::PotType>(m, "PotType")
      .value("UNKNOWN", eonc::PotType::UNKNOWN)
      .value("EMT", eonc::PotType::EMT)
      .value("EXT_POT", eonc::PotType::EXT_POT)
      .value("LJ", eonc::PotType::LJ)
      .value("LJCLUSTER", eonc::PotType::LJCLUSTER)
      .value("MORSE_PT", eonc::PotType::MORSE_PT)
      .value("CUH2", eonc::PotType::CUH2)
      .value("TIP4P", eonc::PotType::TIP4P)
      .value("TIP4P_PT", eonc::PotType::TIP4P_PT)
      .value("TIP4P_H", eonc::PotType::TIP4P_H)
      .value("SPCE", eonc::PotType::SPCE)
      .value("EAM_AL", eonc::PotType::EAM_AL)
      .value("EDIP", eonc::PotType::EDIP)
      .value("FEHE", eonc::PotType::FEHE)
      .value("LENOSKY_SI", eonc::PotType::LENOSKY_SI)
      .value("SW_SI", eonc::PotType::SW_SI)
      .value("TERSOFF_SI", eonc::PotType::TERSOFF_SI)
      .value("VASP", eonc::PotType::VASP)
      .value("LAMMPS", eonc::PotType::LAMMPS)
      .value("MPI", eonc::PotType::MPI)
      .value("AMS", eonc::PotType::AMS)
      .value("AMS_IO", eonc::PotType::AMS_IO)
      .value("GPR", eonc::PotType::GPR)
      .value("CatLearn", eonc::PotType::CatLearn)
      .value("XTB", eonc::PotType::XTB)
      .value("ASE_ORCA", eonc::PotType::ASE_ORCA)
      .value("ASE_POT", eonc::PotType::ASE_POT)
      .value("ASE_NWCHEM", eonc::PotType::ASE_NWCHEM)
      .value("METATOMIC", eonc::PotType::METATOMIC)
      .value("ZBL", eonc::PotType::ZBL)
      .value("SocketNWChem", eonc::PotType::SocketNWChem)
      .value("RGPOT", eonc::PotType::RGPOT);

  m.def(
      "pot_type_from_name",
      [](const std::string &name) {
        auto v = magic_enum::enum_cast<eonc::PotType>(
            name, magic_enum::case_insensitive);
        if (!v)
          throw std::invalid_argument("unknown PotType: " + name);
        return *v;
      },
      nb::arg("name"));
  m.def(
      "pot_type_name",
      [](eonc::PotType t) { return std::string(magic_enum::enum_name(t)); },
      nb::arg("pot_type"));

  nb::enum_<eonc::JobType>(m, "JobType")
      .value("Unknown", eonc::JobType::Unknown)
      .value("Process_Search", eonc::JobType::Process_Search)
      .value("Saddle_Search", eonc::JobType::Saddle_Search)
      .value("Minimization", eonc::JobType::Minimization)
      .value("Point", eonc::JobType::Point)
      .value("Parallel_Replica", eonc::JobType::Parallel_Replica)
      .value("Safe_Hyperdynamics", eonc::JobType::Safe_Hyperdynamics)
      .value("TAD", eonc::JobType::TAD)
      .value("Replica_Exchange", eonc::JobType::Replica_Exchange)
      .value("Basin_Hopping", eonc::JobType::Basin_Hopping)
      .value("Hessian", eonc::JobType::Hessian)
      .value("Finite_Difference", eonc::JobType::Finite_Difference)
      .value("Nudged_Elastic_Band", eonc::JobType::Nudged_Elastic_Band)
      .value("Dynamics", eonc::JobType::Dynamics)
      .value("Prefactor", eonc::JobType::Prefactor)
      .value("Global_Optimization", eonc::JobType::Global_Optimization)
      .value("Structure_Comparison", eonc::JobType::Structure_Comparison)
      .value("Monte_Carlo", eonc::JobType::Monte_Carlo)
      .value("Test", eonc::JobType::Test)
      .value("GP_Surrogate", eonc::JobType::GP_Surrogate);

  m.def(
      "job_type_from_name",
      [](const std::string &name) {
        auto v = magic_enum::enum_cast<eonc::JobType>(
            name, magic_enum::case_insensitive);
        if (!v)
          throw std::invalid_argument("unknown JobType: " + name);
        return *v;
      },
      nb::arg("name"));
  m.def(
      "job_type_name",
      [](eonc::JobType t) { return std::string(magic_enum::enum_name(t)); },
      nb::arg("job_type"));

  nb::enum_<eonc::OptType>(m, "OptType")
      .value("Unknown", eonc::OptType::Unknown)
      .value("None_", eonc::OptType::None) // None is reserved in Python
      .value("QM", eonc::OptType::QM)
      .value("CG", eonc::OptType::CG)
      .value("LBFGS", eonc::OptType::LBFGS)
      .value("FIRE", eonc::OptType::FIRE)
      .value("SD", eonc::OptType::SD);

  nb::enum_<eonc::NEBInit>(m, "NEBInit")
      .value("LINEAR", eonc::NEBInit::LINEAR)
      .value("IDPP", eonc::NEBInit::IDPP)
      .value("IDPP_COLLECTIVE", eonc::NEBInit::IDPP_COLLECTIVE)
      .value("SIDPP", eonc::NEBInit::SIDPP)
      .value("SIDPP_ZBL", eonc::NEBInit::SIDPP_ZBL)
      .value("FILE", eonc::NEBInit::FILE);

  nb::enum_<eonc::RunStatus>(m, "RunStatus")
      .value("GOOD", eonc::RunStatus::GOOD)
      .value("FAIL_MAX_ITERATIONS", eonc::RunStatus::FAIL_MAX_ITERATIONS)
      .value("FAIL_POTENTIAL_FAILED", eonc::RunStatus::FAIL_POTENTIAL_FAILED);
}

} // namespace eonc::pybind
