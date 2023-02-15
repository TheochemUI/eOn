#include "BaseStructures.h"
#include <cassert>
#include <iostream>

namespace helper_functions {
std::string getPotentialName(PotType ptype) {
  switch (ptype) {
  case PotType::EMT: {
    return "emt"s;
    break;
  }
  case PotType::EXT: {
    return "ext"s;
    break;
  }
  case PotType::LJ: {
    return "lj"s;
    break;
  }
  case PotType::LJCLUSTER: {
    return "ljcluster"s;
    break;
  }
  case PotType::MORSE_PT: {
    return "morse_pt"s;
    break;
  }
  case PotType::NEW: {
    return "new"s;
    break;
  }
  case PotType::CUH2: {
    return "cuh2"s;
    break;
  }
  case PotType::IMD: {
    return "imd"s;
    break;
  }
  case PotType::TIP4P: {
    return "tip4p"s;
    break;
  }
  case PotType::TIP4P_PT: {
    return "tip4p_pt"s;
    break;
  }
  case PotType::TIP4P_H: {
    return "tip4p_h"s;
    break;
  }
  case PotType::SPCE: {
    return "spce"s;
    break;
  }
  case PotType::EAM_AL: {
    return "eam_al"s;
    break;
  }
  case PotType::EDIP: {
    return "edip"s;
    break;
  }
  case PotType::FEHE: {
    return "fehe"s;
    break;
  }
  case PotType::LENOSKY_SI: {
    return "lenosky_si"s;
    break;
  }
  case PotType::SW_SI: {
    return "sw_si"s;
    break;
  }
  case PotType::TERSOFF_SI: {
    return "tersoff_si"s;
    break;
  }
  case PotType::VASP: {
    return "vasp"s;
    break;
  }
  case PotType::LAMMPS: {
    return "lammps"s;
    break;
  }
  case PotType::MPI: {
    return "mpi"s;
    break;
  }
  case PotType::PYAMFF: {
    return "pyamff"s;
    break;
  }
  case PotType::QSC: {
    return "qsc"s;
    break;
  }
  case PotType::BOPFOX: {
    return "bopfox"s;
    break;
  }
  case PotType::BOP: {
    return "bop"s;
    break;
  }
  case PotType::AMS: {
    return "ams"s;
    break;
  }
  case PotType::AMS_IO: {
    return "ams_io"s;
    break;
  }
  case PotType::GPR: {
    return "gpr"s;
    break;
  }
  case PotType::PYTHON: {
    return "python"s;
    break;
  }
  default:
    return "unknown potential"s;
    break;
  }
}
PotType getPotentialType(const std::string pname) {
  if (pname == "emt"s) {
    return PotType::EMT;
  }
  else if (pname == "ext"s) {
    return PotType::EXT;
  }
  else if (pname == "lj"s) {
    return PotType::LJ;
  }
  else if (pname == "ljcluster"s) {
    return PotType::LJCLUSTER;
  }
  else if (pname == "morse_pt"s) {
    return PotType::MORSE_PT;
  }
  else if (pname == "new"s) {
    return PotType::NEW;
  }
  else if (pname == "cuh2"s) {
    return PotType::CUH2;
  }
  else if (pname == "imd"s) {
    return PotType::IMD;
  }
  else if (pname == "tip4p"s) {
    return PotType::TIP4P;
  }
  else if (pname == "tip4p_pt"s) {
    return PotType::TIP4P_PT;
  }
  else if (pname == "tip4p_h"s) {
    return PotType::TIP4P_H;
  }
  else if (pname == "spce"s) {
    return PotType::SPCE;
  }
  else if (pname == "eam_al"s) {
    return PotType::EAM_AL;
  }
  else if (pname == "edip"s) {
    return PotType::EDIP;
  }
  else if (pname == "fehe"s) {
    return PotType::FEHE;
  }
  else if (pname == "lenosky_si"s) {
    return PotType::LENOSKY_SI;
  }
  else if (pname == "sw_si"s) {
    return PotType::SW_SI;
  }
  else if (pname == "tersoff_si"s) {
    return PotType::TERSOFF_SI;
  }
  else if (pname == "vasp"s) {
    return PotType::VASP;
  }
  else if (pname == "lammps"s) {
    return PotType::LAMMPS;
  }
  else if (pname == "mpi"s) {
    return PotType::MPI;
  }
  else if (pname == "pyamff"s) {
    return PotType::PYAMFF;
  }
  else if (pname == "qsc"s) {
    return PotType::QSC;
  }
  else if (pname == "bopfox"s) {
    return PotType::BOPFOX;
  }
  else if (pname == "bop"s) {
    return PotType::BOP;
  }
  else if (pname == "python"s) {
    return PotType::PYTHON;
  }
  else {
    return PotType::UNKNOWN;
  }
}
} // namespace helper_functions
