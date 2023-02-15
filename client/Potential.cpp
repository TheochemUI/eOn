#include <time.h>
#include <utility>

#include "HelperFunctions.h"
#include "Log.h"
#include "Parameters.h"
#include "Potential.h"

#include "potentials/EMT/EffectiveMediumTheory.h"
#include "potentials/LJ/LJ.h"
#include "potentials/Morse/Morse.h"

std::pair<double, AtomMatrix> Potential::get_ef(long nAtoms,
                                                const double *positions,
                                                const int *atomicNrs,
                                                const double *box) {
  double energy{0};
  AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
  this->force(nAtoms, positions, atomicNrs, forces.data(), &energy, box);
  return std::make_pair(energy, forces);
}

namespace helper_functions {
Potential *makePotential(Parameters *params) {
  switch (params->potential) {
  case PotType::EMT: {
    return (new EffectiveMediumTheory(params));
    break;
  }
  case PotType::LJ: {
    return (new LJ(params));
    break;
  }
  case PotType::MORSE_PT: {
    return (new Morse(params));
    break;
  }
    // case PotType::EXT: {
    //   return "ext"s;
    //   break;
    // }
    // case PotType::LJ: {
    //   return "lj"s;
    //   break;
    // }
    // case PotType::LJCLUSTER: {
    //   return "ljcluster"s;
    //   break;
    // }
    // case PotType::MORSE_PT: {
    //   return "morse_pt"s;
    //   break;
    // }
    // case PotType::NEW: {
    //   return "new"s;
    //   break;
    // }
    // case PotType::CUH2: {
    //   return "cuh2"s;
    //   break;
    // }
    // case PotType::IMD: {
    //   return "imd"s;
    //   break;
    // }
    // case PotType::TIP4P: {
    //   return "tip4p"s;
    //   break;
    // }
    // case PotType::TIP4P_PT: {
    //   return "tip4p_pt"s;
    //   break;
    // }
    // case PotType::TIP4P_H: {
    //   return "tip4p_h"s;
    //   break;
    // }
    // case PotType::SPCE: {
    //   return "spce"s;
    //   break;
    // }
    // case PotType::EAM_AL: {
    //   return "eam_al"s;
    //   break;
    // }
    // case PotType::EDIP: {
    //   return "edip"s;
    //   break;
    // }
    // case PotType::FEHE: {
    //   return "fehe"s;
    //   break;
    // }
    // case PotType::LENOSKY_SI: {
    //   return "lenosky_si"s;
    //   break;
    // }
    // case PotType::SW_SI: {
    //   return "sw_si"s;
    //   break;
    // }
    // case PotType::TERSOFF_SI: {
    //   return "tersoff_si"s;
    //   break;
    // }
    // case PotType::VASP: {
    //   return "vasp"s;
    //   break;
    // }
    // case PotType::LAMMPS: {
    //   return "lammps"s;
    //   break;
    // }
    // case PotType::MPI: {
    //   return "mpi"s;
    //   break;
    // }
    // case PotType::PYAMFF: {
    //   return "pyamff"s;
    //   break;
    // }
    // case PotType::QSC: {
    //   return "qsc"s;
    //   break;
    // }
    // case PotType::BOPFOX: {
    //   return "bopfox"s;
    //   break;
    // }
    // case PotType::BOP: {
    //   return "bop"s;
    //   break;
    // }
    // case PotType::AMS: {
    //   return "ams"s;
    //   break;
    // }
    // case PotType::AMS_IO: {
    //   return "ams_io"s;
    //   break;
    // }
    // case PotType::GPR: {
    //   return "gpr"s;
    //   break;
    // }
    // case PotType::PYTHON: {
    //   return "python"s;
    //   break;
    // }
    // default:
    //   return "unknown potential"s;
    //   break;
  }
}
} // namespace helper_functions
