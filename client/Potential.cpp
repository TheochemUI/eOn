#include <csignal>
#include <limits>
#include <time.h>
#include <utility>

#include "HelperFunctions.h"
#include "Log.h"
#include "Parameters.h"
#include "Potential.h"

#ifdef IMD_POT
#include "potentials/IMD/IMD.h"
#endif

#ifdef WITH_GPRD
#include "potentials/GPRPotential/GPRPotential.h"
#endif

#include "potentials/EAM/EAM.h"
#include "potentials/EMT/EffectiveMediumTheory.h"
#include "potentials/ExtPot/ExtPot.h"
#include "potentials/LJ/LJ.h"
#include "potentials/LJCluster/LJCluster.h"
#include "potentials/Morse/Morse.h"

#ifdef WITH_FORTRAN
#include "potentials/Aluminum/Aluminum.h"
#include "potentials/EDIP/EDIP.h"
#include "potentials/FeHe/FeHe.h"
#include "potentials/Lenosky/Lenosky.h"
#include "potentials/SW/SW.h"
#include "potentials/Tersoff/Tersoff.h"
#endif

#ifdef WITH_PYTHON
#ifdef PYAMFF_POT
#include "potentials/PyAMFF/PyAMFF.h"
#endif
#include "potentials/QSC/QSC.h"
#endif

#ifdef EONMPI
#include "potentials/MPIPot/MPIPot.h"
#endif

#ifdef LAMMPS_POT
#include "potentials/LAMMPS/LAMMPS.h"
#endif

#ifdef NEW_POT
#include "potentials/NewPot/NewPot.h"
#endif

// TODO: This should be guarded by WITH_FORTRAN as well
#ifdef CUH2_POT
#include "potentials/CuH2/CuH2.h"
#endif

#ifndef WIN32
#ifdef WITH_VASP
#include "potentials/VASP/VASP.h"
#endif
#endif

#ifdef WITH_AMS
#include "potentials/AMS/AMS.h"
#include "potentials/AMS_IO/AMS_IO.h"
#endif

#ifdef WITH_WATER
#include "potentials/Water/Water.hpp"
#ifdef WITH_FORTRAN
#include "potentials/Water_H/Tip4p_H.h"
#endif
#include "potentials/Water_Pt/Tip4p_Pt.hpp"
#endif

#include <limits>

namespace helper_functions {
std::unique_ptr<Potential> makePotential(Parameters *params) {
  switch(params->potential) {
  case PotType::EMT: {
    return (std::make_unique<EffectiveMediumTheory>(params));
    break;
  }
  case PotType::EXT: {
    return (std::make_unique<ExtPot>(params));
    break;
  }
  case PotType::LJ: {
    return (std::make_unique<LJ>(params));
    break;
  }
  case PotType::LJCLUSTER: {
    return (std::make_unique<LJCluster>(params));
    break;
  }
  case PotType::MORSE_PT: {
    return (std::make_unique<Morse>(params));
    break;
  }
#ifdef NEW_POT
  case PotType::NEW: {
    return (std::make_unique<NewPot>(params));
    break;
  }
#endif
#ifdef CUH2_POT
  case PotType::CUH2: {
    return (std::make_unique<CuH2>(params));
    break;
  }
#endif
#ifdef IMD_POT
  case PotType::IMD: {
    return (std::make_unique<IMD>(params));
    break;
  }
#endif
#ifdef WITH_WATER
  case PotType::TIP4P: {
    return (std::make_unique<Tip4p>(params));
    break;
  }
  case PotType::SPCE: {
    return (std::make_unique<SpceCcl>(params));
    break;
  }
#ifdef WITH_FORTRAN
  case PotType::TIP4P_PT: {
    return (std::make_unique<Tip4p_Pt>(params));
    break;
  }
  case PotType::TIP4P_H: {
    return (std::make_unique<Tip4p_H>(params));
    break;
  }
#endif
#endif
#ifdef WITH_FORTRAN
  case PotType::EAM_AL: {
    return (std::make_unique<Aluminum>(params));
    break;
  }
  case PotType::EDIP: {
    return (std::make_unique<EDIP>(params));
    break;
  }
  case PotType::FEHE: {
    return (std::make_unique<FeHe>(params));
    break;
  }
  case PotType::LENOSKY_SI: {
    return (std::make_unique<Lenosky>(params));
    break;
  }
  case PotType::SW_SI: {
    return (std::make_unique<SW>(params));
    break;
  }
  case PotType::TERSOFF_SI: {
    return (std::make_unique<Tersoff>(params));
    break;
  }
#endif
#ifndef WIN32
#ifdef WITH_VASP
  case PotType::VASP: {
    return (std::make_unique<VASP());
    break;
  }
#endif
#endif
#ifdef LAMMPS_POT
  case PotType::LAMMPS: {
    return (std::make_unique<lammps>(params));
    break;
  }
#endif
#ifdef EONMPI
  case PotType::MPI: {
    return (std::make_unique<MPIPot>(params));
    break;
  }
#endif
#ifdef WITH_PYTHON
#ifdef PYAMFF_POT
  case PotType::PYAMFF: {
    return (std::make_unique<PyAMFF());
    break;
  }
#endif
  case PotType::QSC: {
    return (std::make_unique<QSC());
    break;
  }
#endif
  // Unused
  // case PotType::BOPFOX: {
  //   return "bopfox"s;
  //   break;
  // }
  // case PotType::BOP: {
  //   return "bop"s;
  //   break;
  // }
#ifdef WITH_AMS
  case PotType::AMS: {
    return (std::make_unique<AMS>(params));
    break;
  }
  case PotType::AMS_IO: {
    return (std::make_unique<AMS_IO>(params));
    break;
  }
#endif
#ifdef WITH_GPRD
  case PotType::GPR: {
    return "gpr"s;
    break;
  }
#endif
  // case PotType::PYTHON: {
  //   TODO: Implement
  //   return "python"s;
  //   break;
  // }
  default:
    throw std::runtime_error("No known potential could be constructed");
    break;
  }
}

} // namespace helper_functions
