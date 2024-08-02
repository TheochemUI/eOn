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
#include <csignal>

#include "BaseStructures.h"
#include "Parser.hpp"
#include "Potential.h"
#include "potentials/ParseTOML.hpp"

// #ifdef WITH_CATLEARN
// #include "potentials/CatLearnPot/CatLearnPot.h"
// #endif

// #ifdef IMD_POT
// #include "potentials/IMD/IMD.h"
// #endif

// #ifdef WITH_GPRD
// #include "potentials/GPRPotential/GPRPotential.h"
// #endif

// #include "potentials/EAM/EAM.h"
// #include "potentials/EMT/EffectiveMediumTheory.h"
// #include "potentials/ExtPot/ExtPot.h"
#include "potentials/LJ/LJ.h"
// #include "potentials/LJCluster/LJCluster.h"
#include "potentials/Morse/Morse.h"

#ifdef WITH_FORTRAN
#include "potentials/Aluminum/Aluminum.h"
#include "potentials/CuH2/CuH2.h"
#include "potentials/EDIP/EDIP.h"
#include "potentials/FeHe/FeHe.h"
#include "potentials/Lenosky/Lenosky.h"
#include "potentials/SW/SW.h"
#include "potentials/Tersoff/Tersoff.h"
#endif

// #ifdef WITH_PYTHON

// #ifdef PYAMFF_POT
// #include "potentials/PyAMFF/PyAMFF.h"
// #endif
// #ifdef ASE_POT
// #include "potentials/ASE/ASE.h"
// #endif

// #include "potentials/QSC/QSC.h"
// #endif

// #ifdef EONMPI
// #include "potentials/MPIPot/MPIPot.h"
// #endif

// #ifdef LAMMPS_POT
// #include "potentials/LAMMPS/LAMMPSPot.h"
// #endif

// #ifdef NEW_POT
// #include "potentials/NewPot/NewPot.h"
// #endif

// #ifndef WIN32
// #ifdef WITH_VASP
// #include "potentials/VASP/VASP.h"
// #endif
// #endif

// #ifdef WITH_AMS
// #include "potentials/AMS/AMS.h"
// #include "potentials/AMS_IO/AMS_IO.h"
// #endif

#ifdef WITH_ASE_ORCA
#include "potentials/ASE_ORCA/ASE_ORCA.h"
#endif

// #ifdef WITH_WATER
// #include "potentials/Water/Water.hpp"
// #ifdef WITH_FORTRAN
// #include "potentials/Water_H/Tip4p_H.h"
// #endif
// #include "potentials/Water_Pt/Tip4p_Pt.hpp"
// #endif

// // Should respect Fortran availability

// #ifdef WITH_XTB
// #include "potentials/XTBPot/XTBPot.h"
// #endif

namespace eonc {

std::shared_ptr<PotBase> makePotential(const toml::table &config) {
  config_section(config, "Potential");
  auto ptype = get_enum_toml<PotType>(config["Potential"]["potential"]);
  switch (ptype) {
  // case PotType::EMT: {
  //   return (std::make_shared<EffectiveMediumTheory>(a_p.pot.EMTRasmussen,
  //                                                   a_p.main.usePBC));
  //   break;
  // }
  // case PotType::EXT: {
  //   return (std::make_shared<ExtPot>(a_p.pot.extPotPath));
  //   break;
  // }
  case PotType::LJ: {
    auto params = LJ::Params();
    eonc::pot::from_toml(params, config["Potential"]["LJ"]);
    return (std::make_shared<LJ>(params));
    break;
  }
    //   case PotType::LJCLUSTER: {
    //     return (std::make_shared<LJCluster>(a_p.pot.lj));
    //     break;
    //   }
  case PotType::MORSE_PT: {
    auto params = Morse::Params();
    eonc::pot::from_toml(params, config["Potential"]["Morse"]);
    return (std::make_shared<Morse>(params));
    break;
  }
// #ifdef NEW_POT
//   case PotType::NEW: {
//     return (std::make_shared<NewPot>(a_p));
//     break;
//   }
// #endif
// #ifdef IMD_POT
//   case PotType::IMD: {
//     return (std::make_shared<IMD>());
//     break;
//   }
// #endif
// #ifdef WITH_WATER
//   case PotType::TIP4P: {
//     return (std::make_shared<Tip4p>());
//     break;
//   }
//   case PotType::SPCE: {
//     return (std::make_shared<SpceCcl>(a_p));
//     break;
//   }
// #ifdef WITH_FORTRAN
//   case PotType::TIP4P_PT: {
//     return (std::make_shared<Tip4p_Pt>(a_p));
//     break;
//   }
//   case PotType::TIP4P_H: {
//     return (std::make_shared<Tip4p_H>(a_p));
//     break;
//   }
// #endif
// #endif
#ifdef WITH_FORTRAN
  case PotType::EAM_AL: {
    return (std::make_shared<Aluminum>());
    break;
  }
  case PotType::CUH2: {
    return (std::make_shared<CuH2>());
    break;
  }
  case PotType::EDIP: {
    return (std::make_shared<EDIP>());
    break;
  }
  case PotType::FEHE: {
    return (std::make_shared<FeHe>());
    break;
  }
  case PotType::LENOSKY_SI: {
    return (std::make_shared<Lenosky>());
    break;
  }
  case PotType::SW_SI: {
    return (std::make_shared<SW>());
    break;
  }
  case PotType::TERSOFF_SI: {
    return (std::make_shared<Tersoff>());
    break;
  }
#endif
// #ifndef WIN32
// #ifdef WITH_VASP
//   case PotType::VASP: {
//     return (std::make_shared<VASP>(a_p));
//     break;
//   }
// #endif
// #endif
// #ifdef LAMMPS_POT
//   case PotType::LAMMPS: {
//     return (std::make_shared<LAMMPSPot>(params));
//     break;
//   }
// #endif
// #ifdef EONMPI
//   case PotType::MPI: {
//     return (
//         std::make_shared<MPIPot>(a_p.MPIPotentialRank,
//         a_p.pot.MPIPollPeriod));
//     break;
//   }
// #endif
// #ifdef WITH_PYTHON
// #ifdef PYAMFF_POT
//   case PotType::PYAMFF: {
//     return (std::make_shared<PyAMFF>());
//     break;
//   }
// #endif
// #ifdef ASE_POT
//   case PotType::ASE_POT: {
//     return (std::make_shared<ASE_POT>(a_p.pot.extPotPath));
//     break;
//   }
// #endif
//   case PotType::QSC: {
//     return (std::make_shared<QSC>());
//     break;
//   }
// #endif
//   // Unused
//   // case PotType::BOPFOX: {
//   //   return "bopfox"s;
//   //   break;
//   // }
//   // case PotType::BOP: {
//   //   return "bop"s;
//   //   break;
//   // }
// #ifdef WITH_AMS
//   case PotType::AMS: {
//     return (std::make_shared<AMS>(a_p.ams));
//     break;
//   }
//   case PotType::AMS_IO: {
//     return (std::make_shared<AMS_IO>(a_p.ams));
//     break;
//   }
// #endif
// #ifdef WITH_GPRD
//   // case PotType::GPR: {
//   //   return "gpr"s;
//   //   break;
//   // }
// #endif
//   // case PotType::PYTHON: {
//   //   TODO: Implement
//   //   return "python"s;
//   //   break;
//   // }
// #ifdef WITH_CATLEARN
//   case PotType::CatLearn: {
//     return (std::make_shared<CatLearnPot>(a_p.catl));
//     break;
//   }
// #endif
// // TODO: Handle Fortran interaction
// #ifdef WITH_XTB
//   case PotType::XTB: {
//     return (std::make_shared<XTBPot>(a_p));
//     break;
//   }
// #endif
#ifdef WITH_ASE_ORCA
  case PotType::ASE_ORCA: {
    auto params = ASEOrcaPot::Params();
    eonc::pot::from_toml(params, config["Potential"]["ASE_ORCA"]);
    return (std::make_shared<ASEOrcaPot>(params));
    break;
  }
#endif
  default:
    SPDLOG_ERROR("No known potential could be constructed from {}",
                 magic_enum::enum_name(ptype));
    throw std::runtime_error("Terminating");
    break;
  }
}

} // namespace eonc
