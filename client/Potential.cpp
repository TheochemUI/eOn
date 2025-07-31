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

#include "client/BaseStructures.h"
#include "client/Parser.hpp"
#include "client/Potential.h"
#include "client/potentials/ParseTOML.hpp"

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
// #include "client/potentials/EMT/EffectiveMediumTheory.h"
#include "client/potentials/ExtPot/ExtPot.h"
#include "client/potentials/LJ/LJ.h"
// #include "potentials/LJCluster/LJCluster.h"
#include "client/potentials/Morse/Morse.h"

#ifdef WITH_FORTRAN
#include "client/potentials/Aluminum/Aluminum.h"
#include "client/potentials/CuH2/CuH2.h"
#include "client/potentials/EDIP/EDIP.h"
#include "client/potentials/FeHe/FeHe.h"
#include "client/potentials/Lenosky/Lenosky.h"
#include "client/potentials/SW/SW.h"
#include "client/potentials/Tersoff/Tersoff.h"
#endif

// #ifdef WITH_PYTHON

// #ifdef PYAMFF_POT
// #include "client/potentials/PyAMFF/PyAMFF.h"
// #endif
// #ifdef ASE_POT
// #include "client/potentials/ASE/ASE.h"
// #endif

// #include "client/potentials/QSC/QSC.h"
// #endif

// #ifdef EONMPI
// #include "client/potentials/MPIPot/MPIPot.h"
// #endif

#ifdef LAMMPS_POT
#include "client/potentials/LAMMPS/LAMMPSPot.h"
#endif

// #ifdef NEW_POT
// #include "client/potentials/NewPot/NewPot.h"
// #endif

#ifndef WIN32
#ifdef WITH_VASP
#include "client/potentials/VASP/VASP.h"
#endif
#endif

// #ifdef WITH_AMS
// #include "client/potentials/AMS/AMS.h"
// #include "client/potentials/AMS_IO/AMS_IO.h"
// #endif

#ifdef WITH_ASE_ORCA
#include "client/potentials/ASE_ORCA/ASE_ORCA.h"
#endif

#ifdef WITH_ASE_NWCHEM
#include "potentials/ASE_NWCHEM/ASE_NWCHEM.h"
#endif

#ifdef WITH_METATOMIC
#include "potentials/Metatomic/MetatomicPotential.h"
#endif

#ifdef WITH_WATER
#include "client/potentials/Water/Water.hpp"
#ifdef WITH_FORTRAN
#include "client/potentials/Water_H/Tip4p_H.h"
#endif
#include "client/potentials/Water_Pt/Tip4p_Pt.hpp"
#endif

// // Should respect Fortran availability

// #ifdef WITH_XTB
// #include "client/potentials/XTBPot/XTBPot.h"
// #endif

namespace eonc {

std::shared_ptr<PotBase> makePotential(const toml::table &config) {
  config_section(config, "Potential");
  auto ptype = get_enum_toml<PotType>(config["Potential"]["potential"]);
  switch (ptype) {
  // case PotType::EMT: {
  //   return (std::make_shared<EffectiveMediumTheory>(false));
  //   break;
  // }
  // case PotType::EMT_RAS: {
  //   return (std::make_shared<EffectiveMediumTheory>(true));
  //   break;
  // }
  case PotType::EXT: {
    return (std::make_shared<ExtPot>(
        config["Potential"]["ext_pot_path"].value_or("ext_pot")));
    break;
  }
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
#ifdef WITH_WATER
  case PotType::TIP4P: {
    // TODO(rg): SVN returns different results for the forces of
    // client/gtests/data/systems/cuh2_neb_test/init.con
    return (std::make_shared<Tip4p>());
    break;
  }
  case PotType::SPCE: {
    // TODO(rg): SVN returns different results for the forces of
    // client/gtests/data/systems/cuh2_neb_test/init.con
    // client/example/pos.con
    return (std::make_shared<SpceCcl>());
    break;
  }
#ifdef WITH_FORTRAN
  case PotType::TIP4P_PT: {
    return (std::make_shared<Tip4p_Pt>());
    break;
  }
  case PotType::TIP4P_H: {
    return (std::make_shared<Tip4p_H>());
    break;
  }
#endif
#endif
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
#ifndef WIN32
#ifdef WITH_VASP
  case PotType::VASP: {
    // TODO(rg):: This should probably take parameters
    return (std::make_shared<VASP>());
    break;
  }
#endif
#endif
#ifdef LAMMPS_POT
  case PotType::LAMMPS: {
    auto params = LAMMPSPot::Params();
    eonc::pot::from_toml(params, config["Potential"]["LAMMPS"]);
    return (std::make_shared<LAMMPSPot>(params));
    break;
  }
#endif
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
#ifdef WITH_ASE_NWCHEM
  case PotType::ASE_NWCHEM: {
    return (std::make_shared<ASENwchemPot>(params));
    break;
  }
#endif
#ifdef WITH_METATOMIC
  case PotType::METATOMIC: {
    return (std::make_shared<MetatomicPotential>(params));
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
