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
#include "eon/EonLogger.h"
#include <csignal>
#include <ctime>
#include <limits>
#include <utility>

#include "eon/HelperFunctions.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#ifdef WITH_CATLEARN
#include "eon/potentials/CatLearnPot/CatLearnPot.h"
#endif

#ifdef IMD_POT
#include "eon/potentials/IMD/IMD.h"
#endif

#ifdef WITH_GPRD
#include "eon/potentials/GPRPotential/GPRPotential.h"
#endif

#include "eon/potentials/EAM/EAM.h"
#include "eon/potentials/EMT/EffectiveMediumTheory.h"
#include "eon/potentials/ExtPot/ExtPot.h"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJClusterPot.hpp"
#include "rgpot/LennardJones/LJPot.hpp"
#include "rgpot/Morse/MorsePot.hpp"
#include "rgpot/ZBL/ZBLPot.hpp"
#ifndef IS_WINDOWS
#include "eon/potentials/SocketNWChem/SocketNWChemPot.h"
#ifdef WITH_RGPOT
#include "eon/potentials/Rgpot/RgpotPot.h"
#endif
#endif

// Fortran potentials: always compiled, loaded at runtime via dlopen
#include "eon/potentials/Aluminum/Aluminum.h"
#include "eon/potentials/EDIP/EDIP.h"
#include "eon/potentials/FeHe/FeHe.h"
#include "eon/potentials/FortranPotLoader.h"
#include "eon/potentials/Lenosky/Lenosky.h"
#include "eon/potentials/SW/SW.h"
#include "eon/potentials/Tersoff/Tersoff.h"

#ifdef EMBED_PYTHON

#ifdef PYAMFF_POT
#include "eon/potentials/PyAMFF/PyAMFF.h"
#endif
#ifdef WITH_ASE_POT
#include "eon/potentials/ASE/ASE.h"
#endif

#include "eon/potentials/QSC/QSC.h"
#endif

#ifdef EONMPI
#include "eon/potentials/MPIPot/MPIPot.h"
#endif

#include "eon/potentials/LAMMPS/LAMMPSPot.h"

#ifdef NEW_POT
#include "eon/potentials/NewPot/NewPot.h"
#endif

// TODO: This should be guarded by WITH_FORTRAN as well
#ifdef CUH2_POT
#include "eon/potentials/CuH2/CuH2.h"
#endif

#ifndef _WIN32
#ifdef WITH_VASP
#include "eon/potentials/VASP/VASP.h"
#endif
#endif

#ifdef WITH_AMS
#include "eon/potentials/AMS/AMS.h"
#include "eon/potentials/AMS_IO/AMS_IO.h"
#endif

#ifdef WITH_ASE_ORCA
#include "eon/potentials/ASE_ORCA/ASE_ORCA.h"
#endif

#ifdef WITH_ASE_NWCHEM
#include "eon/potentials/ASE_NWCHEM/ASE_NWCHEM.h"
#endif

#ifdef WITH_METATOMIC
#include "eon/potentials/Metatomic/MetatomicPotential.h"
#endif

#ifdef WITH_WATER
#include "eon/potentials/Water/Water.hpp"
#ifdef WITH_FORTRAN
#include "eon/potentials/Water_H/Tip4p_H.h"
#endif
#include "eon/potentials/Water_Pt/Tip4p_Pt.hpp"
#endif

// Should respect Fortran availability

#ifdef WITH_XTB
#include "eon/potentials/XTBPot/XTBPot.h"
#endif

#include <limits>

std::tuple<double, AtomMatrix> Potential::get_ef(const AtomMatrix &pos,
                                                 const VectorXi &atmnrs,
                                                 const Matrix3d &box) {
  double energy{std::numeric_limits<double>::infinity()};
  long nAtoms = static_cast<long>(pos.rows());
  AtomMatrix forces{MatrixXd::Zero(nAtoms, 3)};
  double var{0}; // no variance for true potentials
  this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy, &var,
              box.data());
  forceCallCounter++;
  PotRegistry::get().on_force_call(ptype);

  return std::make_tuple(energy, forces);
}

namespace eonc::helpers {
std::shared_ptr<Potential> makePotential(const Parameters &params) {
  // Inject config-file path before any potential constructor runs
  FortranPotLoader::instance().add_config_paths(
      params.potential_options.potentialsPath);
  return makePotential(params.potential_options.potential, params);
}
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         const Parameters &params) {
  // Inject config-file path before any potential constructor runs.
  // Called on every code path including Job::Job which uses this overload.
  FortranPotLoader::instance().add_config_paths(
      params.potential_options.potentialsPath);
  switch (ptype) {
  // TODO: Every potential must know their own type
  case PotType::EMT: {
    return (std::make_shared<EffectiveMediumTheory>(params));
    break;
  }
  case PotType::EXT_POT: {
    return (std::make_shared<ExtPot>(params));
    break;
  }
  case PotType::LJ: {
    return makeRgpot<rgpot::LJPot>(PotType::LJ, params, rgpot::LJConfig{});
    break;
  }
  case PotType::LJCLUSTER: {
    return makeRgpot<rgpot::LJClusterPot>(PotType::LJCLUSTER, params,
                                          rgpot::LJClusterConfig{});
    break;
  }
  case PotType::MORSE_PT: {
    return makeRgpot<rgpot::MorsePot>(PotType::MORSE_PT, params,
                                      rgpot::MorseConfig{});
    break;
  }
#ifdef NEW_POT
  case PotType::NEW: {
    return (std::make_shared<NewPot>(params));
    break;
  }
#endif
#ifdef CUH2_POT
  case PotType::CUH2: {
    return (std::make_shared<CuH2>(params));
    break;
  }
#endif
#ifdef IMD_POT
  case PotType::IMD: {
    return (std::make_shared<IMD>(params));
    break;
  }
#endif
#ifdef WITH_WATER
  case PotType::TIP4P: {
    return (std::make_shared<Tip4p>(params));
    break;
  }
  case PotType::SPCE: {
    return (std::make_shared<SpceCcl>(params));
    break;
  }
#ifdef WITH_FORTRAN
  case PotType::TIP4P_PT: {
    return (std::make_shared<Tip4p_Pt>(params));
    break;
  }
  case PotType::TIP4P_H: {
    return (std::make_shared<Tip4p_H>(params));
    break;
  }
#endif
#endif
  // Fortran potentials: always available, loaded at runtime via dlopen
  case PotType::EAM_AL: {
    return (std::make_shared<Aluminum>(params));
    break;
  }
  case PotType::EDIP: {
    return (std::make_shared<EDIP>(params));
    break;
  }
  case PotType::FEHE: {
    return (std::make_shared<FeHe>(params));
    break;
  }
  case PotType::LENOSKY_SI: {
    return (std::make_shared<Lenosky>(params));
    break;
  }
  case PotType::SW_SI: {
    return (std::make_shared<SW>(params));
    break;
  }
  case PotType::TERSOFF_SI: {
    return (std::make_shared<Tersoff>(params));
    break;
  }
#ifndef _WIN32
#ifdef WITH_VASP
  case PotType::VASP: {
    return (std::make_shared<VASP>(params));
    break;
  }
#endif
#endif
  case PotType::LAMMPS: {
    return std::make_shared<LAMMPSPot>(params);
  }
#ifdef EONMPI
  case PotType::MPI: {
    return (std::make_shared<MPIPot>(params));
    break;
  }
#endif
#ifdef EMBED_PYTHON
#ifdef PYAMFF_POT
  case PotType::PYAMFF: {
    return (std::make_shared<PyAMFF>());
    break;
  }
#endif
#ifdef WITH_ASE_POT
  case PotType::ASE_POT: {
    return (std::make_shared<ASE>(params));
    break;
  }
#endif
  // case PotType::QSC: {
  //   return (std::make_shared<QSC>());
  //   break;
  // }
#endif
#ifdef WITH_AMS
  case PotType::AMS: {
    return (std::make_shared<AMS>(params));
    break;
  }
  case PotType::AMS_IO: {
    return (std::make_shared<AMS_IO>(params));
    break;
  }
#endif
#ifdef WITH_GPRD
  // case PotType::GPR: {
  //   return "gpr"s;
  //   break;
  // }
#endif
  // case PotType::PYTHON: {
  //   TODO: Implement
  //   return "python"s;
  //   break;
  // }
#ifdef WITH_CATLEARN
  case PotType::CatLearn: {
    return (std::make_shared<CatLearnPot>(params));
    break;
  }
#endif
// TODO: Handle Fortran interaction
#ifdef WITH_XTB
  case PotType::XTB: {
    return (std::make_shared<XTBPot>(params));
    break;
  }
#endif
#ifdef WITH_ASE_ORCA
  case PotType::ASE_ORCA: {
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
  case PotType::ZBL: {
    return makeRgpot<rgpot::ZBLPot>(
        PotType::ZBL, params,
        rgpot::ZBLConfig{
            .cut_inner = params.zbl_options.cut_inner,
            .cut_global = params.zbl_options.cut_global,
        });
    break;
  }
#ifndef IS_WINDOWS
  case PotType::SocketNWChem: {
    return (std::make_shared<SocketNWChemPot>(params));
    break;
  }
#endif
#ifdef WITH_RGPOT
  case PotType::RGPOT: {
    return (std::make_shared<RgpotPot>(params));
    break;
  }
#endif
  default:
    EONC_LOG_ERROR("No known potential could be constructed from {}",
                   magic_enum::enum_name(ptype));
    eonc::log::get()->flush_log();
    throw std::runtime_error("Terminating");
    break;
  }
}

} // namespace eonc::helpers
