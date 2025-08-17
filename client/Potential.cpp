#include <csignal>
#include <limits>
#include <time.h>
#include <utility>

#include "HelperFunctions.h"
#include "Parameters.h"
#include "Potential.h"
#ifdef WITH_CATLEARN
#include "potentials/CatLearnPot/CatLearnPot.h"
#endif

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
#include "potentials/ZBL/ZBLPot.h"
#include "potentials/SocketNWChem/SocketNWChemPot.h"

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
#ifdef ASE_POT
#include "potentials/ASE/ASE.h"
#endif

#include "potentials/QSC/QSC.h"
#endif

#ifdef EONMPI
#include "potentials/MPIPot/MPIPot.h"
#endif

#ifdef LAMMPS_POT
#include "potentials/LAMMPS/LAMMPSPot.h"
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

#ifdef WITH_ASE_ORCA
#include "potentials/ASE_ORCA/ASE_ORCA.h"
#endif

#ifdef WITH_ASE_NWCHEM
#include "potentials/ASE_NWCHEM/ASE_NWCHEM.h"
#endif

#ifdef WITH_METATOMIC
#include "potentials/Metatomic/MetatomicPotential.h"
#endif

#ifdef WITH_WATER
#include "potentials/Water/Water.hpp"
#ifdef WITH_FORTRAN
#include "potentials/Water_H/Tip4p_H.h"
#endif
#include "potentials/Water_Pt/Tip4p_Pt.hpp"
#endif

// Should respect Fortran availability

#ifdef WITH_XTB
#include "potentials/XTBPot/XTBPot.h"
#endif

#include <limits>

// TODO(rg): These aren't really used anymore, just there for eyecandy
int Potential::fcalls = 0;
int Potential::fcallsTotal = 0;
int Potential::wu_fcallsTotal = 0;
double Potential::totalUserTime = 0.0;

std::tuple<double, AtomMatrix> Potential::get_ef(const AtomMatrix pos,
                                                 const VectorXi atmnrs,
                                                 const Matrix3d box) {
  double energy{std::numeric_limits<double>::infinity()};
  long nAtoms{pos.rows()};
  AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
  double var{0}; // no variance for true potentials
  this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy, &var,
              box.data());
  forceCallCounter++;
  m_log->trace("[{}] {} so far", magic_enum::enum_name<PotType>(getType()),
               forceCallCounter);

  return std::make_tuple(energy, forces);
};

namespace helper_functions {
std::shared_ptr<Potential> makePotential(std::shared_ptr<Parameters> params) {
  return makePotential(params->potential, params);
}
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         std::shared_ptr<Parameters> params) {
  switch (ptype) {
  // TODO: Every potential must know their own type
  case PotType::EMT: {
    return (std::make_shared<EffectiveMediumTheory>(params));
    break;
  }
  case PotType::EXT: {
    return (std::make_shared<ExtPot>(params));
    break;
  }
  case PotType::LJ: {
    return (std::make_shared<LJ>(params));
    break;
  }
  case PotType::LJCLUSTER: {
    return (std::make_shared<LJCluster>(params));
    break;
  }
  case PotType::MORSE_PT: {
    return (std::make_shared<Morse>(params));
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
#ifdef WITH_FORTRAN
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
#endif
#ifndef WIN32
#ifdef WITH_VASP
  case PotType::VASP: {
    return (std::make_shared<VASP>(params));
    break;
  }
#endif
#endif
#ifdef LAMMPS_POT
  case PotType::LAMMPS: {
    return (std::make_shared<LAMMPSPot>(params));
    break;
  }
#endif
#ifdef EONMPI
  case PotType::MPI: {
    return (std::make_shared<MPIPot>(params));
    break;
  }
#endif
#ifdef WITH_PYTHON
#ifdef PYAMFF_POT
  case PotType::PYAMFF: {
    return (std::make_shared<PyAMFF>());
    break;
  }
#endif
#ifdef ASE_POT
  case PotType::ASE_POT: {
    return (std::make_shared<ASE_POT>(params));
    break;
  }
#endif
  // case PotType::QSC: {
  //   return (std::make_shared<QSC>());
  //   break;
  // }
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
    return (std::make_shared<ZBLPot>(params));
    break;
  }
  case PotType::SocketNWChem: {
    return (std::make_shared<SocketNWChemPot>(params));
    break;
  }
  default:
    SPDLOG_ERROR("No known potential could be constructed from {}",
                 magic_enum::enum_name(ptype));
    throw std::runtime_error("Terminating");
    break;
  }
}

} // namespace helper_functions
