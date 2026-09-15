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
#include <cctype>
#include <csignal>
#include <ctime>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "eon/HelperFunctions.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"

namespace eonc {
Potential::Potential(PotType a_ptype, const Parameters &p) : Potential(a_ptype) {
  force_serial_ = !p.potential_options.thread_safe;
}
Potential::Potential(const Parameters &a_params)
    : Potential(a_params.potential_options.potential, a_params) {}
} // namespace eonc
#ifdef WITH_CATLEARN
#include "eon/potentials/CatLearnPot/CatLearnPot.h"
#endif

#ifdef WITH_GPRD
#include "eon/potentials/GPRPotential/GPRPotential.h"
#endif

#include "eon/potentials/EAM/EAM.h"
#include "eon/potentials/EMT/EffectiveMediumTheory.h"
#include "eon/potentials/ExtPot/ExtPot.h"
#include "eon/potentials/PluginLoader.h"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJClusterPot.hpp"
#include "rgpot/LennardJones/LJPot.hpp"
#include "rgpot/Morse/MorsePot.hpp"
#include "rgpot/ZBL/ZBLPot.hpp"
#ifdef RGPOT_HAS_DFTD3
#include "rgpot/D3Pot/D3Pot.hpp"
#endif
#ifdef RGPOT_HAS_DFTD4
#include "rgpot/D4Pot/D4Pot.hpp"
#endif
#ifdef RGPOT_HAS_EXPR
#include "rgpot/ExprPot/ExprPot.hpp"
#endif
#include "rgpot/MOPACPot/MOPACPot.hpp"
#include "rgpot/fortran/FortranPots.hpp"
#ifndef IS_WINDOWS
#include "eon/potentials/SocketNWChem/SocketNWChemPot.h"
#ifdef WITH_RGPOT
#include "eon/potentials/Rgpot/RgpotPot.h"
#endif
#endif

// Fortran potentials: always compiled, loaded at runtime via dlopen

#ifdef EMBED_PYTHON
#ifdef WITH_ASE_POT
#include "eon/potentials/ASE/ASE.h"
#endif
#endif

#ifdef EONMPI
#include "eon/potentials/MPIPot/MPIPot.h"
#endif

#include "eon/potentials/LAMMPS/LAMMPSPot.h"

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
#endif
#include "eon/potentials/Water_Pt/Tip4p_Pt.hpp"
#endif

// Should respect Fortran availability

#ifdef WITH_XTB
#include "eon/potentials/XTBPot/XTBPot.h"
#endif

#include <cmath>
#include <limits>
#include <stdexcept>

std::tuple<double, AtomMatrix> eonc::Potential::get_ef(const AtomMatrix &pos,
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
  if (!std::isfinite(energy) || !forces.allFinite()) {
    throw std::runtime_error("Potential::get_ef: non-finite energy or forces");
  }

  return std::make_tuple(energy, forces);
}

namespace eonc::helpers {
namespace {

std::string lower_copy(std::string s) {
  for (char &c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

#ifdef RGPOT_HAS_EXPR
#ifdef RGPOT_HAS_DFTD3
rgpot::D3Damping d3_damping_from_params(const Parameters &params) {
  rgpot::D3Damping damp = rgpot::D3Damping::BJ;
  if (lower_copy(params.dftd_options.d3_damping) == "zero") {
    damp = rgpot::D3Damping::Zero;
  }
  return damp;
}
#endif

std::unique_ptr<rgpot::PotentialBase> make_expr_term(const std::string &raw,
                                                     const Parameters &params) {
  const std::string name = lower_copy(raw);
  if (name == "lj") {
    return std::make_unique<rgpot::LJPot>(rgpot::LJConfig{});
  }
  if (name == "ljcluster") {
    return std::make_unique<rgpot::LJClusterPot>(rgpot::LJClusterConfig{});
  }
  if (name == "morse" || name == "morse_pt") {
    return std::make_unique<rgpot::MorsePot>(rgpot::MorseConfig{});
  }
  if (name == "zbl") {
    return std::make_unique<rgpot::ZBLPot>(rgpot::ZBLConfig{
        .cut_inner = params.zbl_options.cut_inner,
        .cut_global = params.zbl_options.cut_global,
    });
  }
#ifdef RGPOT_HAS_DFTD3
  if (name == "d3" || name == "dftd3") {
    return std::make_unique<rgpot::D3Pot>(rgpot::D3Config{
        .damping = d3_damping_from_params(params),
        .functional = params.dftd_options.functional,
        .atm = params.dftd_options.atm,
    });
  }
#endif
#ifdef RGPOT_HAS_DFTD4
  if (name == "d4" || name == "dftd4") {
    return std::make_unique<rgpot::D4Pot>(rgpot::D4Config{
        .functional = params.dftd_options.functional,
        .charge = params.dftd_options.d4_charge,
        .atm = params.dftd_options.atm,
    });
  }
#endif
  if (name == "mopac") {
    return std::make_unique<rgpot::MOPACPot>(rgpot::MOPACPot::Config{
        .charge = params.mopac_options.charge,
        .spin = params.mopac_options.spin,
        .model = params.mopac_options.model,
        .engine_path = params.mopac_options.engine_path,
    });
  }
  throw std::runtime_error(
      "ExprPot unknown term '" + raw +
      "' (lj, ljcluster, morse, zbl, d3/dftd3, d4/dftd4, mopac)");
}

std::vector<rgpot::ExprPot::Term> parse_expr_terms(const Parameters &params) {
  std::vector<rgpot::ExprPot::Term> terms;
  std::string buf = params.expr_options.terms;
  std::string token;
  auto flush = [&]() {
    while (!token.empty() &&
           std::isspace(static_cast<unsigned char>(token.front()))) {
      token.erase(token.begin());
    }
    while (!token.empty() &&
           std::isspace(static_cast<unsigned char>(token.back()))) {
      token.pop_back();
    }
    if (!token.empty()) {
      terms.emplace_back(token, make_expr_term(token, params));
      token.clear();
    }
  };
  for (char c : buf) {
    if (c == ',') {
      flush();
    } else {
      token.push_back(c);
    }
  }
  flush();
  if (terms.empty()) {
    throw std::runtime_error(
        "ExprPot needs [ExprPot] terms (comma-separated names used in "
        "expression)");
  }
  return terms;
}
#endif

} // namespace

std::shared_ptr<Potential> makePotential(const Parameters &params) {
  // Inject config-file path before any potential constructor runs
  PluginLoader::instance().add_config_paths(
      params.potential_options.potentialsPath);
  return makePotential(params.potential_options.potential, params);
}
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         const Parameters &params) {
  // Inject config-file path before any potential constructor runs.
  // Called on every code path including Job::Job which uses this overload.
  PluginLoader::instance().add_config_paths(
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
#ifdef CUH2_POT
  case PotType::CUH2: {
    return makeRgpotDefault<rgpot::fortranpots::CuH2Pot>(PotType::CUH2, params);
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
    return makeRgpotDefault<rgpot::fortranpots::WaterHPot>(PotType::TIP4P_H,
                                                           params);
    break;
  }
#endif
#endif
  // Fortran potentials: always available, loaded at runtime via dlopen
  case PotType::EAM_AL: {
    return makeRgpotDefault<rgpot::fortranpots::EAMAlPot>(PotType::EAM_AL,
                                                          params);
    break;
  }
  case PotType::EDIP: {
    return makeRgpotDefault<rgpot::fortranpots::EDIPPot>(PotType::EDIP, params);
    break;
  }
  case PotType::FEHE: {
    return makeRgpotDefault<rgpot::fortranpots::FeHePot>(PotType::FEHE, params);
    break;
  }
  case PotType::LENOSKY_SI: {
    return makeRgpotDefault<rgpot::fortranpots::LenoskyPot>(PotType::LENOSKY_SI,
                                                            params);
    break;
  }
  case PotType::SW_SI: {
    return makeRgpotDefault<rgpot::fortranpots::SWPot>(PotType::SW_SI, params);
    break;
  }
  case PotType::TERSOFF_SI: {
    return makeRgpotDefault<rgpot::fortranpots::TersoffPot>(PotType::TERSOFF_SI,
                                                            params);
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
#ifdef WITH_ASE_POT
  case PotType::ASE_POT: {
    return (std::make_shared<ASE>(params));
    break;
  }
#endif
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
#ifdef RGPOT_HAS_DFTD3
  case PotType::DFTD3: {
    rgpot::D3Damping damp = rgpot::D3Damping::BJ;
    std::string dname = params.dftd_options.d3_damping;
    for (char &c : dname) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (dname == "zero") {
      damp = rgpot::D3Damping::Zero;
    }
    return makeRgpot<rgpot::D3Pot>(
        PotType::DFTD3, params,
        rgpot::D3Config{.damping = damp,
                        .functional = params.dftd_options.functional,
                        .atm = params.dftd_options.atm});
    break;
  }
#endif
#ifdef RGPOT_HAS_DFTD4
  case PotType::DFTD4: {
    return makeRgpot<rgpot::D4Pot>(
        PotType::DFTD4, params,
        rgpot::D4Config{.functional = params.dftd_options.functional,
                        .charge = params.dftd_options.d4_charge,
                        .atm = params.dftd_options.atm});
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
  case PotType::MOPAC: {
    return makeRgpot<rgpot::MOPACPot>(
        PotType::MOPAC, params,
        rgpot::MOPACPot::Config{
            .charge = params.mopac_options.charge,
            .spin = params.mopac_options.spin,
            .model = params.mopac_options.model,
            .engine_path = params.mopac_options.engine_path,
        });
    break;
  }
#ifdef RGPOT_HAS_EXPR
  case PotType::EXPR: {
    if (params.expr_options.expression.empty()) {
      throw std::runtime_error(
          "ExprPot needs [ExprPot] expression, e.g. 0.5*lj + d3");
    }
    rgpot::ExprPot expr(params.expr_options.expression,
                        parse_expr_terms(params));
    return std::make_shared<RgpotAdapter<rgpot::ExprPot>>(PotType::EXPR, params,
                                                          std::move(expr));
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
