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

#include "XTBPot.h"
#include "XtbLoader.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

using forcefields::unit_system::BOHR;
using forcefields::unit_system::HARTREE;

namespace {
// libxtb's verbosity sentinel: must match XTB_VERBOSITY_MUTED in xtb.h
// (= 0). We re-define here so this translation unit doesn't pull in
// xtb.h; the loader's runtime ABI is plain int.
constexpr int kXtbVerbosityMuted = 0;
} // namespace

XTBPot::XTBPot(const Parameters &p)
    : Potential(PotType::XTB, p),
      xtb_acc{p.xtb_options.acc},
      xtb_electronic_temperature{p.xtb_options.elec_temperature},
      xtb_max_iter{p.xtb_options.maxiter},
      total_charge{p.xtb_options.charge},
      uhf{p.xtb_options.uhf} {
  auto &xtb = eonc::get_xtb_loader();
  xtb.require_loaded();

  env = xtb.new_environment();
  if (!env) {
    throw std::runtime_error("Failed to create xtb environment");
  }
  xtb.set_verbosity(env, kXtbVerbosityMuted);
  // Release the default output unit to prevent Fortran NEWUNIT
  // conflicts when multiple XTB environments coexist (per-image NEB
  // potentials, e.g.).
  xtb.release_output(env);

  calc = xtb.new_calculator();
  if (!calc) {
    xtb.del_environment(&env);
    throw std::runtime_error("Failed to create xtb calculator");
  }
  res = xtb.new_results();
  if (!res) {
    xtb.del_calculator(&calc);
    xtb.del_environment(&env);
    throw std::runtime_error("Failed to create xtb results");
  }

  if (p.xtb_options.paramset == "GFNFF") {
    xtb_paramset = GFNMethod::GFNFF;
  } else if (p.xtb_options.paramset == "GFN0xTB") {
    xtb_paramset = GFNMethod::GFN0xTB;
  } else if (p.xtb_options.paramset == "GFN1xTB") {
    xtb_paramset = GFNMethod::GFN1xTB;
  } else if (p.xtb_options.paramset == "GFN2xTB") {
    xtb_paramset = GFNMethod::GFN2xTB;
  } else {
    throw std::runtime_error("Parameter set for XTB must be one of GFNFF, "
                             "GFN0xTB, GFN1xTB or GFN2xTB.\n");
  }
}

XTBPot::~XTBPot() {
  // Loader instance is process-singleton; survives every potential.
  auto &xtb = eonc::get_xtb_loader();
  if (!xtb.is_loaded()) {
    return; // dlopen failed at construction; nothing to free.
  }
  if (res) {
    xtb.del_results(&res);
  }
  if (calc) {
    xtb.del_calculator(&calc);
  }
  if (mol) {
    xtb.del_molecule(&mol);
  }
  if (env) {
    xtb.release_output(env);
    xtb.del_environment(&env);
  }
  QUILL_LOG_INFO(eonc::log::get(), "[XTB] called potential {} times",
                 counter++);
}

void XTBPot::cleanMemory() {}

void XTBPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  auto &xtb = eonc::get_xtb_loader();

  int intN = static_cast<int>(N);
  // TODO: Periodicity shouldn't crash
  const bool periodicity[3]{false, false, false};
  double box_bohr[3 * 3];

  // Convert positions and lattice from Angstrom to Bohr.
  std::vector<double> R_bohr(3 * N);
  for (long idx = 0; idx < 3 * N; ++idx) {
    R_bohr[idx] = R[idx] / BOHR;
  }
  for (long idx = 0; idx < 9; ++idx) {
    box_bohr[idx] = box[idx] / BOHR;
  }

  if (!initialized) {
    mol = xtb.new_molecule(env, &intN, atomicNrs, R_bohr.data(), &total_charge,
                           &uhf, box_bohr, periodicity);

    switch (xtb_paramset) {
    case GFNMethod::GFNFF:
      if (!xtb.load_gfnff) {
        throw std::runtime_error(
            "libxtb installed but xtb_loadGFNFF is not exported (older or "
            "stripped build); pick a different paramset or upgrade xtb.");
      }
      xtb.load_gfnff(env, mol, calc, nullptr);
      break;
    case GFNMethod::GFN0xTB:
      if (!xtb.load_gfn0) {
        throw std::runtime_error(
            "libxtb installed but xtb_loadGFN0xTB is not exported.");
      }
      xtb.load_gfn0(env, mol, calc, nullptr);
      break;
    case GFNMethod::GFN1xTB:
      if (!xtb.load_gfn1) {
        throw std::runtime_error(
            "libxtb installed but xtb_loadGFN1xTB is not exported.");
      }
      xtb.load_gfn1(env, mol, calc, nullptr);
      break;
    case GFNMethod::GFN2xTB:
      if (!xtb.load_gfn2) {
        throw std::runtime_error(
            "libxtb installed but xtb_loadGFN2xTB is not exported.");
      }
      xtb.load_gfn2(env, mol, calc, nullptr);
      break;
    }

    xtb.set_accuracy(env, calc, xtb_acc);
    xtb.set_electronic_temp(env, calc, xtb_electronic_temperature);
    xtb.set_max_iter(env, calc, static_cast<int>(xtb_max_iter));
    initialized = true;
  } else {
    xtb.update_molecule(env, mol, R_bohr.data(), box_bohr);
  }

  xtb.singlepoint(env, mol, calc, res);

  if (xtb.check_environment(env) != 0) {
    char err_msg[512]{};
    int buf = sizeof(err_msg);
    xtb.get_error(env, err_msg, &buf);
    throw std::runtime_error(std::string("xTB Error: ") + err_msg);
  }

  xtb.get_energy(env, res, U);
  xtb.get_gradient(env, res, F);

  // Convert Hartree / Bohr -> eV / Angstrom and flip sign (xtb returns
  // gradient, eOn wants force).
  for (long i = 0; i < 3 * N; ++i) {
    F[i] *= -1.0 * (HARTREE / BOHR);
  }
  *U *= HARTREE;
  counter++;
}
