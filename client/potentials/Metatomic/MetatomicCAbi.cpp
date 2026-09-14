/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** C ABI surface of libmetatomic_pot.so.
*/
#define EON_MTA_BUILD
#include "eon/potentials/Metatomic/metatomic_c_abi.h"

#include "eon/Parameters.h"
#include "eon/potentials/Metatomic/MetatomicPotential.h"

#include <cstdio>
#include <exception>
#include <memory>
#include <string>

struct EonMtaPot {
  std::unique_ptr<MetatomicPotential> impl;
};

static const char *nz(const char *s) { return s ? s : ""; }

static void set_err(char *errbuf, size_t errlen, const char *msg) {
  if (!errbuf || errlen == 0)
    return;
  std::snprintf(errbuf, errlen, "%s", msg ? msg : "unknown error");
}

extern "C" {

int eon_mta_abi_version(void) { return EON_MTA_ABI_VERSION; }

EonMtaPot *eon_mta_pot_create(const EonMtaConfig *cfg, char *errbuf,
                              size_t errlen) {
  if (!cfg || !cfg->model_path || cfg->model_path[0] == '\0') {
    set_err(errbuf, errlen, "eon_mta_pot_create: model_path is required");
    return nullptr;
  }
  try {
    Parameters params;
    params.potential_options.potential = PotType::METATOMIC;
    auto &o = params.metatomic_options;
    o.model_path = cfg->model_path;
    o.device = nz(cfg->device)[0] ? cfg->device : "cpu";
    o.length_unit = nz(cfg->length_unit)[0] ? cfg->length_unit : "angstrom";
    o.extensions_directory = nz(cfg->extensions_directory);
    o.check_consistency = cfg->check_consistency != 0;
    o.uncertainty_threshold = cfg->uncertainty_threshold;
    o.energy_output = nz(cfg->energy_output);
    o.energy_uncertainty_output = nz(cfg->energy_uncertainty_output);
    o.force_output = nz(cfg->force_output);
    o.non_conservative = cfg->non_conservative != 0;
    o.random_rotation = cfg->random_rotation != 0;
    o.n_symmetry_rotations = cfg->n_symmetry_rotations;
    o.deterministic = cfg->deterministic != 0;
    o.deterministic_strict = cfg->deterministic_strict != 0;
    o.variant.base = nz(cfg->variant_base);
    o.variant.energy = nz(cfg->variant_energy);
    o.variant.energy_uncertainty = nz(cfg->variant_energy_uncertainty);
    o.variant.force = nz(cfg->variant_force);

    auto pot = std::make_unique<EonMtaPot>();
    pot->impl = std::make_unique<MetatomicPotential>(params);
    return pot.release();
  } catch (const std::exception &e) {
    set_err(errbuf, errlen, e.what());
    return nullptr;
  } catch (...) {
    set_err(errbuf, errlen, "eon_mta_pot_create: unknown exception");
    return nullptr;
  }
}

void eon_mta_pot_destroy(EonMtaPot *pot) { delete pot; }

int eon_mta_pot_force(EonMtaPot *pot, long nAtoms, const double *positions,
                      const int *atomicNrs, double *forces, double *energy,
                      double *variance, const double *box) {
  if (!pot || !pot->impl || !positions || !atomicNrs || !forces || !energy ||
      !box) {
    return 1;
  }
  try {
    pot->impl->force(nAtoms, positions, atomicNrs, forces, energy, variance,
                     box);
    return 0;
  } catch (...) {
    return 2;
  }
}

} // extern "C"
