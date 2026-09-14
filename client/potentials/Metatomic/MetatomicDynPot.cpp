/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/potentials/Metatomic/MetatomicDynPot.h"
#include "eon/potentials/Metatomic/MetatomicLoader.h"

#include <array>
#include <stdexcept>
#include <string>

namespace eonc {

static EonMtaConfig config_from_params(const Parameters &params) {
  const auto &o = params.metatomic_options;
  EonMtaConfig c{};
  c.model_path = o.model_path.c_str();
  c.device = o.device.c_str();
  c.length_unit = o.length_unit.c_str();
  c.extensions_directory = o.extensions_directory.c_str();
  c.check_consistency = o.check_consistency ? 1 : 0;
  c.uncertainty_threshold = o.uncertainty_threshold;
  c.energy_output = o.energy_output.c_str();
  c.energy_uncertainty_output = o.energy_uncertainty_output.c_str();
  c.force_output = o.force_output.c_str();
  c.non_conservative = o.non_conservative ? 1 : 0;
  c.random_rotation = o.random_rotation ? 1 : 0;
  c.n_symmetry_rotations = o.n_symmetry_rotations;
  c.deterministic = o.deterministic ? 1 : 0;
  c.deterministic_strict = o.deterministic_strict ? 1 : 0;
  c.variant_base = o.variant.base.c_str();
  c.variant_energy = o.variant.energy.c_str();
  c.variant_energy_uncertainty = o.variant.energy_uncertainty.c_str();
  c.variant_force = o.variant.force.c_str();
  return c;
}

MetatomicDynPot::MetatomicDynPot(const Parameters &params)
    : Potential(PotType::METATOMIC) {
  auto &loader = MetatomicLoader::instance();
  loader.require_loaded();
  std::array<char, 1024> err{};
  auto cfg = config_from_params(params);
  m_handle = loader.create(&cfg, err.data(), err.size());
  if (!m_handle) {
    throw std::runtime_error(std::string("MetatomicDynPot: create failed: ") +
                             err.data());
  }
}

MetatomicDynPot::~MetatomicDynPot() {
  if (m_handle) {
    MetatomicLoader::instance().destroy(m_handle);
    m_handle = nullptr;
  }
}

void MetatomicDynPot::force(long nAtoms, const double *positions,
                            const int *atomicNrs, double *forces,
                            double *energy, double *variance,
                            const double *box) {
  auto &loader = MetatomicLoader::instance();
  const int rc = loader.force(m_handle, nAtoms, positions, atomicNrs, forces,
                              energy, variance, box);
  if (rc != 0) {
    throw std::runtime_error("MetatomicDynPot: eon_mta_pot_force failed rc=" +
                             std::to_string(rc));
  }
}

} // namespace eonc
