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
#pragma once

#include "eon/Parameters.h"

#include <string>

namespace eonc {

struct Parameters::Impl {
  std::string last_source;
  int last_error{0};
  constants_t constants_{};
  main_options_t main_options_{};
  potential_options_t potential_options_{};
  ams_options_t ams_options_{};
  xtb_options_t xtb_options_{};
  zbl_options_t zbl_options_{};
  dftd_options_t dftd_options_{};
  expr_options_t expr_options_{};
  mopac_options_t mopac_options_{};
  socket_nwchem_options_t socket_nwchem_options_{};
  rgpot_options_t rgpot_options_{};
  structure_comparison_options_t structure_comparison_options_{};
  process_search_options_t process_search_options_{};
  saddle_search_options_t saddle_search_options_{};
  optimizer_options_t optimizer_options_{};
  dimer_options_t dimer_options_{};
  gpr_dimer_options_t gpr_dimer_options_{};
  gp_surrogate_options_t gp_surrogate_options_{};
  catlearn_options_t catlearn_options_{};
  ase_orca_options_t ase_orca_options_{};
  ase_nwchem_options_t ase_nwchem_options_{};
  metatomic_options_t metatomic_options_{};
  lanczos_options_t lanczos_options_{};
  davidson_options_t davidson_options_{};
  prefactor_options_t prefactor_options_{};
  hessian_options_t hessian_options_{};
  neb_options_t neb_options_{};
  dynamics_options_t dynamics_options_{};
  parallel_replica_options_t parallel_replica_options_{};
  tad_options_t tad_options_{};
  thermostat_options_t thermostat_options_{};
  replica_exchange_options_t replica_exchange_options_{};
  hyperdynamics_options_t hyperdynamics_options_{};
  basin_hopping_options_t basin_hopping_options_{};
  global_optimization_options_t global_optimization_options_{};
  monte_carlo_options_t monte_carlo_options_{};
  bgsd_options_t bgsd_options_{};
  serve_options_t serve_options_{};
  artn_options_t artn_options_{};
  ira_options_t ira_options_{};
  debug_options_t debug_options_{};
  oh_tst_options_t oh_tst_options_{};
};

} // namespace eonc
