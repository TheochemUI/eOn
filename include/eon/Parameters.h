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

#include "BaseStructures.h"
#include "ParametersOptions.h"
#include <cstdio>
#include <string>

/** Contains all runtime parameters and results. No functionality just
 * bookkeeping.*/
namespace eonc {

class Parameters {

public:
  Parameters();
  ~Parameters() = default;
  Parameters(const Parameters &) = default;
  int load(std::string filename);
  int load(FILE *file);
  int load_ini_text(const std::string &ini_text);
  int load_json(const std::string &json_str);
  std::string to_json() const;

  using constants_t = eonc::constants_t;
  using main_options_t = eonc::main_options_t;
  using potential_options_t = eonc::potential_options_t;
  using ams_options_t = eonc::ams_options_t;
  using xtb_options_t = eonc::xtb_options_t;
  using zbl_options_t = eonc::zbl_options_t;
  using dftd_options_t = eonc::dftd_options_t;
  using expr_options_t = eonc::expr_options_t;
  using mopac_options_t = eonc::mopac_options_t;
  using socket_nwchem_options_t = eonc::socket_nwchem_options_t;
  using rgpot_options_t = eonc::rgpot_options_t;
  using structure_comparison_options_t = eonc::structure_comparison_options_t;
  using process_search_options_t = eonc::process_search_options_t;
  using saddle_search_options_t = eonc::saddle_search_options_t;
  using optimizer_options_t = eonc::optimizer_options_t;
  using dimer_options_t = eonc::dimer_options_t;
  using gpr_dimer_options_t = eonc::gpr_dimer_options_t;
  using gp_surrogate_options_t = eonc::gp_surrogate_options_t;
  using catlearn_options_t = eonc::catlearn_options_t;
  using ase_orca_options_t = eonc::ase_orca_options_t;
  using ase_nwchem_options_t = eonc::ase_nwchem_options_t;
  using metatomic_options_t = eonc::metatomic_options_t;
  using lanczos_options_t = eonc::lanczos_options_t;
  using davidson_options_t = eonc::davidson_options_t;
  using prefactor_options_t = eonc::prefactor_options_t;
  using hessian_options_t = eonc::hessian_options_t;
  using neb_options_t = eonc::neb_options_t;
  using dynamics_options_t = eonc::dynamics_options_t;
  using parallel_replica_options_t = eonc::parallel_replica_options_t;
  using tad_options_t = eonc::tad_options_t;
  using thermostat_options_t = eonc::thermostat_options_t;
  using replica_exchange_options_t = eonc::replica_exchange_options_t;
  using hyperdynamics_options_t = eonc::hyperdynamics_options_t;
  using basin_hopping_options_t = eonc::basin_hopping_options_t;
  using global_optimization_options_t = eonc::global_optimization_options_t;
  using monte_carlo_options_t = eonc::monte_carlo_options_t;
  using bgsd_options_t = eonc::bgsd_options_t;
  using serve_options_t = eonc::serve_options_t;
  using artn_options_t = eonc::artn_options_t;
  using ira_options_t = eonc::ira_options_t;
  using debug_options_t = eonc::debug_options_t;
  using oh_tst_options_t = eonc::oh_tst_options_t;

  constants_t constants;
  main_options_t main_options;
  potential_options_t potential_options;
  ams_options_t ams_options;
  xtb_options_t xtb_options;
  zbl_options_t zbl_options;
  dftd_options_t dftd_options;
  expr_options_t expr_options;
  mopac_options_t mopac_options;
  socket_nwchem_options_t socket_nwchem_options;
  rgpot_options_t rgpot_options;
  structure_comparison_options_t structure_comparison_options;
  process_search_options_t process_search_options;
  saddle_search_options_t saddle_search_options;
  optimizer_options_t optimizer_options;
  dimer_options_t dimer_options;
  gpr_dimer_options_t gpr_dimer_options;
  gp_surrogate_options_t gp_surrogate_options;
  catlearn_options_t catlearn_options;
  ase_orca_options_t ase_orca_options;
  ase_nwchem_options_t ase_nwchem_options;
  metatomic_options_t metatomic_options;
  lanczos_options_t lanczos_options;
  davidson_options_t davidson_options;
  prefactor_options_t prefactor_options;
  hessian_options_t hessian_options;
  neb_options_t neb_options;
  dynamics_options_t dynamics_options;
  parallel_replica_options_t parallel_replica_options;
  tad_options_t tad_options;
  thermostat_options_t thermostat_options;
  replica_exchange_options_t replica_exchange_options;
  hyperdynamics_options_t hyperdynamics_options;
  basin_hopping_options_t basin_hopping_options;
  global_optimization_options_t global_optimization_options;
  monte_carlo_options_t monte_carlo_options;
  bgsd_options_t bgsd_options;
  serve_options_t serve_options;
  artn_options_t artn_options;
  ira_options_t ira_options;
  debug_options_t debug_options;
  oh_tst_options_t oh_tst_options;
};

} // namespace eonc
