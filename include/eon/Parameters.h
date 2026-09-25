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
#include <memory>
#include <string>
#include <string_view>

class INIReader;

/** Contains all runtime parameters and results. No functionality just
 * bookkeeping.*/
namespace eonc {

class Parameters;
namespace config {
int load_ini(::INIReader &, Parameters &);
void apply_ssot_defaults(Parameters &);
void validate_and_link(Parameters &);
} // namespace config
struct ParametersLoadAccess;

class Parameters {

public:
  Parameters();
  ~Parameters();
  Parameters(const Parameters &);
  Parameters(Parameters &&) noexcept;
  Parameters &operator=(const Parameters &);
  Parameters &operator=(Parameters &&) noexcept;
  int load(std::string_view filename);
  int load(FILE *file);
  int load_ini_text(std::string_view ini_text);
  int load_json(std::string_view json_str);
  std::string to_json() const;

  /// Last load source recorded by load / load_ini_text / load_json.
  /// Empty before the first load. Moved-from objects return empty.
  [[nodiscard]] std::string_view last_load_source() const;
  /// Last load return code (0 ok). Zero before the first load.
  [[nodiscard]] int last_load_error() const;

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

  void set_mpi_client_comm(std::uintptr_t raw);
  [[nodiscard]] std::uintptr_t mpi_client_comm() const;
  void set_mpi_potential_rank(int rank);

  const constants_t &constants() const;
  const main_options_t &main_options() const;
  const potential_options_t &potential_options() const;
  const ams_options_t &ams_options() const;
  const xtb_options_t &xtb_options() const;
  const zbl_options_t &zbl_options() const;
  const dftd_options_t &dftd_options() const;
  const expr_options_t &expr_options() const;
  const mopac_options_t &mopac_options() const;
  const socket_nwchem_options_t &socket_nwchem_options() const;
  const rgpot_options_t &rgpot_options() const;
  const structure_comparison_options_t &structure_comparison_options() const;
  const process_search_options_t &process_search_options() const;
  const saddle_search_options_t &saddle_search_options() const;
  const optimizer_options_t &optimizer_options() const;
  const dimer_options_t &dimer_options() const;
  const gpr_dimer_options_t &gpr_dimer_options() const;
  const gp_surrogate_options_t &gp_surrogate_options() const;
  const catlearn_options_t &catlearn_options() const;
  const ase_orca_options_t &ase_orca_options() const;
  const ase_nwchem_options_t &ase_nwchem_options() const;
  const metatomic_options_t &metatomic_options() const;
  const lanczos_options_t &lanczos_options() const;
  const davidson_options_t &davidson_options() const;
  const prefactor_options_t &prefactor_options() const;
  const hessian_options_t &hessian_options() const;
  const neb_options_t &neb_options() const;
  const dynamics_options_t &dynamics_options() const;
  const parallel_replica_options_t &parallel_replica_options() const;
  const tad_options_t &tad_options() const;
  const thermostat_options_t &thermostat_options() const;
  const replica_exchange_options_t &replica_exchange_options() const;
  const hyperdynamics_options_t &hyperdynamics_options() const;
  const basin_hopping_options_t &basin_hopping_options() const;
  const global_optimization_options_t &global_optimization_options() const;
  const monte_carlo_options_t &monte_carlo_options() const;
  const bgsd_options_t &bgsd_options() const;
  const serve_options_t &serve_options() const;
  const artn_options_t &artn_options() const;
  const ira_options_t &ira_options() const;
  const debug_options_t &debug_options() const;
  const oh_tst_options_t &oh_tst_options() const;

  friend int config::load_ini(::INIReader &, Parameters &);
  friend void config::apply_ssot_defaults(Parameters &);
  friend void config::validate_and_link(Parameters &);
  friend struct ParametersLoadAccess;

private:
  /// Option groups and load state. sizeof(Parameters) is the
  /// pointer; see include/eon/api.h.
  struct Impl;
  std::unique_ptr<Impl> impl_;
  Impl &ensure_impl();
  void record_load(std::string_view source, int error);
};

/// Write hole for INI/JSON loaders and nanobind property setters.
struct ParametersLoadAccess {
  static constants_t &constants(Parameters &p);
  static const constants_t &constants(const Parameters &p);
  static main_options_t &main_options(Parameters &p);
  static const main_options_t &main_options(const Parameters &p);
  static potential_options_t &potential_options(Parameters &p);
  static const potential_options_t &potential_options(const Parameters &p);
  static ams_options_t &ams_options(Parameters &p);
  static const ams_options_t &ams_options(const Parameters &p);
  static xtb_options_t &xtb_options(Parameters &p);
  static const xtb_options_t &xtb_options(const Parameters &p);
  static zbl_options_t &zbl_options(Parameters &p);
  static const zbl_options_t &zbl_options(const Parameters &p);
  static dftd_options_t &dftd_options(Parameters &p);
  static const dftd_options_t &dftd_options(const Parameters &p);
  static expr_options_t &expr_options(Parameters &p);
  static const expr_options_t &expr_options(const Parameters &p);
  static mopac_options_t &mopac_options(Parameters &p);
  static const mopac_options_t &mopac_options(const Parameters &p);
  static socket_nwchem_options_t &socket_nwchem_options(Parameters &p);
  static const socket_nwchem_options_t &
  socket_nwchem_options(const Parameters &p);
  static rgpot_options_t &rgpot_options(Parameters &p);
  static const rgpot_options_t &rgpot_options(const Parameters &p);
  static structure_comparison_options_t &
  structure_comparison_options(Parameters &p);
  static const structure_comparison_options_t &
  structure_comparison_options(const Parameters &p);
  static process_search_options_t &process_search_options(Parameters &p);
  static const process_search_options_t &
  process_search_options(const Parameters &p);
  static saddle_search_options_t &saddle_search_options(Parameters &p);
  static const saddle_search_options_t &
  saddle_search_options(const Parameters &p);
  static optimizer_options_t &optimizer_options(Parameters &p);
  static const optimizer_options_t &optimizer_options(const Parameters &p);
  static dimer_options_t &dimer_options(Parameters &p);
  static const dimer_options_t &dimer_options(const Parameters &p);
  static gpr_dimer_options_t &gpr_dimer_options(Parameters &p);
  static const gpr_dimer_options_t &gpr_dimer_options(const Parameters &p);
  static gp_surrogate_options_t &gp_surrogate_options(Parameters &p);
  static const gp_surrogate_options_t &
  gp_surrogate_options(const Parameters &p);
  static catlearn_options_t &catlearn_options(Parameters &p);
  static const catlearn_options_t &catlearn_options(const Parameters &p);
  static ase_orca_options_t &ase_orca_options(Parameters &p);
  static const ase_orca_options_t &ase_orca_options(const Parameters &p);
  static ase_nwchem_options_t &ase_nwchem_options(Parameters &p);
  static const ase_nwchem_options_t &ase_nwchem_options(const Parameters &p);
  static metatomic_options_t &metatomic_options(Parameters &p);
  static const metatomic_options_t &metatomic_options(const Parameters &p);
  static lanczos_options_t &lanczos_options(Parameters &p);
  static const lanczos_options_t &lanczos_options(const Parameters &p);
  static davidson_options_t &davidson_options(Parameters &p);
  static const davidson_options_t &davidson_options(const Parameters &p);
  static prefactor_options_t &prefactor_options(Parameters &p);
  static const prefactor_options_t &prefactor_options(const Parameters &p);
  static hessian_options_t &hessian_options(Parameters &p);
  static const hessian_options_t &hessian_options(const Parameters &p);
  static neb_options_t &neb_options(Parameters &p);
  static const neb_options_t &neb_options(const Parameters &p);
  static dynamics_options_t &dynamics_options(Parameters &p);
  static const dynamics_options_t &dynamics_options(const Parameters &p);
  static parallel_replica_options_t &parallel_replica_options(Parameters &p);
  static const parallel_replica_options_t &
  parallel_replica_options(const Parameters &p);
  static tad_options_t &tad_options(Parameters &p);
  static const tad_options_t &tad_options(const Parameters &p);
  static thermostat_options_t &thermostat_options(Parameters &p);
  static const thermostat_options_t &thermostat_options(const Parameters &p);
  static replica_exchange_options_t &replica_exchange_options(Parameters &p);
  static const replica_exchange_options_t &
  replica_exchange_options(const Parameters &p);
  static hyperdynamics_options_t &hyperdynamics_options(Parameters &p);
  static const hyperdynamics_options_t &
  hyperdynamics_options(const Parameters &p);
  static basin_hopping_options_t &basin_hopping_options(Parameters &p);
  static const basin_hopping_options_t &
  basin_hopping_options(const Parameters &p);
  static global_optimization_options_t &
  global_optimization_options(Parameters &p);
  static const global_optimization_options_t &
  global_optimization_options(const Parameters &p);
  static monte_carlo_options_t &monte_carlo_options(Parameters &p);
  static const monte_carlo_options_t &monte_carlo_options(const Parameters &p);
  static bgsd_options_t &bgsd_options(Parameters &p);
  static const bgsd_options_t &bgsd_options(const Parameters &p);
  static serve_options_t &serve_options(Parameters &p);
  static const serve_options_t &serve_options(const Parameters &p);
  static artn_options_t &artn_options(Parameters &p);
  static const artn_options_t &artn_options(const Parameters &p);
  static ira_options_t &ira_options(Parameters &p);
  static const ira_options_t &ira_options(const Parameters &p);
  static debug_options_t &debug_options(Parameters &p);
  static const debug_options_t &debug_options(const Parameters &p);
  static oh_tst_options_t &oh_tst_options(Parameters &p);
  static const oh_tst_options_t &oh_tst_options(const Parameters &p);
};

} // namespace eonc
