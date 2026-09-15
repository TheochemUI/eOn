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

  void set_mpi_client_comm(std::uintptr_t raw) {
    potential_options_.MPIClientComm = raw;
  }
  std::uintptr_t mpi_client_comm() const {
    return potential_options_.MPIClientComm;
  }
  void set_mpi_potential_rank(int rank) {
    potential_options_.MPIPotentialRank = rank;
  }

  friend int config::load_ini(::INIReader &, Parameters &);
  friend void config::apply_ssot_defaults(Parameters &);
  friend void config::validate_and_link(Parameters &);
  friend struct ParametersLoadAccess;

  const constants_t &constants() const { return constants_; }
  const main_options_t &main_options() const { return main_options_; }
  const potential_options_t &potential_options() const {
    return potential_options_;
  }
  const ams_options_t &ams_options() const { return ams_options_; }
  const xtb_options_t &xtb_options() const { return xtb_options_; }
  const zbl_options_t &zbl_options() const { return zbl_options_; }
  const dftd_options_t &dftd_options() const { return dftd_options_; }
  const expr_options_t &expr_options() const { return expr_options_; }
  const mopac_options_t &mopac_options() const { return mopac_options_; }
  const socket_nwchem_options_t &socket_nwchem_options() const {
    return socket_nwchem_options_;
  }
  const rgpot_options_t &rgpot_options() const { return rgpot_options_; }
  const structure_comparison_options_t &structure_comparison_options() const {
    return structure_comparison_options_;
  }
  const process_search_options_t &process_search_options() const {
    return process_search_options_;
  }
  const saddle_search_options_t &saddle_search_options() const {
    return saddle_search_options_;
  }
  const optimizer_options_t &optimizer_options() const {
    return optimizer_options_;
  }
  const dimer_options_t &dimer_options() const { return dimer_options_; }
  const gpr_dimer_options_t &gpr_dimer_options() const {
    return gpr_dimer_options_;
  }
  const gp_surrogate_options_t &gp_surrogate_options() const {
    return gp_surrogate_options_;
  }
  const catlearn_options_t &catlearn_options() const {
    return catlearn_options_;
  }
  const ase_orca_options_t &ase_orca_options() const {
    return ase_orca_options_;
  }
  const ase_nwchem_options_t &ase_nwchem_options() const {
    return ase_nwchem_options_;
  }
  const metatomic_options_t &metatomic_options() const {
    return metatomic_options_;
  }
  const lanczos_options_t &lanczos_options() const { return lanczos_options_; }
  const davidson_options_t &davidson_options() const {
    return davidson_options_;
  }
  const prefactor_options_t &prefactor_options() const {
    return prefactor_options_;
  }
  const hessian_options_t &hessian_options() const { return hessian_options_; }
  const neb_options_t &neb_options() const { return neb_options_; }
  const dynamics_options_t &dynamics_options() const {
    return dynamics_options_;
  }
  const parallel_replica_options_t &parallel_replica_options() const {
    return parallel_replica_options_;
  }
  const tad_options_t &tad_options() const { return tad_options_; }
  const thermostat_options_t &thermostat_options() const {
    return thermostat_options_;
  }
  const replica_exchange_options_t &replica_exchange_options() const {
    return replica_exchange_options_;
  }
  const hyperdynamics_options_t &hyperdynamics_options() const {
    return hyperdynamics_options_;
  }
  const basin_hopping_options_t &basin_hopping_options() const {
    return basin_hopping_options_;
  }
  const global_optimization_options_t &global_optimization_options() const {
    return global_optimization_options_;
  }
  const monte_carlo_options_t &monte_carlo_options() const {
    return monte_carlo_options_;
  }
  const bgsd_options_t &bgsd_options() const { return bgsd_options_; }
  const serve_options_t &serve_options() const { return serve_options_; }
  const artn_options_t &artn_options() const { return artn_options_; }
  const ira_options_t &ira_options() const { return ira_options_; }
  const debug_options_t &debug_options() const { return debug_options_; }
  const oh_tst_options_t &oh_tst_options() const { return oh_tst_options_; }

private:
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

  /// Load-state pimpl (source + last error). Option-group layout stays in
  /// this header and is not ABI-stable; see include/eon/api.h.
  struct Impl;
  std::unique_ptr<Impl> impl_;
  void record_load(std::string_view source, int error);
};

/// Write hole for INI/JSON loaders and nanobind property setters.
struct ParametersLoadAccess {
  static constants_t &constants(Parameters &p) { return p.constants_; }
  static const constants_t &constants(const Parameters &p) {
    return p.constants_;
  }
  static main_options_t &main_options(Parameters &p) { return p.main_options_; }
  static const main_options_t &main_options(const Parameters &p) {
    return p.main_options_;
  }
  static potential_options_t &potential_options(Parameters &p) {
    return p.potential_options_;
  }
  static const potential_options_t &potential_options(const Parameters &p) {
    return p.potential_options_;
  }
  static ams_options_t &ams_options(Parameters &p) { return p.ams_options_; }
  static const ams_options_t &ams_options(const Parameters &p) {
    return p.ams_options_;
  }
  static xtb_options_t &xtb_options(Parameters &p) { return p.xtb_options_; }
  static const xtb_options_t &xtb_options(const Parameters &p) {
    return p.xtb_options_;
  }
  static zbl_options_t &zbl_options(Parameters &p) { return p.zbl_options_; }
  static const zbl_options_t &zbl_options(const Parameters &p) {
    return p.zbl_options_;
  }
  static dftd_options_t &dftd_options(Parameters &p) { return p.dftd_options_; }
  static const dftd_options_t &dftd_options(const Parameters &p) {
    return p.dftd_options_;
  }
  static expr_options_t &expr_options(Parameters &p) { return p.expr_options_; }
  static const expr_options_t &expr_options(const Parameters &p) {
    return p.expr_options_;
  }
  static mopac_options_t &mopac_options(Parameters &p) {
    return p.mopac_options_;
  }
  static const mopac_options_t &mopac_options(const Parameters &p) {
    return p.mopac_options_;
  }
  static socket_nwchem_options_t &socket_nwchem_options(Parameters &p) {
    return p.socket_nwchem_options_;
  }
  static const socket_nwchem_options_t &
  socket_nwchem_options(const Parameters &p) {
    return p.socket_nwchem_options_;
  }
  static rgpot_options_t &rgpot_options(Parameters &p) {
    return p.rgpot_options_;
  }
  static const rgpot_options_t &rgpot_options(const Parameters &p) {
    return p.rgpot_options_;
  }
  static structure_comparison_options_t &
  structure_comparison_options(Parameters &p) {
    return p.structure_comparison_options_;
  }
  static const structure_comparison_options_t &
  structure_comparison_options(const Parameters &p) {
    return p.structure_comparison_options_;
  }
  static process_search_options_t &process_search_options(Parameters &p) {
    return p.process_search_options_;
  }
  static const process_search_options_t &
  process_search_options(const Parameters &p) {
    return p.process_search_options_;
  }
  static saddle_search_options_t &saddle_search_options(Parameters &p) {
    return p.saddle_search_options_;
  }
  static const saddle_search_options_t &
  saddle_search_options(const Parameters &p) {
    return p.saddle_search_options_;
  }
  static optimizer_options_t &optimizer_options(Parameters &p) {
    return p.optimizer_options_;
  }
  static const optimizer_options_t &optimizer_options(const Parameters &p) {
    return p.optimizer_options_;
  }
  static dimer_options_t &dimer_options(Parameters &p) {
    return p.dimer_options_;
  }
  static const dimer_options_t &dimer_options(const Parameters &p) {
    return p.dimer_options_;
  }
  static gpr_dimer_options_t &gpr_dimer_options(Parameters &p) {
    return p.gpr_dimer_options_;
  }
  static const gpr_dimer_options_t &gpr_dimer_options(const Parameters &p) {
    return p.gpr_dimer_options_;
  }
  static gp_surrogate_options_t &gp_surrogate_options(Parameters &p) {
    return p.gp_surrogate_options_;
  }
  static const gp_surrogate_options_t &
  gp_surrogate_options(const Parameters &p) {
    return p.gp_surrogate_options_;
  }
  static catlearn_options_t &catlearn_options(Parameters &p) {
    return p.catlearn_options_;
  }
  static const catlearn_options_t &catlearn_options(const Parameters &p) {
    return p.catlearn_options_;
  }
  static ase_orca_options_t &ase_orca_options(Parameters &p) {
    return p.ase_orca_options_;
  }
  static const ase_orca_options_t &ase_orca_options(const Parameters &p) {
    return p.ase_orca_options_;
  }
  static ase_nwchem_options_t &ase_nwchem_options(Parameters &p) {
    return p.ase_nwchem_options_;
  }
  static const ase_nwchem_options_t &ase_nwchem_options(const Parameters &p) {
    return p.ase_nwchem_options_;
  }
  static metatomic_options_t &metatomic_options(Parameters &p) {
    return p.metatomic_options_;
  }
  static const metatomic_options_t &metatomic_options(const Parameters &p) {
    return p.metatomic_options_;
  }
  static lanczos_options_t &lanczos_options(Parameters &p) {
    return p.lanczos_options_;
  }
  static const lanczos_options_t &lanczos_options(const Parameters &p) {
    return p.lanczos_options_;
  }
  static davidson_options_t &davidson_options(Parameters &p) {
    return p.davidson_options_;
  }
  static const davidson_options_t &davidson_options(const Parameters &p) {
    return p.davidson_options_;
  }
  static prefactor_options_t &prefactor_options(Parameters &p) {
    return p.prefactor_options_;
  }
  static const prefactor_options_t &prefactor_options(const Parameters &p) {
    return p.prefactor_options_;
  }
  static hessian_options_t &hessian_options(Parameters &p) {
    return p.hessian_options_;
  }
  static const hessian_options_t &hessian_options(const Parameters &p) {
    return p.hessian_options_;
  }
  static neb_options_t &neb_options(Parameters &p) { return p.neb_options_; }
  static const neb_options_t &neb_options(const Parameters &p) {
    return p.neb_options_;
  }
  static dynamics_options_t &dynamics_options(Parameters &p) {
    return p.dynamics_options_;
  }
  static const dynamics_options_t &dynamics_options(const Parameters &p) {
    return p.dynamics_options_;
  }
  static parallel_replica_options_t &parallel_replica_options(Parameters &p) {
    return p.parallel_replica_options_;
  }
  static const parallel_replica_options_t &
  parallel_replica_options(const Parameters &p) {
    return p.parallel_replica_options_;
  }
  static tad_options_t &tad_options(Parameters &p) { return p.tad_options_; }
  static const tad_options_t &tad_options(const Parameters &p) {
    return p.tad_options_;
  }
  static thermostat_options_t &thermostat_options(Parameters &p) {
    return p.thermostat_options_;
  }
  static const thermostat_options_t &thermostat_options(const Parameters &p) {
    return p.thermostat_options_;
  }
  static replica_exchange_options_t &replica_exchange_options(Parameters &p) {
    return p.replica_exchange_options_;
  }
  static const replica_exchange_options_t &
  replica_exchange_options(const Parameters &p) {
    return p.replica_exchange_options_;
  }
  static hyperdynamics_options_t &hyperdynamics_options(Parameters &p) {
    return p.hyperdynamics_options_;
  }
  static const hyperdynamics_options_t &
  hyperdynamics_options(const Parameters &p) {
    return p.hyperdynamics_options_;
  }
  static basin_hopping_options_t &basin_hopping_options(Parameters &p) {
    return p.basin_hopping_options_;
  }
  static const basin_hopping_options_t &
  basin_hopping_options(const Parameters &p) {
    return p.basin_hopping_options_;
  }
  static global_optimization_options_t &
  global_optimization_options(Parameters &p) {
    return p.global_optimization_options_;
  }
  static const global_optimization_options_t &
  global_optimization_options(const Parameters &p) {
    return p.global_optimization_options_;
  }
  static monte_carlo_options_t &monte_carlo_options(Parameters &p) {
    return p.monte_carlo_options_;
  }
  static const monte_carlo_options_t &monte_carlo_options(const Parameters &p) {
    return p.monte_carlo_options_;
  }
  static bgsd_options_t &bgsd_options(Parameters &p) { return p.bgsd_options_; }
  static const bgsd_options_t &bgsd_options(const Parameters &p) {
    return p.bgsd_options_;
  }
  static serve_options_t &serve_options(Parameters &p) {
    return p.serve_options_;
  }
  static const serve_options_t &serve_options(const Parameters &p) {
    return p.serve_options_;
  }
  static artn_options_t &artn_options(Parameters &p) { return p.artn_options_; }
  static const artn_options_t &artn_options(const Parameters &p) {
    return p.artn_options_;
  }
  static ira_options_t &ira_options(Parameters &p) { return p.ira_options_; }
  static const ira_options_t &ira_options(const Parameters &p) {
    return p.ira_options_;
  }
  static debug_options_t &debug_options(Parameters &p) {
    return p.debug_options_;
  }
  static const debug_options_t &debug_options(const Parameters &p) {
    return p.debug_options_;
  }
  static oh_tst_options_t &oh_tst_options(Parameters &p) {
    return p.oh_tst_options_;
  }
  static const oh_tst_options_t &oh_tst_options(const Parameters &p) {
    return p.oh_tst_options_;
  }
};

} // namespace eonc
