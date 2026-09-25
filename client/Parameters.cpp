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
#include "eon/Parameters.h"
#include "eon/ParametersINI.h"
#include "eon/ParametersJSON.h"
#include "eon/ParametersSSOT.h"
#include "magic_enum/magic_enum.hpp"

#include <INIReader.h>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <utility>

#include "eon/EonLogger.h"

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

static_assert(sizeof(Parameters) == sizeof(std::unique_ptr<int>),
              "Parameters layout is the Impl pointer");

Parameters::Parameters()
    : impl_(std::make_unique<Impl>()) {
  // Covered groups: defaults originate from schema/eon_params.capnp via
  // apply_ssot_defaults (codegen). Uncovered groups still use NSDMI.
  eonc::config::apply_ssot_defaults(*this);
  // Resolve computed fields via validate_and_link.
  eonc::config::validate_and_link(*this);
}

Parameters::~Parameters() = default;
Parameters::Parameters(Parameters &&) noexcept = default;
Parameters &Parameters::operator=(Parameters &&) noexcept = default;

Parameters::Parameters(const Parameters &other)
    : impl_(std::make_unique<Impl>(other.impl_ ? *other.impl_ : Impl{})) {}

Parameters &Parameters::operator=(const Parameters &other) {
  if (this == &other) {
    return *this;
  }
  Parameters tmp(other);
  *this = std::move(tmp);
  return *this;
}

std::string_view Parameters::last_load_source() const {
  return impl_ ? std::string_view{impl_->last_source} : std::string_view{};
}

int Parameters::last_load_error() const {
  return impl_ ? impl_->last_error : 0;
}

void Parameters::record_load(std::string_view source, int error) {
  Impl &impl = ensure_impl();
  impl.last_source.assign(source);
  impl.last_error = error;
}

int Parameters::load(std::string_view filename) {
  INIReader ini{std::string(filename)};
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Can't load INI file: {}", filename);
    record_load(filename, 1);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);

  // Sanity Checks
  if (impl_->parallel_replica_options_.state_check_interval >
          impl_->dynamics_options_.time &&
      magic_enum::enum_name<JobType>(impl_->main_options_.job) ==
          "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }

  if (!impl_->neb_options_.initialization.input_path.empty() &&
      impl_->neb_options_.initialization.method == NEBInit::LINEAR) {
    EONC_LOG_WARNING(
        "[Nudged Elastic Band] 'initial_path_in' is provided, but "
        "'initializer' defaults to linear. "
        "Ensure this is intentional, as the loaded path will not be "
        "used without initializer set to file.");
  }

  if (impl_->saddle_search_options_.dynamics.record_interval_input >
      impl_->saddle_search_options_.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }

  if (impl_->potential_options_.potential == PotType::AMS ||
      impl_->potential_options_.potential == PotType::AMS_IO) {
    if (impl_->ams_options_.forcefield.empty() &&
        impl_->ams_options_.model.empty() && impl_->ams_options_.xc.empty()) {
      EONC_LOG_ERROR("[AMS] Must provide atleast forcefield or model or xc");
      error = 1;
    }

    if (!impl_->ams_options_.forcefield.empty() &&
        !impl_->ams_options_.model.empty() && !impl_->ams_options_.xc.empty()) {
      EONC_LOG_ERROR("[AMS] Must provide either forcefield or model");
      error = 1;
    }
  }

  record_load(filename, error);
  return error;
}

int Parameters::load(FILE *file) {
  constexpr std::string_view kFileSource{"<FILE*>"};
  if (!file) {
    EONC_LOG_ERROR("Can't load INI from a null FILE*");
    record_load(kFileSource, 1);
    return 1;
  }
  if (fseek(file, 0, SEEK_END) != 0) {
    EONC_LOG_ERROR("Can't seek INI FILE*");
    record_load(kFileSource, 1);
    return 1;
  }
  const long size = ftell(file);
  if (size < 0) {
    EONC_LOG_ERROR("Can't tell INI FILE* size");
    record_load(kFileSource, 1);
    return 1;
  }
  if (fseek(file, 0, SEEK_SET) != 0) {
    EONC_LOG_ERROR("Can't rewind INI FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  std::string buffer(static_cast<size_t>(size), '\0');
  if (fread(buffer.data(), 1, static_cast<size_t>(size), file) !=
      static_cast<size_t>(size)) {
    EONC_LOG_ERROR("Couldn't read the ini file from FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  INIReader ini(buffer.c_str(), buffer.size());
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Couldn't parse the ini file from FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);
  // Same validation as the filename overload
  // (duplicated intentionally for now; will be extracted to validate())
  if (impl_->parallel_replica_options_.state_check_interval >
          impl_->dynamics_options_.time &&
      magic_enum::enum_name<JobType>(impl_->main_options_.job) ==
          "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }
  if (impl_->saddle_search_options_.dynamics.record_interval_input >
      impl_->saddle_search_options_.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }
  record_load(kFileSource, error);
  return error;
}

int Parameters::load_ini_text(std::string_view ini_text) {
  constexpr std::string_view kIniSource{"<ini>"};
  INIReader ini(ini_text.data(), ini_text.size());
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Couldn't parse INI from memory");
    record_load(kIniSource, 1);
    return 1;
  }
  const int error = eonc::config::load_ini(ini, *this);
  record_load(kIniSource, error);
  return error;
}

int Parameters::load_json(std::string_view json_str) {
  constexpr std::string_view kJsonSource{"<json>"};
  const int error = eonc::config::load_json(json_str, *this);
  record_load(kJsonSource, error);
  return error;
}

std::string Parameters::to_json() const {
  return eonc::config::to_json(*this).dump(2);
}

Parameters::Impl &Parameters::ensure_impl() {
  if (!impl_) {
    impl_ = std::make_unique<Impl>();
  }
  return *impl_;
}

const Parameters::constants_t &Parameters::constants() const {
  return impl_->constants_;
}

const Parameters::main_options_t &Parameters::main_options() const {
  return impl_->main_options_;
}

const Parameters::potential_options_t &Parameters::potential_options() const {
  return impl_->potential_options_;
}

const Parameters::ams_options_t &Parameters::ams_options() const {
  return impl_->ams_options_;
}

const Parameters::xtb_options_t &Parameters::xtb_options() const {
  return impl_->xtb_options_;
}

const Parameters::zbl_options_t &Parameters::zbl_options() const {
  return impl_->zbl_options_;
}

const Parameters::dftd_options_t &Parameters::dftd_options() const {
  return impl_->dftd_options_;
}

const Parameters::expr_options_t &Parameters::expr_options() const {
  return impl_->expr_options_;
}

const Parameters::mopac_options_t &Parameters::mopac_options() const {
  return impl_->mopac_options_;
}

const Parameters::socket_nwchem_options_t &
Parameters::socket_nwchem_options() const {
  return impl_->socket_nwchem_options_;
}

const Parameters::rgpot_options_t &Parameters::rgpot_options() const {
  return impl_->rgpot_options_;
}

const Parameters::structure_comparison_options_t &
Parameters::structure_comparison_options() const {
  return impl_->structure_comparison_options_;
}

const Parameters::process_search_options_t &
Parameters::process_search_options() const {
  return impl_->process_search_options_;
}

const Parameters::saddle_search_options_t &
Parameters::saddle_search_options() const {
  return impl_->saddle_search_options_;
}

const Parameters::optimizer_options_t &Parameters::optimizer_options() const {
  return impl_->optimizer_options_;
}

const Parameters::dimer_options_t &Parameters::dimer_options() const {
  return impl_->dimer_options_;
}

const Parameters::gpr_dimer_options_t &Parameters::gpr_dimer_options() const {
  return impl_->gpr_dimer_options_;
}

const Parameters::gp_surrogate_options_t &
Parameters::gp_surrogate_options() const {
  return impl_->gp_surrogate_options_;
}

const Parameters::catlearn_options_t &Parameters::catlearn_options() const {
  return impl_->catlearn_options_;
}

const Parameters::ase_orca_options_t &Parameters::ase_orca_options() const {
  return impl_->ase_orca_options_;
}

const Parameters::ase_nwchem_options_t &Parameters::ase_nwchem_options() const {
  return impl_->ase_nwchem_options_;
}

const Parameters::metatomic_options_t &Parameters::metatomic_options() const {
  return impl_->metatomic_options_;
}

const Parameters::lanczos_options_t &Parameters::lanczos_options() const {
  return impl_->lanczos_options_;
}

const Parameters::davidson_options_t &Parameters::davidson_options() const {
  return impl_->davidson_options_;
}

const Parameters::prefactor_options_t &Parameters::prefactor_options() const {
  return impl_->prefactor_options_;
}

const Parameters::hessian_options_t &Parameters::hessian_options() const {
  return impl_->hessian_options_;
}

const Parameters::neb_options_t &Parameters::neb_options() const {
  return impl_->neb_options_;
}

const Parameters::dynamics_options_t &Parameters::dynamics_options() const {
  return impl_->dynamics_options_;
}

const Parameters::parallel_replica_options_t &
Parameters::parallel_replica_options() const {
  return impl_->parallel_replica_options_;
}

const Parameters::tad_options_t &Parameters::tad_options() const {
  return impl_->tad_options_;
}

const Parameters::thermostat_options_t &Parameters::thermostat_options() const {
  return impl_->thermostat_options_;
}

const Parameters::replica_exchange_options_t &
Parameters::replica_exchange_options() const {
  return impl_->replica_exchange_options_;
}

const Parameters::hyperdynamics_options_t &
Parameters::hyperdynamics_options() const {
  return impl_->hyperdynamics_options_;
}

const Parameters::basin_hopping_options_t &
Parameters::basin_hopping_options() const {
  return impl_->basin_hopping_options_;
}

const Parameters::global_optimization_options_t &
Parameters::global_optimization_options() const {
  return impl_->global_optimization_options_;
}

const Parameters::monte_carlo_options_t &
Parameters::monte_carlo_options() const {
  return impl_->monte_carlo_options_;
}

const Parameters::bgsd_options_t &Parameters::bgsd_options() const {
  return impl_->bgsd_options_;
}

const Parameters::serve_options_t &Parameters::serve_options() const {
  return impl_->serve_options_;
}

const Parameters::artn_options_t &Parameters::artn_options() const {
  return impl_->artn_options_;
}

const Parameters::ira_options_t &Parameters::ira_options() const {
  return impl_->ira_options_;
}

const Parameters::debug_options_t &Parameters::debug_options() const {
  return impl_->debug_options_;
}

const Parameters::oh_tst_options_t &Parameters::oh_tst_options() const {
  return impl_->oh_tst_options_;
}

void Parameters::set_mpi_client_comm(std::uintptr_t raw) {
  ensure_impl().potential_options_.MPIClientComm = raw;
}

std::uintptr_t Parameters::mpi_client_comm() const {
  return impl_->potential_options_.MPIClientComm;
}

void Parameters::set_mpi_potential_rank(int rank) {
  ensure_impl().potential_options_.MPIPotentialRank = rank;
}

constants_t &ParametersLoadAccess::constants(Parameters &p) {
  return p.ensure_impl().constants_;
}

const constants_t &ParametersLoadAccess::constants(const Parameters &p) {
  return p.impl_->constants_;
}

main_options_t &ParametersLoadAccess::main_options(Parameters &p) {
  return p.ensure_impl().main_options_;
}

const main_options_t &ParametersLoadAccess::main_options(const Parameters &p) {
  return p.impl_->main_options_;
}

potential_options_t &ParametersLoadAccess::potential_options(Parameters &p) {
  return p.ensure_impl().potential_options_;
}

const potential_options_t &
ParametersLoadAccess::potential_options(const Parameters &p) {
  return p.impl_->potential_options_;
}

ams_options_t &ParametersLoadAccess::ams_options(Parameters &p) {
  return p.ensure_impl().ams_options_;
}

const ams_options_t &ParametersLoadAccess::ams_options(const Parameters &p) {
  return p.impl_->ams_options_;
}

xtb_options_t &ParametersLoadAccess::xtb_options(Parameters &p) {
  return p.ensure_impl().xtb_options_;
}

const xtb_options_t &ParametersLoadAccess::xtb_options(const Parameters &p) {
  return p.impl_->xtb_options_;
}

zbl_options_t &ParametersLoadAccess::zbl_options(Parameters &p) {
  return p.ensure_impl().zbl_options_;
}

const zbl_options_t &ParametersLoadAccess::zbl_options(const Parameters &p) {
  return p.impl_->zbl_options_;
}

dftd_options_t &ParametersLoadAccess::dftd_options(Parameters &p) {
  return p.ensure_impl().dftd_options_;
}

const dftd_options_t &ParametersLoadAccess::dftd_options(const Parameters &p) {
  return p.impl_->dftd_options_;
}

expr_options_t &ParametersLoadAccess::expr_options(Parameters &p) {
  return p.ensure_impl().expr_options_;
}

const expr_options_t &ParametersLoadAccess::expr_options(const Parameters &p) {
  return p.impl_->expr_options_;
}

mopac_options_t &ParametersLoadAccess::mopac_options(Parameters &p) {
  return p.ensure_impl().mopac_options_;
}

const mopac_options_t &
ParametersLoadAccess::mopac_options(const Parameters &p) {
  return p.impl_->mopac_options_;
}

socket_nwchem_options_t &
ParametersLoadAccess::socket_nwchem_options(Parameters &p) {
  return p.ensure_impl().socket_nwchem_options_;
}

const socket_nwchem_options_t &
ParametersLoadAccess::socket_nwchem_options(const Parameters &p) {
  return p.impl_->socket_nwchem_options_;
}

rgpot_options_t &ParametersLoadAccess::rgpot_options(Parameters &p) {
  return p.ensure_impl().rgpot_options_;
}

const rgpot_options_t &
ParametersLoadAccess::rgpot_options(const Parameters &p) {
  return p.impl_->rgpot_options_;
}

structure_comparison_options_t &
ParametersLoadAccess::structure_comparison_options(Parameters &p) {
  return p.ensure_impl().structure_comparison_options_;
}

const structure_comparison_options_t &
ParametersLoadAccess::structure_comparison_options(const Parameters &p) {
  return p.impl_->structure_comparison_options_;
}

process_search_options_t &
ParametersLoadAccess::process_search_options(Parameters &p) {
  return p.ensure_impl().process_search_options_;
}

const process_search_options_t &
ParametersLoadAccess::process_search_options(const Parameters &p) {
  return p.impl_->process_search_options_;
}

saddle_search_options_t &
ParametersLoadAccess::saddle_search_options(Parameters &p) {
  return p.ensure_impl().saddle_search_options_;
}

const saddle_search_options_t &
ParametersLoadAccess::saddle_search_options(const Parameters &p) {
  return p.impl_->saddle_search_options_;
}

optimizer_options_t &ParametersLoadAccess::optimizer_options(Parameters &p) {
  return p.ensure_impl().optimizer_options_;
}

const optimizer_options_t &
ParametersLoadAccess::optimizer_options(const Parameters &p) {
  return p.impl_->optimizer_options_;
}

dimer_options_t &ParametersLoadAccess::dimer_options(Parameters &p) {
  return p.ensure_impl().dimer_options_;
}

const dimer_options_t &
ParametersLoadAccess::dimer_options(const Parameters &p) {
  return p.impl_->dimer_options_;
}

gpr_dimer_options_t &ParametersLoadAccess::gpr_dimer_options(Parameters &p) {
  return p.ensure_impl().gpr_dimer_options_;
}

const gpr_dimer_options_t &
ParametersLoadAccess::gpr_dimer_options(const Parameters &p) {
  return p.impl_->gpr_dimer_options_;
}

gp_surrogate_options_t &
ParametersLoadAccess::gp_surrogate_options(Parameters &p) {
  return p.ensure_impl().gp_surrogate_options_;
}

const gp_surrogate_options_t &
ParametersLoadAccess::gp_surrogate_options(const Parameters &p) {
  return p.impl_->gp_surrogate_options_;
}

catlearn_options_t &ParametersLoadAccess::catlearn_options(Parameters &p) {
  return p.ensure_impl().catlearn_options_;
}

const catlearn_options_t &
ParametersLoadAccess::catlearn_options(const Parameters &p) {
  return p.impl_->catlearn_options_;
}

ase_orca_options_t &ParametersLoadAccess::ase_orca_options(Parameters &p) {
  return p.ensure_impl().ase_orca_options_;
}

const ase_orca_options_t &
ParametersLoadAccess::ase_orca_options(const Parameters &p) {
  return p.impl_->ase_orca_options_;
}

ase_nwchem_options_t &ParametersLoadAccess::ase_nwchem_options(Parameters &p) {
  return p.ensure_impl().ase_nwchem_options_;
}

const ase_nwchem_options_t &
ParametersLoadAccess::ase_nwchem_options(const Parameters &p) {
  return p.impl_->ase_nwchem_options_;
}

metatomic_options_t &ParametersLoadAccess::metatomic_options(Parameters &p) {
  return p.ensure_impl().metatomic_options_;
}

const metatomic_options_t &
ParametersLoadAccess::metatomic_options(const Parameters &p) {
  return p.impl_->metatomic_options_;
}

lanczos_options_t &ParametersLoadAccess::lanczos_options(Parameters &p) {
  return p.ensure_impl().lanczos_options_;
}

const lanczos_options_t &
ParametersLoadAccess::lanczos_options(const Parameters &p) {
  return p.impl_->lanczos_options_;
}

davidson_options_t &ParametersLoadAccess::davidson_options(Parameters &p) {
  return p.ensure_impl().davidson_options_;
}

const davidson_options_t &
ParametersLoadAccess::davidson_options(const Parameters &p) {
  return p.impl_->davidson_options_;
}

prefactor_options_t &ParametersLoadAccess::prefactor_options(Parameters &p) {
  return p.ensure_impl().prefactor_options_;
}

const prefactor_options_t &
ParametersLoadAccess::prefactor_options(const Parameters &p) {
  return p.impl_->prefactor_options_;
}

hessian_options_t &ParametersLoadAccess::hessian_options(Parameters &p) {
  return p.ensure_impl().hessian_options_;
}

const hessian_options_t &
ParametersLoadAccess::hessian_options(const Parameters &p) {
  return p.impl_->hessian_options_;
}

neb_options_t &ParametersLoadAccess::neb_options(Parameters &p) {
  return p.ensure_impl().neb_options_;
}

const neb_options_t &ParametersLoadAccess::neb_options(const Parameters &p) {
  return p.impl_->neb_options_;
}

dynamics_options_t &ParametersLoadAccess::dynamics_options(Parameters &p) {
  return p.ensure_impl().dynamics_options_;
}

const dynamics_options_t &
ParametersLoadAccess::dynamics_options(const Parameters &p) {
  return p.impl_->dynamics_options_;
}

parallel_replica_options_t &
ParametersLoadAccess::parallel_replica_options(Parameters &p) {
  return p.ensure_impl().parallel_replica_options_;
}

const parallel_replica_options_t &
ParametersLoadAccess::parallel_replica_options(const Parameters &p) {
  return p.impl_->parallel_replica_options_;
}

tad_options_t &ParametersLoadAccess::tad_options(Parameters &p) {
  return p.ensure_impl().tad_options_;
}

const tad_options_t &ParametersLoadAccess::tad_options(const Parameters &p) {
  return p.impl_->tad_options_;
}

thermostat_options_t &ParametersLoadAccess::thermostat_options(Parameters &p) {
  return p.ensure_impl().thermostat_options_;
}

const thermostat_options_t &
ParametersLoadAccess::thermostat_options(const Parameters &p) {
  return p.impl_->thermostat_options_;
}

replica_exchange_options_t &
ParametersLoadAccess::replica_exchange_options(Parameters &p) {
  return p.ensure_impl().replica_exchange_options_;
}

const replica_exchange_options_t &
ParametersLoadAccess::replica_exchange_options(const Parameters &p) {
  return p.impl_->replica_exchange_options_;
}

hyperdynamics_options_t &
ParametersLoadAccess::hyperdynamics_options(Parameters &p) {
  return p.ensure_impl().hyperdynamics_options_;
}

const hyperdynamics_options_t &
ParametersLoadAccess::hyperdynamics_options(const Parameters &p) {
  return p.impl_->hyperdynamics_options_;
}

basin_hopping_options_t &
ParametersLoadAccess::basin_hopping_options(Parameters &p) {
  return p.ensure_impl().basin_hopping_options_;
}

const basin_hopping_options_t &
ParametersLoadAccess::basin_hopping_options(const Parameters &p) {
  return p.impl_->basin_hopping_options_;
}

global_optimization_options_t &
ParametersLoadAccess::global_optimization_options(Parameters &p) {
  return p.ensure_impl().global_optimization_options_;
}

const global_optimization_options_t &
ParametersLoadAccess::global_optimization_options(const Parameters &p) {
  return p.impl_->global_optimization_options_;
}

monte_carlo_options_t &
ParametersLoadAccess::monte_carlo_options(Parameters &p) {
  return p.ensure_impl().monte_carlo_options_;
}

const monte_carlo_options_t &
ParametersLoadAccess::monte_carlo_options(const Parameters &p) {
  return p.impl_->monte_carlo_options_;
}

bgsd_options_t &ParametersLoadAccess::bgsd_options(Parameters &p) {
  return p.ensure_impl().bgsd_options_;
}

const bgsd_options_t &ParametersLoadAccess::bgsd_options(const Parameters &p) {
  return p.impl_->bgsd_options_;
}

serve_options_t &ParametersLoadAccess::serve_options(Parameters &p) {
  return p.ensure_impl().serve_options_;
}

const serve_options_t &
ParametersLoadAccess::serve_options(const Parameters &p) {
  return p.impl_->serve_options_;
}

artn_options_t &ParametersLoadAccess::artn_options(Parameters &p) {
  return p.ensure_impl().artn_options_;
}

const artn_options_t &ParametersLoadAccess::artn_options(const Parameters &p) {
  return p.impl_->artn_options_;
}

ira_options_t &ParametersLoadAccess::ira_options(Parameters &p) {
  return p.ensure_impl().ira_options_;
}

const ira_options_t &ParametersLoadAccess::ira_options(const Parameters &p) {
  return p.impl_->ira_options_;
}

debug_options_t &ParametersLoadAccess::debug_options(Parameters &p) {
  return p.ensure_impl().debug_options_;
}

const debug_options_t &
ParametersLoadAccess::debug_options(const Parameters &p) {
  return p.impl_->debug_options_;
}

oh_tst_options_t &ParametersLoadAccess::oh_tst_options(Parameters &p) {
  return p.ensure_impl().oh_tst_options_;
}

const oh_tst_options_t &
ParametersLoadAccess::oh_tst_options(const Parameters &p) {
  return p.impl_->oh_tst_options_;
}

} // namespace eonc
