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
// Parameters storage and accessors. Built into eoncbase, so Potential
// constructors and plugin libraries read options without linking eonclib.
// Construction with schema defaults and the INI/JSON loaders stay in
// Parameters.cpp.
#include "ParametersImpl.h"

#include <memory>
#include <string_view>
#include <utility>

namespace eonc {

static_assert(sizeof(Parameters) == sizeof(std::unique_ptr<int>),
              "Parameters layout is the Impl pointer");

Parameters::Parameters(FieldDefaults)
    : impl_(std::make_unique<Impl>()) {}

Parameters ParametersLoadAccess::field_defaults() {
  return Parameters(Parameters::FieldDefaults{});
}

Parameters::~Parameters() = default;
Parameters::Parameters(Parameters &&) noexcept = default;
Parameters &Parameters::operator=(Parameters &&) noexcept = default;

Parameters::Parameters(const Parameters &other)
    : impl_(std::make_unique<Impl>(other.impl_ ? *other.impl_ : Impl{})) {}

// Copy in place, so option-group references taken before the assignment
// stay valid.
Parameters &Parameters::operator=(const Parameters &other) {
  if (this == &other) {
    return *this;
  }
  if (!other.impl_) {
    impl_.reset();
  } else if (impl_) {
    *impl_ = *other.impl_;
  } else {
    impl_ = std::make_unique<Impl>(*other.impl_);
  }
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

const Parameters::instanton_options_t &Parameters::instanton_options() const {
  return impl_->instanton_options_;
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

instanton_options_t &ParametersLoadAccess::instanton_options(Parameters &p) {
  return p.ensure_impl().instanton_options_;
}

const instanton_options_t &
ParametersLoadAccess::instanton_options(const Parameters &p) {
  return p.impl_->instanton_options_;
}

} // namespace eonc
