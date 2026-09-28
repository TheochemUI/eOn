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
#include "ParametersImpl.h"

namespace eonc {

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

} // namespace eonc
