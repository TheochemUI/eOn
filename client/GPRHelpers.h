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
#include "Matter.h"
#include "Parameters.h"

#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gpr_optim/structures/Structures.h"

namespace helper_functions {
/**
 * \brief Create a parameters object for gpr_dimer
 *
 * TODO: Get the cell size from a Matter object
 *
 * @param *parameters An EON parameters object
 */
gpr::InputParameters eon_parameters_to_gpr(Parameters *parameters);

/**
 * \brief Create a configuration of atoms for gpr_dimer
 *
 * @param *Matter An EON Matter object
 */
gpr::AtomsConfiguration eon_matter_to_atmconf(Matter *matter);

/**
 * \brief Create an initial Observation object for gpr_dimer
 *
 * Note that this is essentially only for the setup, and it does not
 * actually append or handle the Observation structure other than for
 * initialization of atomic gp dimer
 *
 * @param *Matter An EON Matter object
 */
gpr::Observation eon_matter_to_init_obs(Matter *matter);

} // namespace helper_functions
