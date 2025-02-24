#ifndef GPRHELPERS_H
#define GPRHELPERS_H

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
#endif /* GPRHELPERS_H */
