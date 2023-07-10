#ifndef GPRHELPERS_H
#define GPRHELPERS_H

#include "Matter.h"
#include "Parameters.h"

#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gpr_optim/structures/Structures.h"

namespace helpers::gproptim::input {
/**
 * \brief Create a parameters object for gpr_dimer
 *
 * TODO: Get the cell size from a Matter object
 *
 * @param *parameters An EON parameters object
 */
gpr::InputParameters
eon_parameters_to_gpr(std::shared_ptr<Parameters> a_params);

/**
 * \brief Create a configuration of atoms for gpr_dimer
 *
 * @param *Matter An EON Matter object
 */
gpr::AtomsConfiguration eon_matter_to_atmconf(std::shared_ptr<Matter> a_matter);

/**
 * \brief Create an initial Observation object for gpr_dimer
 *
 * Note that this is essentially only for the setup, and it does not
 * actually append or handle the Observation structure other than for
 * initialization of atomic gp dimer
 *
 * @param *Matter An EON Matter object
 */
gpr::Observation eon_matter_to_init_obs(std::shared_ptr<Matter> a_matter);

/**
 * \brief Prepare the atoms configuration data
 */
std::pair<gpr::AtomsConfiguration, gpr::Coord>
eon_matter_to_frozen_conf_info(std::shared_ptr<Matter> a_matter,
                               double a_activeRadius);

gpr::Field<double> generateAtomsConfigField(const Matter &mat);

} // namespace helpers::gproptim::input
#endif /* GPRHELPERS_H */
