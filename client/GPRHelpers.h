#ifndef GPRHELPERS_H
#define GPRHELPERS_H

#include "Matter.h"
#include "Parameters.h"
#include "Log.h"

#include "subprojects/gprdimer/structures/Structures.h"
#include "subprojects/gprdimer/gpr/auxiliary/ProblemSetUp.h"

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

    /**
     * \brief Call an EON potential from a matter object
     *
     *  \param matter The object with coordinates
     *  \param pot The potential energy to use
     * */
    std::pair<double, AtomMatrix> energy_and_forces(Matter *matter, Potential *pot);

    /**
     * \brief Generate a conf_info object from a matter object
     *
     * This actually calls eon_matter_to_atmconf under the hood
     * The difference being that this function sets the frozen atoms
     *  \param matter The object to use as a constructor
     * */
    std::pair<gpr::AtomsConfiguration, gpr::Coord> eon_matter_to_frozen_conf_info(Matter *matter, double activeRadius);

    } // namespace helper_functions
#endif /* GPRHELPERS_H */
