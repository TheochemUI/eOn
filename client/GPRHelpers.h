#ifndef GPRHELPERS_H
#define GPRHELPERS_H

#include<memory>
#include<utility>

#include "Matter.h"
#include "Parameters.h"
#include "Log.h"
#include "HelperFunctions.h"
#include "potentials/GPRPotential/GPRPotential.h"

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
    gpr::InputParameters eon_parameters_to_gprd(Parameters *parameters);

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
     * \brief Call a GPR-EON potential from a matter object
     *
     *  \param matter The object with coordinates
     *  \param pot The trained GPR to use
     * */
    std::pair<double, AtomMatrix> gpr_energy_and_forces(Matter* matter, GPRPotential* gprpot);

    /**
     * \brief Generate a conf_info object from a matter object
     *
     * This actually calls eon_matter_to_atmconf under the hood
     * The difference being that this function sets the frozen atoms
     *  \param matter The object to use as a constructor
     * */
    std::pair<gpr::AtomsConfiguration, gpr::Coord> eon_matter_to_frozen_conf_info(Matter *matter, double activeRadius);

    // SCG Helpers
    struct MatterHolder {
            // TODO: Make private
            Matter* mrr;
            void getEnergyGradient(const Eigen::VectorXd& w,
                                   const gpr::EigenMatrix& x,
                                   const Eigen::VectorXd& x_ind,
                                   const Eigen::VectorXd& y,
                                   gpr::EnergyAndGradient& energy_and_gradient);
    };

    /**
     * \brief Generate a conf_info object from a matter object
     *
     * This actually calls eon_matter_to_atmconf under the hood
     * The difference being that this function sets the frozen atoms
     *  \param matter The object to use as a constructor
     * */

    /**
     * \brief Setup initial path
     *
     * This is essentially what the NEB initialization does, however, we need
     * this here to prepare the initial observations for the GPR surface
     * */
    std::pair<std::vector<Matter>,
               std::vector<AtomMatrix> > prepInitialPath(
                   Parameters *params,
                   std::string fname_reactant="reactant.con"s,
                   std::string fname_product="product.con"s);
    /**
     * \brief Setup initial observations from path
     *
     * This starts with the construction of an observation object and then uses
     * the linearly interpolated NEB path from prepInitialPath to populate the
     * observation object
     * */
    gpr::Observation prepInitialObs(std::vector<Matter> &vecmat);
    } // namespace helper_functions
#endif /* GPRHELPERS_H */
