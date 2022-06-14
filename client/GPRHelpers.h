#ifndef GPRHELPERS_H
#define GPRHELPERS_H

#include<memory>
#include<utility>
#include <functional>

#include "Parameters.h"
#include "Log.h"
#include "HelperFunctions.h"
#include "NudgedElasticBand.h"
#include "potentials/Morse/Morse.h"
#include"potentials/GPRPotential/GPRPotential.h"

#include "subprojects/gprdimer/data_types/Coord.h"
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
     * \brief Create a parameters object for gpr_pot
     *
     * TODO: Get the cell size from a Matter object
     *
     * @param *parameters An EON parameters object
     */
    gpr::InputParameters eon_parameters_to_gprpot(Parameters *parameters);

    /**
     * \brief Create a configuration of atoms for gpr_dimer
     *
     * @param *Matter An EON Matter object
     */
    gpr::AtomsConfiguration eon_matter_to_atmconf(Matter *matter);

    /**
     * \brief Recreate the fields representing a .con file
     *
     * @param mat An EON Matter object
     */
    gpr::Field<double> generateAtomsConfigField(const Matter& mat);

    /**
     * \brief Create an initial Observation object for gpr_dimer
     *
     * Note that this is essentially only for the setup, and it does not
     * actually append or handle the Observation structure other than for
     * initialization of atomic gp dimer
     *
     * @param *Matter An EON Matter object
     */
    gpr::Observation eon_matter_to_init_obs(Matter& matter);

    /**
     * \brief Call an EON potential from a matter object
     *
     *  \param matter The object with coordinates
     *  \param pot The potential energy to use
     * */
    std::pair<double, AtomMatrix> energy_and_forces(Matter *matter, Potential *pot);
    std::pair<double, AtomMatrix> energy_and_forces_free(Matter *matter, Potential *pot);

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
    std::vector<Matter> prepInitialPath(
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

    /**
     *  \brief Checks the NEB object and returns a new observation for the GPR if required
     *
     *  This function is meant to determine if additional rounds of the GPR are
     *  required and if so, to update the observations for the same.
     *
     *  \param neb
     *  \return bool
     */
    bool maybeUpdateObs(NudgedElasticBand& neb, gpr::Observation& prevObs, Parameters& params);

    /**
     * \brief Returns a unique pointer to an NEB for the GPR-NEB
     *
     * This uses the MATLAB code logic of restarting each round with a linear interpolation
     * Effectively this ensures that the NEB is constructed on the GPR surface.
     * The unique pointer ensures in a loop, the NEB objects are destroyed
     * \param
     */
    std::unique_ptr<NudgedElasticBand> prepGPRNEBround(GPRPotential& trainedGPR, Matter& reactant, Matter& product, Parameters& params);

    /**
     * \brief Initializes a GPR potential
     *
     * The logic here is to simply prepare the GPR for training, it is a
     * super-set of EON parameter conversions
     *
     * Note that it is **initialized** but NOT trained.
     * To train this run setHyperParameters and optimize on the returned object
     * */
    gpr::GaussianProcessRegression& initializeGPR(gpr::GaussianProcessRegression& gprfunc,
                                                  gpr::AtomsConfiguration& atoms_config,
                                                  gpr::Observation& obsPath,
                                                  std::pair<Parameters, Matter>& eon_matter_params);
    /**
     * \brief Get the length of a path constructed with Matter objects
     * */
    double get_path_length(std::vector<Matter>& path);

    /**
     * \brief Generate a single coord object
     *
     * This is the conversion of one matter object to a single coordinate object
     * \note these are *all FREE* positions only
     */
     gpr::Coord single_img(const Matter& spoint);

    /**
     * \brief Generate a coord object
     *
     * The point is to get a path into a single coordinate object. One of the
     * helpers used to eventually pass / check the log based per-system image
     * ball
     * \note these are *all FREE* positions only
     */
     gpr::Coord prev_path(const std::vector<Matter>& ppath);

    /**
     * \brief Check early 1dMaxDist stopping
     *
     * Essentially a thin wrapper over the equivalent member function defined in
     * GPR dimer.
     *
     * Relevant function signatures from the GPRD are:
     *
     * bool AtomicDimer::isInterAtomicDistanceTooBig(const gpr::Coord& R_new,
     *                                        const gpr::Coord& R_all,
     *                                        gpr::Index_t& num_es1)
     *                                        ^ this one is unused
     *
     * Which in turn will call:
     *
     * void Distance::dist_max1Dlog(const gpr::Coord& x1, const gpr::Coord& x2,
     *                       const gpr::AtomsConfiguration& conf_info,
     *                       gpr::Field<double>& dist)
     **/
     bool hasEarly1DmaxStopping(const Matter& cpoint,
                                const std::vector<Matter>& ppath,
                                const gpr::AtomsConfiguration& atoms_config,
                                const double ratioAtLimit);
    } // namespace helper_functions
#endif /* GPRHELPERS_H */
