#include "Matter.h"
#include "Potential.h"
#include "potentials/Morse/Morse.h"
#include "potentials/GPRPotential/GPRPotential.h"

#include "HelperFunctions.h"
#include "GPRHelpers.h"
#include "Parameters.h"

#include <memory>

class GPRobj {
    private:
        Parameters eonp;
        gpr::AtomsConfiguration atmconf;
        gpr::GPRSetup gpr_parameters;
        gpr::Observation prepobs(std::vector<Matter>& matvec, Potential* pot);
    public:
        gpr::GaussianProcessRegression gprfunc; // TODO: fix yieldGPRPot to make this private
        gpr::Observation curpath;
        GPRobj(Matter initMatter, Parameters eonp);
        ~GPRobj(); // Destructor
        void trainGPR(std::vector<Matter>& initialPoints);
        void retrainGPR(std::vector<Matter>& newPoints);
        GPRPotential yieldGPRPot();
};

class GPRMatter {
    private:
        AtomMatrix trueForcesFree;
        AtomMatrix trueForces;
        double truePotEnergy;
        std::shared_ptr<GPRobj> gprobj;
    public:
        Matter truePotMatter; // TODO: Make private
        GPRMatter(Matter& initMatter, std::shared_ptr<GPRobj> gpf);
        ~GPRMatter(); // Destructor
        std::pair<double, AtomMatrix> gpr_energy_forces();
        std::pair<double, AtomMatrix> true_free_energy_forces();
        std::pair<double, AtomMatrix> true_energy_forces();
        bool areEnergiesCloseToTrue(double eps = 1e-3);
        bool areForcesCloseToTrue(double eps = 1e-3);
        void updateMatter(Matter& otherMatter);
};

namespace helper_functions {
        bool maybeUpdateGPRobj(NudgedElasticBand& neb, std::shared_ptr<GPRobj>& gpf);
        /**
         * \brief Compares energies of GPRMatter objects
         *
         * This is one of a few comparison functions to use with <algorithm>
         */
        bool max_gpr_energy(GPRMatter matOne, GPRMatter matTwo);
} // namespace helper_functions
