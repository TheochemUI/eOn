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
        GPRPotential trainedGPR;
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
        Matter truePotMatter;
        std::shared_ptr<GPRobj> gprobj;
    public:
        GPRMatter(Matter initMatter, std::shared_ptr<GPRobj> gpf);
        ~GPRMatter(); // Destructor
        std::pair<double, AtomMatrix> gpr_energy_forces();
        std::pair<double, AtomMatrix> true_free_energy_forces();
        bool isCloseToTrue(double eps = 1e-6);
        void updateMatter(const Matter otherMatter);
};
