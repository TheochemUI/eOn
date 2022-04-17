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
        gpr::GaussianProcessRegression gprfunc;
        GPRPotential trainedGPR;
        gpr::Observation prepobs(std::vector<Matter>& matvec, Potential* pot);
    public:
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
        Parameters eonp;
        gpr::AtomsConfiguration atmconf;
        gpr::GPRSetup gpr_parameters;
        gpr::GaussianProcessRegression gprfunc;
        std::unique_ptr<GPRPotential> trainedGPR;
        gpr::Observation curpath;
    public:
        GPRMatter(Matter initMatter);
        ~GPRMatter(); // Destructor
        void trainGPR(gpr::Observation obspath);
        void retrainGPR(gpr::Observation& newpath);
        std::pair<double, AtomMatrix> gpr_energy_forces(Matter& getat);
        std::pair<double, AtomMatrix> true_free_energy_forces(Matter& getat);
        bool isCloseTo(Matter& testat, double eps = 1e-6);
};
