#include "Matter.h"
#include "Potential.h"
#include "potentials/Morse/Morse.h"
#include "potentials/GPRPotential/GPRPotential.h"

#include "HelperFunctions.h"
#include "GPRHelpers.h"
#include "Parameters.h"

#include <memory>

class GPRMatter {
    private:
        Matter truePotMatter;
        Parameters eonp;
        gpr::AtomsConfiguration atmconf;
        gpr::GPRSetup gpr_parameters;
        gpr::GaussianProcessRegression gprfunc;
    public:
        GPRMatter(Matter initMatter);
        ~GPRMatter(); // Destructor
};
