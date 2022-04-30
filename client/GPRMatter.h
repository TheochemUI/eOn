#include "Matter.h"
#include "Potential.h"
#include "potentials/Morse/Morse.h"
#include "potentials/GPRPotential/GPRPotential.h"

#include "HelperFunctions.h"
#include "GPRHelpers.h"
#include "Parameters.h"

#include <memory>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/os.h>

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
        bool isForceMaxElementLow(double eps = 1e-3);
        void updateMatter(Matter& otherMatter);
};

namespace helper_functions {
        bool maybeUpdateGPRobj(NudgedElasticBand& neb, std::shared_ptr<GPRobj>& gpf);
        /**
         * \brief Compares energies of GPRMatter objects
         *
         * This is one of a few comparison functions to use with <algorithm>
         */
        bool max_gpr_energy(GPRMatter& matOne, GPRMatter& matTwo);
        bool max_true_energy(GPRMatter& matOne, GPRMatter& matTwo);
        /**
         * \brief Compares energy of GPRMatter to a precision
         *
         * This is one of a few comparison functions to use with <algorithm>
         */
        bool true_force_norm_converged(GPRMatter& mat, double eps = 1e-3);
        /**
         * \brief Constructs inputs for a GP-NEB round
         */
        std::vector<GPRMatter> prepGPRMatterVec(std::vector<Matter>& newPath, std::shared_ptr<GPRobj>& gpf);

        /**
         * \brief Helper to use as a comparator
         */
        template <typename T>
        bool abs_compare(const T& a, const T& b){
                return std::abs(a) < std::abs(b);
        }

        /**
         * \brief Drag a ghost atom accross a level surface, getting a slice
         *
         * For a 1 Pt on a heptamer (frozen) surface:
         * xrange -> {7.545, 11.660}
         * yrange -> {7.523, 13.074}
         * elemXY -> {50, 50}
         * zlvl -> {14.5845}
         *
         * Also, this function is meant to be used **with** a GPR object.
         */
        void peSliceSurface(const std::pair<double, double> xrange,
                            const std::pair<double, double> yrange,
                            const double zlvl,
                            const std::pair<size_t, size_t> elemXY,
                            const Matter baseObject,
                            const std::shared_ptr<GPRobj> gpf,
                            const size_t gpr_num);
        /**
         * \brief Generates a vector of Matter objects
         */
         std::vector<Matter> getMatterVector(const std::vector<GPRMatter>& gpmatvec);

} // namespace helper_functions
