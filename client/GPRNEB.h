#ifndef GPRNEB_H_
#define GPRNEB_H_

#include "GPRMatter.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "Optimizer.h"
#include <fmt/core.h>
#include <fmt/format.h>

class GPRNEB{
    private:
        void getTrueNEBForces();
        void updateForces();
        friend class GPRNEBObjectiveFunction;
    public:
    enum{
        STATUS_GOOD, //0
        STATUS_INIT, //1
        STATUS_BAD_MAX_ITERATIONS, //2
    };
    GPRNEB(std::vector<GPRMatter> initPath, Parameters params);
    ~GPRNEB();
    void clean();
    int compute(const std::vector<Matter> ppoints, bool& isWithin);
    double convergenceForce();
    // Tuple representing convergence values for the saddlepoint and the rest of the path
    std::pair<double, double> getConvergenceTrue();
    void findExtrema();
    void printImageData(bool writeToFile = false, size_t neb_id = 0);
    void findExtremaTrue();
    void printImageDataTrue(bool writeToFile = false, size_t neb_id = 0);

    size_t natoms, nimages, totImages, nfree;
    size_t maxEnergyImage, climbingImage, numExtrema;
    bool movedAfterForceCall;
    double threshold, init_path_length;
    Parameters params;
    std::vector<std::reference_wrapper<GPRMatter> > nebImages;
    std::vector<GPRMatter> imageArray;
    std::vector<AtomMatrix> tangentArray;
    std::vector<AtomMatrix> projectedForceArray;
    std::vector<AtomMatrix> projectedForceArrayTrue;
    std::vector<AtomMatrix> tangentArrayTrue;
    std::vector<double> extremumEnergies;
    std::vector<double> extremumPositions;
    std::vector<double> extremumCurvatures;
    bool needsRetraining(double eps = 1e-3);
    bool stoppedEarly(std::vector<Matter> prevPath, double max_dist_factor = 0.5);
    bool notStoppedEarly(std::vector<Matter> prevPath, double max_dist_factor = 0.5);
    // This returns a Matter vector for retraining
    std::vector<Matter> getCurPath();
    std::vector<Matter> getCurPathFull();
    // These are outer loop helpers
    double getTrueConvForce();
};

#endif // GPRNEB_H_
