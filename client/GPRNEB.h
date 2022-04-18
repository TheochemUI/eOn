#ifndef GPRNEB_H_
#define GPRNEB_H_

#include "GPRMatter.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "Optimizer.h"

class GPRNEB{
    public:
    enum{
        STATUS_GOOD, //0
        STATUS_INIT, //1
        STATUS_BAD_MAX_ITERATIONS, //2
    };
    GPRNEB(std::vector<GPRMatter> initPath, Parameters params);
    ~GPRNEB();
    void clean();
    int compute();
    void updateForces();
    double convergenceForce();
    void findExtrema();
    void printImageData(bool writeToFile = false);

    size_t natoms, nimages, totImages;
    size_t maxEnergyImage, climbingImage, numExtrema;
    bool movedAfterForceCall;
    double threshold;
    Parameters params;
    std::vector<GPRMatter> imageArray;
    std::vector<AtomMatrix> tangentArray;
    std::vector<AtomMatrix> projectedForceArray;
    std::vector<double> extremumEnergies;
    std::vector<double> extremumPositions;
    std::vector<double> extremumCurvatures;
    bool needsRetraining();
    // This returns a Matter vector for retraining
    std::vector<Matter> getCurPath();
};

#endif // GPRNEB_H_
