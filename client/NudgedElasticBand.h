#ifndef NudgedElasticBand_H
#define NudgedElasticBand_H

#include <math.h>
#include <cmath>

#include "Eigen.h"

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

class Matter;
class Parameters;

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {

public:

    enum{
        STATUS_GOOD, //0
        STATUS_INIT, //1
        STATUS_BAD_MAX_ITERATIONS, //2
    };

    NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed);
    ~NudgedElasticBand();

    void clean(void);
    int compute(void);
    void updateForces(void);
    double convergenceForce(void);
    void findExtrema(void);
    void printImageData(bool writeToFile=false);

    int atoms;
    long images, climbingImage, numExtrema;
    Matter **image; // NEB images
    AtomMatrix **tangent;
    AtomMatrix **projectedForce;
    bool movedAfterForceCall;
    double *extremumEnergy;
    double *extremumPosition;
    double *extremumCurvature;

    long maxEnergyImage;

private:

    Parameters *parameters;

};

#endif
