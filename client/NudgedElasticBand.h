#ifndef NudgedElasticBand_H
#define NudgedElasticBand_H

#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"

#include <cmath>
#include <math.h>

class Matter;
class Parameters;

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {

public:
    enum {
        STATUS_GOOD,               // 0
        STATUS_INIT,               // 1
        STATUS_BAD_MAX_ITERATIONS, // 2
    };

    NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed);
    ~NudgedElasticBand();

    void clean(void);
    int compute(void);
    void updateForces(void);
    double convergenceForce(void);
    void findExtrema(void);
    void printImageData(bool writeToFile = false);

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

class NEBObjectiveFunction : public ObjectiveFunction {
public:
    NEBObjectiveFunction(NudgedElasticBand *nebPassed, Parameters *parametersPassed)
        : neb{nebPassed}, parameters{parametersPassed} {}
    ~NEBObjectiveFunction(){};
    double getEnergy();
    double getConvergence();
    bool isConverged();
    int degreesOfFreedom();
    VectorXd getPositions();
    VectorXd getGradient(bool fdstep = false);
    VectorXd difference(VectorXd a, VectorXd b);
    void setPositions(VectorXd x);

private:
    NudgedElasticBand *neb;
    Parameters *parameters;
};
#endif
