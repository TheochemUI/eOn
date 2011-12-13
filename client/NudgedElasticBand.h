//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

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

    int atoms;
    long images;
    long climbingImage;
    Matter **image; // NEB images
    AtomMatrix **tangent;

private:

    Parameters *parameters;

};

#endif
