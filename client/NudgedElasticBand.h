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
#include <cassert>

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
        STATUS_GOOD,
        STATUS_INIT
    };

    // Optimization for the neb
    static const char OPT_QM[];
    static const char OPT_CG[];
    static const char OPT_LBFGS[];

    NudgedElasticBand();
    NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed);
    ~NudgedElasticBand();

    void clean(void);
    void initialize(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed);
    int compute(void);
    void updateForces(void);
    double convergenceForce(void);

private:

    Matter *initial;
    Matter *final;
    Parameters *parameters;
    Matter *neb[]; // NEB images
    int nAtoms;
    long images;
    long climbingImage;
    AtomMatrix *tangent;

};

#endif
