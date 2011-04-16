//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef NEB_H
#define NEB_H

#include <math.h>
#include <cmath>
#include <cassert>

#include "Eigen.h"

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

// NEB method for determining a minimum energy path between two matter objects
class NEB {

public:

    // Optimization for the neb
//    static const string OPT_QM;
//    static const string OPT_CG;
//    static const string OPT_LBFGS;
    static const char OPT_QM[];
    static const char OPT_CG[];
    static const char OPT_LBFGS[];

    NEB(Matter const *matterInitial, Matter const *matterFinal, Parameters *parameters);
    ~NEB();

    void compute(void);
    void updateForces(void);
    double convergenceForce(void);

private:

    Parameters *parameters;
    Matter *neb[]; // NEB images
    int nAtoms;
    long images;
    long climbingImage;
    AtomMatrix *tangent;
    

};

#endif
