//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LANCZOS_H
#define LANCZOS_H

#include <math.h>
#include <cmath>
#include <cassert>
#include "debug.h"

#include "Eigen.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "LowestEigenmodeInterface.h"

// dimer method to find the lowest curvature mode
class Lanczos : public LowestEigenmodeInterface
{

    public:

    Lanczos(Matter const *matter, Parameters *parameters);
    ~Lanczos();

    void initialize(Matter const *matter, AtomMatrix displacement);
    void compute(Matter const *matter); 

    double getEigenvalue();
    AtomMatrix getEigenvector();
    void setEigenvector(AtomMatrix const eigenvector);

    Parameters *parameters;

    private:
    AtomMatrix lowestEv;
    double lowestEw;
    VectorXd r;
};

#endif



