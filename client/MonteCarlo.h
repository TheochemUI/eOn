//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

// dimer method to find the lowest curvature mode
class MonteCarlo {

public:

    MonteCarlo(Matter const *matter, Parameters *parameters);
    ~MonteCarlo();

    void run(int numSteps, double temperature, double stepSize);

private:
    Parameters *parameters;
    Matter *matter;
};

#endif
