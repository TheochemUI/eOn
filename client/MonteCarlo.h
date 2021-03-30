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
