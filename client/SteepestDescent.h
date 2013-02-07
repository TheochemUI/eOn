//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef SteepestDescent_H
#define SteepestDescent_H

#include "Eigen.h"
#include "Matter.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include "HelperFunctions.h"
#include <vector>

class SteepestDescent : public Optimizer
{

public:
    SteepestDescent(ObjectiveFunction *objf, Parameters *parameters);
    ~SteepestDescent();

    bool step(double maxMove);
    bool run(int maxIterations, double maxMove);

private:
    VectorXd getStep(VectorXd f);
    Parameters *parameters;
    ObjectiveFunction *objf;

    int iteration;

    VectorXd rPrev;
    VectorXd fPrev;
};

#endif
