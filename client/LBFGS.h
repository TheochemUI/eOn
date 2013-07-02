//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LBFGS_H
#define LBFGS_H

#include "Eigen.h"
#include "Matter.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include "HelperFunctions.h"
#include <vector>

class LBFGS : public Optimizer
{

public:
    LBFGS(ObjectiveFunction *objf, Parameters *parameters);
    ~LBFGS();

    bool step(double maxMove);
    bool run(int maxIterations, double maxMove);
    void update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0);
    void reset(void);

private:
    VectorXd getStep(double maxMove, VectorXd f);
    Parameters *parameters;
    ObjectiveFunction *objf;

    int iteration;
    int memory;

    std::vector<VectorXd> s;
    std::vector<VectorXd> y;
    std::vector<double> rho;

    VectorXd rPrev;
    VectorXd fPrev;
};

#endif
