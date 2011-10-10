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
#include "Minimizer.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include "HelperFunctions.h"
#include <vector>

class LBFGS : public Minimizer 
{

public:
    LBFGS(Matter *matter, Parameters *parameters);
    LBFGS(Matter *matter, Parameters *parameters, AtomMatrix forces);

    ~LBFGS();

    void oneStep();
    long fullRelax();
    void update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0);
    bool isItConverged(double convergeCriterion);
    void setOutput(int level);

private:
    int outputLevel;

    Parameters *parameters;

    ObjectiveFunction *objf;

    double ePrev;

    int iteration;
    int memory;
    std::vector<VectorXd> s;
    std::vector<VectorXd> y;
    std::vector<double> rho;
    VectorXd rPrev;
    VectorXd fPrev;

    VectorXd getDescentDirection();

};

#endif
