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
    bool isItConverged(double convergeCriterion);
    void setOutput(int level);
    AtomMatrix makeInfinitesimalStepModifiedForces(AtomMatrix pos);
    AtomMatrix getNewPosModifiedForces(AtomMatrix pos, AtomMatrix forceBeforeStep, AtomMatrix forceAfterStep, double maxStep,bool saddleSearch=false);

    void setForces(AtomMatrix forces); // enables the use of modified forces

private:
    int outputLevel;

    Matter *matter;
    Parameters *parameters;

    AtomMatrix force;

    int degreesOfFreedom;
    int iteration;
    int memory;
    std::vector<VectorXd> s;
    std::vector<VectorXd> y;
    std::vector<double> rho;
    VectorXd r0;
    VectorXd f0;
    

    void initialize(Matter *matter, Parameters *parameters);

    void determineSearchDirection();

    AtomMatrix getStep(AtomMatrix forceBeforeStep, AtomMatrix forceAfterStep, double maxStep, bool saddleSearch=false);
};

#endif
