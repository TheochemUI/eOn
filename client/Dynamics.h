//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "MinimizersInterface.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Constants.h"
#include "Parameters.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

class Dynamics {

public:

    enum{
        ANDERSEN,
        NOSE_HOVER
    };

    Dynamics(Matter *matter, Parameters *parameters);

    ~Dynamics();

    void oneStep(double temperature);	
    void andersenVerlet();
    void fullSteps(double temperature); 
    void andersen(double temperature);
    void velocityScale(double temperature);
    void noseHoverVerlet(double temperature);


private:
    long nAtoms;

    Matter *matter;
    Parameters *parameters;
  
    double dtScale;
    double kb;
    bool init;
};

#endif
