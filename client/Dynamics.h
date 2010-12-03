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

/** Functionality relying on the conjugate gradients algorithm. The object minimizes a Matter object or modified forces passed in.*/
class Dynamics {

public:

    enum{
        ANDERSEN,
        NOSE_HOVER
    };

    /** Constructor to be used when a structure is minimized.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    Dynamics(Matter *matter, Parameters *parameters);

    ~Dynamics();///< Destructor.

    void oneStep(double temperature);	
    void verletStep1();
    void verletStep2();
    void fullSteps(double temperature); 
    void andersen(double temperature);
    void velocityScale(double temperature);
    void noseHoverVerlet(double temperature);


private:
    long nAtoms;///< Number of free coordinates.

    Matter *matter;///< Pointer to atom object \b outside the scope of the class.
    Parameters *parameters;///< Pointer to a structure outside the scope of the class containing runtime parameters.
  
    double dtScale;
    double kb;
    bool init;
};

#endif
