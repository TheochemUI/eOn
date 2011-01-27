//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LOWESTEIGENMODEINTERFACE_H
#define LOWESTEIGENMODEINTERFACE_H

#include "Matter.h"

#include "Parameters.h"

#include "Eigen.h"

/* Define the interface for the lowest eigenvalue determination algorithm */
class LowestEigenmodeInterface{

public:

    // stats information
    long totalForceCalls;
    double statsTorque;
    double statsCurvature;
    double statsAngle;
    long statsRotations;

    virtual ~LowestEigenmodeInterface() {}
    void virtual initialize(Matter const *matter, AtomMatrix displacement) = 0;
    void virtual compute(Matter const *matter) = 0;
    double virtual getEigenvalue() = 0;
    virtual Matrix<double, Eigen::Dynamic, 3> getEigenvector() = 0;
    virtual void setEigenvector(AtomMatrix const eigenvector) = 0;

};

#endif
