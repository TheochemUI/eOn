//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef DIMER_H
#define DIMER_H

#include <math.h>
#include <cmath>
#include <cassert>
#include "debug.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "LowestEigenmodeInterface.h"

// dimer method to find the lowest curvature mode
class Dimer : public LowestEigenmodeInterface{

public:

    Dimer(Matter const *matter, Parameters *parameters);
    ~Dimer();
        
    void initialize(Matter const *matter, Matrix<double, Eigen::Dynamic, 3>); // initialize the dimer
    void compute(Matter const *matter); // compute the lowest eigenmode
    double getEigenvalue(); // return the current eigenvalue
    void setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector); // set the dimer direction (eigenvector)
    Matrix<double, Eigen::Dynamic, 3>  getEigenvector();  // return the current eigenvector

    // stats information
    long totalForceCalls;
    double statsTorque;
    double statsCurvature;
    double statsAngle;
    double statsRotations;

private:

    Matter *matterCenter; // center of the dimer
    Matter *matterDimer; //one configuration of the dimer
    Matrix<double, Eigen::Dynamic, 3> direction; // direction along the dimer
    Matrix<double, Eigen::Dynamic, 3> rotationalPlane; // direction normal to the plane of dimer rotation 
    double eigenvalue; // current curvature along the dimer
    int nAtoms;
    Parameters *parameters;

    // The rotational plane that is going to be used is determined with the conjugate gradient method
    void determineRotationalPlane(Matrix<double, Eigen::Dynamic, 3> rotationalForce, 
                                  Matrix<double, Eigen::Dynamic, 3> &rotationalForceOld, 
                                  Matrix<double, Eigen::Dynamic, 3> rotationalPlaneNormOld,
                                  double *lengthRotationalForceOld);
    
    void rotate(double rotationAngle); // rotate the dimer by rotationAngle (radians)
    double calcRotationalForce(Matrix<double, Eigen::Dynamic, 3> &forceDiffOrthogonalToDimer); // determine the rotational force on the dimer

};

#endif
