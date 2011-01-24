//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef IMPROVEDDIMER_H
#define IMPROVEDDIMER_H

#include <math.h>
#include <cmath>
#include <cassert>
#include "debug.h"

#include "Eigen.h"

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "LowestEigenmodeInterface.h"

// dimer method to find the lowest curvature mode
class ImprovedDimer : public LowestEigenmodeInterface
{

    public:

    // Optimization for the dimer
    enum{
        OPT_SD,
        OPT_CG,
        OPT_LBFGS
    };

    ImprovedDimer(Matter const *matter, Parameters *parameters);
    ~ImprovedDimer();

    void initialize(Matter const *matter, Matrix<double, Eigen::Dynamic, 3>);
    void compute(Matter const *matter); 
    double getEigenvalue();
    void setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector);
    Matrix<double, Eigen::Dynamic, 3>  getEigenvector();

    Parameters *parameters;

    Matter *x0;                                 // Center Image
    Matter *x1;                                 // Forward image.
    Matrix<double, Eigen::Dynamic, 3> tau;      // Dimer direction.
    Matrix<double, Eigen::Dynamic, 3> theta;    // Dimer rotation direction.
    Matrix<double, Eigen::Dynamic, 3> F_R;      // Dimer rotational force.
    double C_tau;                               // Curvature along tau.

    // parameters used for conjugate gradients
    Matrix<double, Eigen::Dynamic, 3> F_R_Old;
    Matrix<double, Eigen::Dynamic, 3> thetaOld;
    double a, b, gamma;
    bool init_cg;
};

#endif



