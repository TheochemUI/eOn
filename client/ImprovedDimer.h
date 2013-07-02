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

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"
#include "LowestEigenmode.h"
#include <vector>

// dimer method to find the lowest curvature mode
class ImprovedDimer : public LowestEigenmode
{

    public:

    // Optimization for the dimer
//    static const string OPT_SD;
//    static const string OPT_CG;
//    static const string OPT_LBFGS;
    static const char OPT_SD[];
    static const char OPT_CG[];
    static const char OPT_LBFGS[];

    ImprovedDimer(Matter *matter, Parameters *parameters);
    ~ImprovedDimer();

    void compute(Matter *matter, AtomMatrix initialDirection);
    double getEigenvalue();
    AtomMatrix getEigenvector();

    Parameters *parameters;

    Matter *x0;          // Center image
    Matter *x1;          // Forward image
    VectorXd tau;      // Dimer direction
    VectorXd theta;    // Dimer rotation direction
    VectorXd F_R;      // Dimer rotational force
    double C_tau;        // Curvature along tau

    // parameters used for conjugate gradients
    VectorXd F_R_Old;
    VectorXd thetaOld;
    double a, b, gamma;
    bool init_cg;

    // variables for LBFGS
    std::vector<VectorXd> s,y;
    std::vector<double> rho;
    bool init_lbfgs;
    VectorXd rPrev;

    std::vector<VectorXd> gradients;
    std::vector<VectorXd> positions;
};

#endif

