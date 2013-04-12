//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef HESSIAN_H
#define HESSIAN_H

#include "Eigen.h"

#include "Matter.h"
#include "Parameters.h"

class Hessian
{
public:
    Hessian(Parameters *params, Matter *matter);
    ~Hessian();

    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getHessian(Matter *matterIn, VectorXi atomsIn);
    VectorXd getFreqs(Matter *matterIn, VectorXi atomsIn);
//    VectorXd getModes(Matter *matterIn, VectorXi atomsIn);
    VectorXd removeZeroFreqs(VectorXd freqs);

private:
    Matter *matter;
    Parameters *parameters;

    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian;
//    VectorXd modes;
    VectorXd freqs;

    VectorXi atoms;
    bool calculate(void);
};

#endif
