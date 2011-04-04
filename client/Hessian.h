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
    static const string REACTANT;
    static const string SADDLE;
    static const string PRODUCT;

    Hessian(Matter *reactant, Matter *saddle, Matter *product, Parameters *params);
    ~Hessian();

    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getHessian(string which);
    VectorXd getModes(string which);

private:
    Matter *reactant;
    Matter *product;
    Matter *saddle;
    Parameters *parameters;

    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessians[3];
    VectorXd modes[3];

    VectorXi movedAtoms(const double distance);
    bool calculate(string which);
    int whichNum(string which);
};

#endif
