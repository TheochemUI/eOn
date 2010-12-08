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

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

#include "Matter.h"
#include "Parameters.h"

class Hessian
{
public:
    enum
    {
        REACTANT = 0,
        SADDLE,
        PRODUCT
    };
    Hessian(Matter *reactant, Matter *saddle, Matter *product, Parameters *params);
    ~Hessian();

    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getHessian(int which);
    VectorXd getModes(int which);

private:
    Matter *reactant;
    Matter *product;
    Matter *saddle;
    Parameters *parameters;

    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessians[3];
    VectorXd modes[3];

    VectorXi movedAtoms(const double distance);
    bool calculate(int which);
};

#endif
