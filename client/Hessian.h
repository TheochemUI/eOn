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

    double getModeProduct(int which);
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getHessian(int which);
    VectorXd getModes(int which);

private:
    Matter *reactant;
    Matter *product;
    Matter *saddle;
    Parameters *parameters;

    double modeProducts[3];
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessians[3];
    VectorXd modes[3];

    VectorXi movedAtoms(const double distance);
    void calculate(int which);
};

#endif
