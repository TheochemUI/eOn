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
private:
    Matter *reactant;
    Matter *product;
    Matter *saddle;
    Parameters *parameters;

    VectorXi movedAtoms(const double distance);
};

#endif
