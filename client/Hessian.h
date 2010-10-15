#ifndef HESSIAN_H
#define HESSIAN_H

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

class Hessian
{
public:
    enum
    {
        REACTANT = 0,
        SADDLE,
        PRODUCT
    };
    Hessian(const Matter* reactant, const Matter* saddle, const Matter* product, Parameters* params);
    ~Hessian();

    double getModeProduct(int which);
private:
    Matter *reactant;
    Matter *product;
    Matter *saddle;

    VectorXi movedAtoms(const double distance) const;
}

#endif
