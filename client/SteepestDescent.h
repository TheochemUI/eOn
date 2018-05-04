
#ifndef SteepestDescent_H
#define SteepestDescent_H

#include "Eigen.h"
#include "Matter.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include "HelperFunctions.h"
#include <vector>

class SteepestDescent : public Optimizer
{

public:
    SteepestDescent(ObjectiveFunction *objf, Parameters *parameters);
    ~SteepestDescent();

    int step(double maxMove);
    int run(int maxIterations, double maxMove);

private:
    VectorXd getStep(VectorXd f);
    Parameters *parameters;
    ObjectiveFunction *objf;

    int iteration;

    VectorXd rPrev;
    VectorXd fPrev;
};

#endif
