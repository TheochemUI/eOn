#ifndef CG_H
#define CG_H

#include "Eigen.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"

class ConjugateGradients : public Optimizer {
public:
    ConjugateGradients(ObjectiveFunction *objf, Parameters *parameters);
    ~ConjugateGradients();

    int step(double maxMove);
    int run(int maxIterations, double maxMove);
    VectorXd getStep();

private:
    ObjectiveFunction *objf;
    Parameters *parameters;

    VectorXd direction;
    VectorXd directionOld;
    VectorXd directionNorm;
    VectorXd force;
    VectorXd forceOld;

    int cg_i;
    int single_step(double maxMove);
    int line_search(double maxMove);
};

#endif
