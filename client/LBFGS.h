#ifndef LBFGS_H
#define LBFGS_H

#include "Eigen.h"
#include "Matter.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include "HelperFunctions.h"
#include <vector>

#define LBFGS_EPS 1e-30

class LBFGS : public Optimizer
{

public:
    LBFGS(ObjectiveFunction *objf, Parameters *parameters);
    ~LBFGS();

    int step(double maxMove);
    int run(int maxIterations, double maxMove);
    int update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0);
    void reset(void);

private:
    VectorXd getStep(double maxMove, VectorXd f);
    Parameters *parameters;
    ObjectiveFunction *objf;

    int iteration;
    int memory;

    std::vector<VectorXd> s;
    std::vector<VectorXd> y;
    std::vector<double> rho;

    VectorXd rPrev;
    VectorXd fPrev;
};

#endif
