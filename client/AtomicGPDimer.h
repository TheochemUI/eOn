#ifndef ATOMICGPDIMER_H
#define ATOMICGPDIMER_H

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"
#include "LowestEigenmode.h"
#include <vector>

#include "gprdimer/gpr/AtomicDimer.h"


#include "gprdimer/gpr/Enums.h"
#include "gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "gprdimer/gpr/covariance_functions/ConstantCF.h"
#include "gprdimer/gpr/covariance_functions/SexpatCF.h"
#include "gprdimer/gpr/ml/GaussianProcessRegression.h"
#include "gprdimer/managers/io/FileManager.h"

// dimer method to find the lowest curvature mode
class AtomicGPDimer : public LowestEigenmode
{

    public:

    // Optimization for the dimer
    static const char OPT_SCG[];
    static const char OPT_LBFGS[];

    AtomicGPDimer(Matter *matter, Parameters *parameters);
    ~AtomicGPDimer();

    void compute(Matter *matter, AtomMatrix initialDirectionMatrix);
    double getEigenvalue();
    AtomMatrix getEigenvector();

private:

    Matter *matterCenter; // center of the dimer
    Matter *matterDimer; //one configuration of the dimer
    AtomMatrix direction; // direction along the dimer
    AtomMatrix rotationalPlane; // direction normal to the plane of dimer rotation
    Parameters *parameters;

    InputParameters p;
    atmd::AtomicDimer atomic_dimer;
    aux::ProblemSetUp problem_setup;
};

#endif
