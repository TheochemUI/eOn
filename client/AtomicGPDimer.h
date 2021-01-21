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
#include "gprdimer/gpr/auxiliary/Setup.h"
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

    void compute(Matter *matter, AtomMatrix initialDirection);
    double getEigenvalue();
    AtomMatrix getEigenvector();

    Parameters *parameters;

    atmd::AtomicDimer *atomic_dimer;

    Matter *x0;          // Center image
    Matter *x1;          // Forward image
    VectorXd tau;      // Dimer direction
    VectorXd theta;    // Dimer rotation direction
    VectorXd F_R;      // Dimer rotational force
    double C_tau;        // Curvature along tau

    // parameters used for conjugate gradients
    VectorXd F_R_Old;
    VectorXd thetaOld;
    double a, b, gamma;
    bool init_cg;

    // variables for LBFGS
    std::vector<VectorXd> s,y;
    std::vector<double> rho;
    bool init_lbfgs;
    VectorXd rPrev;

    std::vector<VectorXd> gradients;
    std::vector<VectorXd> positions;
};

#endif
