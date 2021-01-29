// An interface to the GPDimer library

#include "HelperFunctions.h"
#include "AtomicGPDimer.h"
#include "Log.h"
#include <cmath>
#include <cassert>

#include "gprdimer/gpr/Enums.h"
#include "gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "gprdimer/gpr/AtomicDimer.h"
#include "gprdimer/gpr/covariance_functions/ConstantCF.h"
#include "gprdimer/gpr/covariance_functions/SexpatCF.h"
#include "gprdimer/gpr/ml/GaussianProcessRegression.h"

using namespace helper_functions;

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(Matter *matter, Parameters *params)
{
    parameters    = params;
    x0            = new Matter(parameters);
    *x0           = *matter;
    InputParameters p = eon_parameters_to_gpr(params);
    atmd::AtomicDimer atomic_dimer;
    aux::ProblemSetUp problem_setup;

    typedef gpr::Field<double> FieldDbl;
    FieldDbl E_all_init;    // no initial data
    gpr::Coord R_all_init;  // no initial data points
    gpr::Coord G_all_init;  // no initial data
    E_all_init.assignFromEigenMatrix(x0->getPotentialEnergy());
}

AtomicGPDimer::~AtomicGPDimer()
{
    delete x0;
}

void AtomicGPDimer::compute()
{
    atomic_dimer->execute();
}

double AtomicGPDimer::getEigenvalue()
{
    atomic_dimer->getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector()
{
   gpr::Coord *finVec = atomic_dimer->getFinalOrientation();
    // return
    return finVec->extractEigenVector();
}
