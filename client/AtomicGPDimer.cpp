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
    matterCenter  = new Matter(parameters);
    matterDimer   = new Matter(parameters);
    *matterCenter = *matter;
    *matterDimer  = *matter;
    InputParameters p = eon_parameters_to_gpr(params);
    atmd::AtomicDimer atomic_dimer;
    aux::ProblemSetUp problem_setup;
}

AtomicGPDimer::~AtomicGPDimer()
{
    delete matterCenter;
    delete matterDimer;
}

void AtomicGPDimer::compute(Matter *matter, AtomMatrix initialDirection)
{
}

double AtomicGPDimer::getEigenvalue()
{
    return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector()
{
    return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
