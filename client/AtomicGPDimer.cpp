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
    AtomsConfiguration a;
    Observation o, init_middle_point;
    gpr::Coord orient_init;
    ConfInfo c,co;
    o = eon_matter_to_init_obs(matter);
    c = eon_matter_to_init_confinfo(matter);
    co = eon_matter_to_init_confinfo(matter);
    a = eon_matter_to_atmconf(matter);
    init_middle_point.clear();
    initialDirection.normalize();
    direction = initialDirection;

    atomic_dimer.initialize(p, o, init_middle_point,
                            orient_init, c, co, a);
    atomic_dimer.execute();
    return;
}

double AtomicGPDimer::getEigenvalue()
{
    return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector()
{
    return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
