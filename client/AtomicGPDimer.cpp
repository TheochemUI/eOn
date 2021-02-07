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
    AtomsConfiguration atoms_config;
    Observation o, init_middle_point;
    gpr::Coord orient_init, R_init;
    aux::ProblemSetUp problem_setup;
    o = eon_matter_to_init_obs(matter);
    atoms_config = eon_matter_to_atmconf(matter);
    init_middle_point.clear();
    initialDirection.normalize();
    direction = initialDirection;
    Potential *potential = Potential::getPotential(parameters);
    // FIXME: Initialize R_init correctly
    // R_init.resize(1, R_sp.getNj());
    // double dist_sp = parameters.dist_sp.value[parameters.i_dist.value];
    // for(Index_t n = 0; n < R_sp.getNj(); ++n) {
    //     R_init[n] = R_sp[n] + dist_sp * orient_start(parameters.i_run.value, n);
    // }
    problem_setup.activateFrozenAtoms(R_init, parameters->gprActiveRadius,
                                      atoms_config);
    atomic_dimer.initialize(p, o, init_middle_point,
                            orient_init, atoms_config);
    atomic_dimer.execute(potential);
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
