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
    // Extracted from the execute function of the AtomicDimer class
    Observation middle_point;   // Observation at the middle point: coordinates,
                                // energy and gradient

    ::gpr::Field<double> E_gp;    // Vector gathering approximated energy of the
                                // middle point of the dimer for each inner iteration
    ::gpr::Coord orient;          // Unit vector along the direction of the dimer
    ::gpr::Coord orient_old;      // Unit vector along the direction of the dimer
                                // form the previous iteration

    Observation middle_point_init = genObsGPR(matterCenter);
    ::gpr::Coord orient_init = genCoordGPR(matterDimer);

    orient_init.normalizeRow(0);
    middle_point.R = middle_point_init.R;
    orient = orient_init;

    // 1) Evaluate accurate energy E0 and force F0 at the middle point R0

    if (middle_point_init.E.isEmpty()) {
        // Calculate energy and gradient at the middle point of the dimer
        atomic_dimer.callGeneraPotentialFromEON(atomic_dimer.atom_config, atomic_dimer.cell_dimensions, middle_point, Potential::getPotential(parameters));

        // Set zero level of biased potential to the energy of the middle point of
        // the initial dimer
        atomic_dimer.E_zero_level = middle_point.E;
        middle_point.E.setZero();
        problem_setup.cutOffEnergy(atomic_dimer.E_zero_level, atomic_dimer.all_obs.E);

        // Assemble coordinates, energy and gradients of all Observation& points
        atomic_dimer.all_obs.append(middle_point);
    }
    else {
        atomic_dimer.E_zero_level = middle_point_init.E;

        atomic_dimer.callGeneraPotentialFromEON(atomic_dimer.atom_config, atomic_dimer.cell_dimensions, middle_point, Potential::getPotential(parameters));

        middle_point.E = middle_point_init.E;
        middle_point.G = middle_point_init.G;
        problem_setup.cutOffEnergy(atomic_dimer.E_zero_level, middle_point.E);

        atomic_dimer.all_obs.R = middle_point_init.R;
        atomic_dimer.all_obs.E = middle_point_init.E;
        atomic_dimer.all_obs.G = middle_point_init.G;
        problem_setup.cutOffEnergy(atomic_dimer.E_zero_level, atomic_dimer.all_obs.E);
    }

    // 2) Check final convergence using accurate force F0
    if (isFinalConvergenceReached(middle_point, atomic_dimer.stop_citeria_dimer.force)) {
        log("Final convergence obtained in the beginning %9s evaluations\n",atomic_dimer.getNumEvaluations());
        return;
    }

    if (atomic_dimer.max_iter_init_rot.outer > 0) {
        // 3) Evaluate energy E1 and force F1 at R1
        updateLocation(orient, middle_point.R, atomic_dimer.image1.R);
        evaluateAccurateEnergyAndForce(atomic_dimer.image1);

        // 5) Repeat initial rotations until rotational convergence
        atomic_dimer.performInitialRotations(middle_point, orient);
    }

    if (atomic_dimer.all_obs.R.getNi() < 2) {
        updateLocation(orient, middle_point.R, atomic_dimer.image1.R);
        evaluateAccurateEnergyAndForce(atomic_dimer.image1);
        log("Evaluated image 1 of the dimer for the initial GP model.\n");
    }

    // 6) Repeat GPR iterations until final convergence
    auto eigen_data = atomic_dimer.performGPRIterations(middle_point, orient);

    atomic_dimer.middle_point_final.R = middle_point.R;
    atomic_dimer.middle_point_final.E = middle_point.E;
    atomic_dimer.middle_point_final.G = middle_point.G;
    atomic_dimer.curvature_final = eigen_data.first;
    atomic_dimer.orient_final = eigen_data.second;
}

double AtomicGPDimer::getEigenvalue()
{
    return atomic_dimer.curvature_final;
}

AtomMatrix AtomicGPDimer::getEigenvector()
{
    return atomic_dimer.orient_final.extractEigenMatrix();
}
