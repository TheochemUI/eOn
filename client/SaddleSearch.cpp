//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "SaddleSearch.h"
#include "ConjugateGradients.h"
#include "HelperFunctions.h"
#include "Lanczos.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "ExactMinMode.h"
#include "EpiCenters.h"
#include "Log.h"

#include <cstdlib>

using namespace helper_functions;

const char SaddleSearch::MINMODE_DIMER[] =           "dimer";
const char SaddleSearch::MINMODE_LANCZOS[] =         "lanczos";
const char SaddleSearch::MINMODE_EXACT[] =           "exact";

SaddleSearch::SaddleSearch(Matter *matterPassed, AtomMatrix modePassed, double reactantEnergyPassed, Parameters *parametersPassed)
{
    reactantEnergy = reactantEnergyPassed;
    matter = matterPassed;
    mode = modePassed;
    parameters = parametersPassed;
    status = STATUS_GOOD;
    iteration = 0;
    objf = NULL;
}

SaddleSearch::~SaddleSearch()
{
    if (objf != NULL) {
        delete objf;
    }
}

int SaddleSearch::run()
{
    log("Saddle point search started from reactant with energy %f eV.\n", reactantEnergy);

    if(parameters->saddleMinmodeMethod == MINMODE_DIMER) {
        log("[Dimer]  %9s   %9s   %16s   %9s   %9s   %9s   %9s   %9s\n", 
            "Step", "Step Size", "Energy", "Force", "Curvature", 
            "Torque", "Angle", "Rotations");
    }else if (parameters->saddleMinmodeMethod == MINMODE_LANCZOS) {
        log("[Lanczos]  %9s  %9s  %16s  %9s  %9s\n", 
            "Step", "Step Size", "Energy", "Force", "Curvature");
    }

    ostringstream climb;
    climb << "climb";
    if(parameters->writeMovies)
    {
        matter->matter2con(climb.str(), false);
    }
     
    objf = new MinModeObjectiveFunction(matter, mode, parameters);
    Optimizer *optimizer = Optimizer::getOptimizer(objf, parameters);
    
    while (!objf->isConverged()) {
        if (iteration >= parameters->saddleMaxIterations) {
            status = STATUS_BAD_MAX_ITERATIONS;
            break;
        }

        AtomMatrix pos = matter->getPositions();

        optimizer->step(parameters->optMaxMove);
        iteration++;

        double de = objf->getEnergy()-reactantEnergy;
        double stepSize = (matter->pbc(matter->getPositions() - pos )).norm();

        if (de > parameters->saddleMaxEnergy) {
            status = STATUS_BAD_HIGH_ENERGY;
            break;
        }
            if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
            {
                log("[Dimer]  %9ld  % 9.3e   %16.4f  % 9.3e  % 9.3e  % 9.3e  % 9.3e   % 9ld\n",
                            iteration, stepSize, matter->getPotentialEnergy(),
                            matter->getForces().norm(),
                            objf->minModeMethod->getEigenvalue(),
                            objf->minModeMethod->statsTorque,
                            objf->minModeMethod->statsAngle,
                            objf->minModeMethod->statsRotations);
            }else if (parameters->saddleMinmodeMethod == MINMODE_LANCZOS) {
                log("[Lanczos]  %9ld  % 9.3f   %16.4f  % 9.3f  % 9.3f\n", 
                    iteration, stepSize, matter->getPotentialEnergy(),
                    matter->getForces().norm(),
                    objf->minModeMethod->getEigenvalue());
            }

        if (parameters->writeMovies) {
            matter->matter2con(climb.str(), true);
        }

        if (parameters->checkpoint) {
            matter->matter2con("displacement_checkpoint.con", false);
            FILE *fileMode = fopen("mode_checkpoint.dat", "wb");
            helper_functions::saveMode(fileMode, matter, 
                                       objf->minModeMethod->getEigenvector());
            fclose(fileMode);
        }
    }

    delete optimizer;

    return status;
}

double SaddleSearch::getEigenvalue()
{
    return objf->eigenvalue;
}

AtomMatrix SaddleSearch::getEigenvector()
{
    return objf->eigenvector;
}
