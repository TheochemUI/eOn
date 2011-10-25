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
#include "ObjectiveFunction.h"
#include "Log.h"

#include <cstdlib>

using namespace helper_functions;

const char SaddleSearch::MINMODE_DIMER[] =           "dimer";
const char SaddleSearch::MINMODE_LANCZOS[] =         "lanczos";
const char SaddleSearch::MINMODE_EXACT[] =           "exact";

class MinModeObjectiveFunction : public ObjectiveFunction
{
    public:
        MinModeObjectiveFunction(Matter *matterPassed, LowestEigenmode *minModeMethodPassed,
                                 AtomMatrix modePassed, Parameters *parametersPassed)
        {
            matter = matterPassed;
            minModeMethod = minModeMethodPassed;
            eigenvector = modePassed;
            parameters = parametersPassed;
        }
        ~MinModeObjectiveFunction(void){};

        VectorXd getGradient() 
        { 
            AtomMatrix proj;
            AtomMatrix force = matter->getForces();

            minModeMethod->compute(matter, eigenvector);
            eigenvector = minModeMethod->getEigenvector();
            double eigenvalue = minModeMethod->getEigenvalue();

            proj = (force.cwise() * eigenvector).sum() * eigenvector.normalized();

            if (0 < eigenvalue) {
                if (parameters->saddlePerpForceRatio > 0.0) {
                    // reverse force parallel to eigenvector, and reduce perpendicular force
                    double const d = parameters->saddlePerpForceRatio;
                    force = d*force - (1.+d)*proj;
                }else{
                    // follow eigenmode
                    force = -proj;
                }
            }else{
                // reversing force parallel to eigenmode
                force += -2.*proj;
            }

            VectorXd forceV = VectorXd::Map(force.data(), 3*matter->numberOfAtoms());
            return -forceV;
        }
        double getEnergy() { return matter->getPotentialEnergy(); }
        void setPositions(VectorXd x) { matter->setPositionsV(x); }
        VectorXd getPositions() { return matter->getPositionsV(); }
        int degreesOfFreedom() { return 3*matter->numberOfAtoms(); }
        bool isConverged() { return getConvergence() < parameters->optConvergedForce; }
        double getConvergence() { return matter->maxForce(); }

    private:
        AtomMatrix eigenvector;
        LowestEigenmode *minModeMethod;
        Matter *matter;
        Parameters *parameters;

};

SaddleSearch::SaddleSearch(Matter *matterPassed, AtomMatrix modePassed, 
                          double reactantEnergyPassed, Parameters *parametersPassed)
{
    reactantEnergy = reactantEnergyPassed;
    matter = matterPassed;
    mode = modePassed;
    parameters = parametersPassed;
    status = STATUS_GOOD;
    iteration = 0;

    if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        if (parameters->dimerImproved) {
            minModeMethod = new ImprovedDimer(matter, parameters);
        }else{
            minModeMethod = new Dimer(matter, parameters);
        }
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
        minModeMethod = new Lanczos(matter, parameters);
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_EXACT) {
        minModeMethod = new ExactMinMode(matter, parameters);
    }
     
}

SaddleSearch::~SaddleSearch()
{
    delete minModeMethod;
}

int SaddleSearch::run()
{
    log("Saddle point search started from reactant with energy %f eV.\n", reactantEnergy);

    if(parameters->saddleMinmodeMethod == MINMODE_DIMER) {
        log("[Dimer]  %9s   %9s   %10s   %9s   %9s   %7s   %6s   %4s\n", 
            "Step", "Step Size", "Delta E", "Force", "Curvature", 
            "Torque", "Angle", "Rots");
    }else if (parameters->saddleMinmodeMethod == MINMODE_LANCZOS) {
        log("[Lanczos]  %9s  %9s  %10s  %9s  %9s\n", 
            "Step", "Step Size", "Delta E", "Force", "Curvature");
    }

    ostringstream climb;
    climb << "climb";
    if(parameters->writeMovies)
    {
        matter->matter2con(climb.str(), false);
    }

    MinModeObjectiveFunction objf(matter, minModeMethod, mode, parameters);
    objf.getGradient();

    if (parameters->saddleNonnegativeDisplacementAbort && minModeMethod->getEigenvalue() > 0) {
        printf("%f\n", minModeMethod->getEigenvalue());
        return STATUS_NONNEGATIVE_ABORT;
    }

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);
    
    while (!objf.isConverged()) {

        if (iteration >= parameters->saddleMaxIterations) {
            status = STATUS_BAD_MAX_ITERATIONS;
            break;
        }

        AtomMatrix pos = matter->getPositions();

        optimizer->step(parameters->optMaxMove);
        iteration++;

        double de = objf.getEnergy()-reactantEnergy;
        double stepSize = (matter->pbc(matter->getPositions() - pos )).norm();

        if (de > parameters->saddleMaxEnergy) {
            status = STATUS_BAD_HIGH_ENERGY;
            break;
        }
            if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
            {
                log("[Dimer]  %9ld   %9.7f   %10.4f   %9.5f   %9.5f   %7.5f   %6.3f   %4ld\n",
                            iteration, stepSize, matter->getPotentialEnergy()-reactantEnergy,
                            matter->maxForce(),
                            minModeMethod->getEigenvalue(),
                            minModeMethod->statsTorque,
                            minModeMethod->statsAngle,
                            minModeMethod->statsRotations);
            }else if (parameters->saddleMinmodeMethod == MINMODE_LANCZOS) {
                log("[Lanczos]  %9ld  % 9.6f   %10.4f  %9.5f  %9.5f\n", 
                    iteration, stepSize, matter->getPotentialEnergy()-reactantEnergy,
                    matter->maxForce(),
                    minModeMethod->getEigenvalue());
            }

        if (parameters->writeMovies) {
            matter->matter2con(climb.str(), true);
        }

        if (parameters->checkpoint) {
            matter->matter2con("displacement_checkpoint.con", false);
            FILE *fileMode = fopen("mode_checkpoint.dat", "wb");
            helper_functions::saveMode(fileMode, matter, 
                                       minModeMethod->getEigenvector());
            fclose(fileMode);
        }
    }

    delete optimizer;

    return status;
}

double SaddleSearch::getEigenvalue()
{
    return minModeMethod->getEigenvalue();
}

AtomMatrix SaddleSearch::getEigenvector()
{
    return minModeMethod->getEigenvector();
}
