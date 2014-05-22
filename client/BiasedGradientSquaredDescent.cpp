//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "BiasedGradientSquaredDescent.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <string.h>
#include <map>

static const char LOG_PREFIX[] = "[BGSD]";

class BGSDObjectiveFunction : public ObjectiveFunction
{
    public:
        BGSDObjectiveFunction(Matter *matterPassed,
                              Parameters *parametersPassed)
        {
            matter = matterPassed;
            parameters = parametersPassed;
        }
        ~BGSDObjectiveFunction(void){};
        double getEnergy()
        { 

            double energy = matter->getPotentialEnergy();
            double energySq = energy * energy;
            return energySq;
        }
        VectorXd getGradient(bool fdstep=false) { return -matter->getForcesFreeV(); }
        void setPositions(VectorXd x) { matter->setPositionsFreeV(x); }
        VectorXd getPositions() { return matter->getPositionsFreeV(); }
        int degreesOfFreedom() { return 3*matter->numberOfFreeAtoms(); }
        bool isConverged() { return getConvergence() < parameters->optConvergedForce; }

        double getConvergence() {
            if (parameters->optConvergenceMetric == "norm") {
                return matter->getForcesFreeV().norm(); 
            } else if (parameters->optConvergenceMetric == "max_atom") {
                return matter->maxForce(); 
            } else if (parameters->optConvergenceMetric == "max_component") {
                return matter->getForces().maxCoeff(); 
            } else {
                log("%s Unknown opt_convergence_metric: %s\n", LOG_PREFIX,
                    parameters->optConvergenceMetric.c_str());
                exit(1);
            }
        }

        VectorXd difference(VectorXd a, VectorXd b) {
            return matter->pbcV(a-b);
        }
    private:
        Matter *matter;
        Parameters *parameters;
};



BiasedGradientSquaredDescent::BiasedGradientSquaredDescent(Matter *matterPassed, double reactantEnergyPassed, Parameters *parametersPassed)
{
    parameters = parametersPassed;
    reactantEnergy = reactantEnergyPassed;
    saddle = matterPassed;
    eigenvector.resize(saddle->numberOfAtoms(), 3);
    eigenvector.setZero();
}

BiasedGradientSquaredDescent::~BiasedGradientSquaredDescent()
{
}

int BiasedGradientSquaredDescent::run()
{
    BGSDObjectiveFunction objf(saddle, parameters);
    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);
    int iteration = 0;
    while (!objf.isConverged() || iteration == 0) {
        optimizer->step(parameters->optMaxMove);
        printf("iteration %i energy: %.8f\n", iteration, saddle->getPotentialEnergy());
        iteration++;
    }

    return 0;
}

double BiasedGradientSquaredDescent::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix BiasedGradientSquaredDescent::getEigenvector()
{
    return eigenvector;
}
