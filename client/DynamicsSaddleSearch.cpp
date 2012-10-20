#include "DynamicsSaddleSearch.h"
#include "Log.h"
#include "Dynamics.h"
#include "BondBoost.h"
#include "NudgedElasticBand.h"
#include "MinModeSaddleSearch.h"
#include "LowestEigenmode.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"

DynamicsSaddleSearch::DynamicsSaddleSearch(Matter *matterPassed, 
                                           Parameters *parametersPassed)
{
    reactant = new Matter(parameters);
    *reactant = *matterPassed;
    parameters = parametersPassed;
    product = new Matter(parameters);
    saddle = matterPassed;
}

DynamicsSaddleSearch::~DynamicsSaddleSearch()
{
    delete reactant;
    delete product;
}

int DynamicsSaddleSearch::run(void)
{
    log("Starting dynamics NEB saddle search\n");

    Dynamics dyn(saddle, parameters);
    dyn.setTemperature(parameters->saddleDynamicsTemperature);
    dyn.setThermalVelocity();

    BondBoost bondBoost(saddle, parameters);
    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        bondBoost.initialize();
    }

    int checkInterval = int(parameters->saddleDynamicsStateCheckInterval/parameters->mdTimeStepInput);

    if (parameters->writeMovies == true) {
        saddle->matter2con("dynamics", false);
    }

    for (int i=0;i<parameters->mdSteps;i++) {
        dyn.oneStep();

        if (parameters->writeMovies == true) {
            saddle->matter2con("dynamics", true);
        }

        if (i%checkInterval == 0 && i > 0) {
            log("Checking, step %i\n", i);

            *product = *saddle;
            product->relax(true, false);

            if (!product->compare(reactant)) {
                log("Found new state\n");

                NudgedElasticBand neb(reactant, product, parameters);
                neb.compute();
                *saddle = *neb.image[neb.climbingImage];

                if (saddle->maxForce() > parameters->optConvergedForce) {
                    log("CI-NEB did not converge to a saddle, force too big\n");
                    return MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY;
                }

                AtomMatrix mode;
                mode = saddle->getPositions()- neb.image[neb.climbingImage-1]->getPositions();
                mode.normalize();

                LowestEigenmode *minModeMethod;

                if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
                    if (parameters->dimerImproved) {
                        minModeMethod = new ImprovedDimer(saddle, parameters);
                    }else{
                        minModeMethod = new Dimer(saddle, parameters);
                    }
                }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
                    minModeMethod = new Lanczos(saddle, parameters);
                }

                minModeMethod->compute(saddle, mode);
                eigenvalue = minModeMethod->getEigenvalue();
                eigenvector = minModeMethod->getEigenvector();
                log("eigenvalue: %.3f\n", eigenvalue);
                delete minModeMethod;

                if (eigenvalue > 0.0) {
                    log("eigenvalue not negative\n");
                    //XXX:error is not meaningful
                    return MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY;
                }

                double barrier = saddle->getPotentialEnergy()-reactant->getPotentialEnergy();
                log("Barrier of %.3f\n", barrier);
                return MinModeSaddleSearch::STATUS_GOOD; 
            }else{
                log("Still in original state\n");
            }
        }
    }
    return MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS; 
}

double DynamicsSaddleSearch::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix DynamicsSaddleSearch::getEigenvector()
{
    return eigenvector;
}
