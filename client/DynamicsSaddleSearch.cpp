#include "DynamicsSaddleSearch.h"
#include "Log.h"
#include "Dynamics.h"
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

    int checkInterval = int(parameters->parrepStateCheckInterval/parameters->mdTimeStepInput);

    for (int i=0;i<parameters->mdSteps;i++) {
        dyn.oneStep();
        if (i%checkInterval == 0) {
            log("Checking, step %i\n", i);

            *product = *saddle;
            product->relax(true, false);

            if (!product->compare(reactant)) {
                log("Found new state\n");

                NudgedElasticBand neb(reactant, product, parameters);
                neb.compute();
                *saddle = *neb.image[neb.climbingImage];

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
