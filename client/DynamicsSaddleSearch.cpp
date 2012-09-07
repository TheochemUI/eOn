#include "DynamicsSaddleSearch.h"
#include "Log.h"
#include "Dynamics.h"
#include "NudgedElasticBand.h"
#include "MinModeSaddleSearch.h"

DynamicsSaddleSearch::DynamicsSaddleSearch(Matter *matterPassed, 
                                           Parameters *parametersPassed)
{
    reactant = matterPassed;
    parameters = parametersPassed;
    saddle = new Matter(parameters);
    product = new Matter(parameters);
    *saddle = *reactant;
    saddleSearch = NULL;
}

DynamicsSaddleSearch::~DynamicsSaddleSearch()
{
    delete saddle;
    delete product;
}

int DynamicsSaddleSearch::run(void)
{
    log("Starting dynamics NEB saddle search\n");

    Dynamics dyn(saddle, parameters);

    int checkInterval = parameters->parrepStateCheckInterval;
    int maxSteps = parameters->mdTime/parameters->mdTimeStepInput;

    for (int i=0;i<maxSteps;i++) {
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

                if (saddleSearch != NULL) {
                    delete saddleSearch;
                }
                saddleSearch = new MinModeSaddleSearch(saddle, mode, 
                                                       reactant->getPotentialEnergy(),
                                                       parameters);
                saddleSearch->run();
                double barrier = saddle->getPotentialEnergy()-reactant->getPotentialEnergy();
                log("Barrier of %.3f\n", barrier);
                return 0; 
            }else{
                log("Still in original state\n");
                AtomMatrix mode;
                mode = saddle->getPositions() - reactant->getPositions();
                saddleSearch = new MinModeSaddleSearch(saddle, mode, 
                                                       reactant->getPotentialEnergy(),
                                                       parameters);
            }
        }
    }
    return 1;
}
