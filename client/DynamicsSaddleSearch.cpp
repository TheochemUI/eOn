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
    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
}

DynamicsSaddleSearch::~DynamicsSaddleSearch()
{
    delete reactant;
    delete product;
}

int DynamicsSaddleSearch::run(void)
{
    std::vector<Matter*> MDSnapshots;
    log("Starting dynamics NEB saddle search\n");

    Dynamics dyn(saddle, parameters);
    dyn.setTemperature(parameters->saddleDynamicsTemperature);
    dyn.setThermalVelocity();

    BondBoost bondBoost(saddle, parameters);
    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        bondBoost.initialize();
    }

    int checkInterval = int(parameters->saddleDynamicsStateCheckInterval/parameters->mdTimeStepInput);
    int recordInterval = int(parameters->saddleDynamicsRecordInterval/parameters->mdTimeStepInput);

    if (parameters->writeMovies == true) {
        saddle->matter2con("dynamics", false);
    }

    for (int i=0;i<parameters->mdSteps;i++) {
        dyn.oneStep();

        if (i % recordInterval == 0) {
            Matter *tmp = new Matter(parameters);    
            *tmp = *saddle;
            MDSnapshots.push_back(tmp);
        }

        //if (parameters->writeMovies == true) {
        //    saddle->matter2con("dynamics", true);
        //}

        if (i%checkInterval == 0 && i > 0) {
            log("Checking, step %i\n", i);

            *product = *saddle;
            product->relax(false, false);

            if (!product->compare(reactant)) {
                log("Found new state\n");
                int image = refineTransition(MDSnapshots, product);
                *saddle = *MDSnapshots[image];
                printf("Found trasition at image %i\n", image);
                printf("Time %.2f fs\n", image * parameters->saddleDynamicsRecordInterval + 
                        parameters->mdTimeStepInput*(i%checkInterval));

                product->relax(false, false);

                NudgedElasticBand neb(reactant, product, parameters);

                if (parameters->saddleDynamicsLinearInterpolation == false) {
                    AtomMatrix reactantToSaddle = saddle->pbc(saddle->getPositions()  - reactant->getPositions());
                    AtomMatrix saddleToProduct  = saddle->pbc(product->getPositions() - saddle->getPositions());
                    neb.image[0]->matter2con("neb_initial_band.con", false);
                    for (int image=1;image<=neb.images;image++) {
                        int mid = neb.images/2 + 1;
                        if (image < mid) {
                            double frac = ((double)image) / ((double)mid);
                            neb.image[image]->setPositions(reactant->getPositions() + frac * reactantToSaddle);
                        }else if (image > mid) {
                            double frac = (double)(image-mid) / (double)(neb.images - mid + 1);
                            neb.image[image]->setPositions(saddle->getPositions() + frac * saddleToProduct);
                        }else if (image == mid) {
                            neb.image[image]->setPositions(saddle->getPositions());
                        }
                        neb.image[image]->matter2con("neb_initial_band.con",true);
                    }
                    neb.image[neb.images+1]->matter2con("neb_initial_band.con",true);
                }

                if (parameters->nebMaxIterations > 0) {
                    neb.compute();
                    neb.printImageData(true);
                    *saddle = *neb.image[neb.maxEnergyImage];
                }else{
                    neb.maxEnergyImage = neb.images/2 + 1;
                }

                AtomMatrix mode;
                mode = saddle->getPositions()- neb.image[neb.maxEnergyImage-1]->getPositions();
                mode.normalize();


                MinModeSaddleSearch search = MinModeSaddleSearch(saddle, 
                        mode, reactant->getPotentialEnergy(), parameters);
                search.run();

                if (saddle->maxForce() > parameters->optConvergedForce) {
                    log("did not converge to a saddle, force too big\n");
                    return MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY;
                }


                eigenvalue = search.getEigenvalue();
                eigenvector = search.getEigenvector();
                log("eigenvalue: %.3f\n", eigenvalue);

                if (eigenvalue > 0.0) {
                    log("eigenvalue not negative\n");
                    //XXX:error is not meaningful
                    return MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY;
                }

                double barrier = saddle->getPotentialEnergy()-reactant->getPotentialEnergy();
                log("found barrier of %.3f\n", barrier);
                MDSnapshots.clear();
                return MinModeSaddleSearch::STATUS_GOOD; 
            }else{
                log("Still in original state\n");
                MDSnapshots.clear();
            }
        }
    }
    MDSnapshots.clear();
    return MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS; 
}

int DynamicsSaddleSearch::refineTransition(std::vector<Matter*> MDSnapshots, Matter *product)
{
    int min, max, mid;
    bool midTest;
    min = 0;
    max = MDSnapshots.size() - 1;

    while( (max-min) > 1 ) {
        mid = min + (max-min)/2;
        Matter snapshot(parameters);
        snapshot = *MDSnapshots[mid];

        snapshot.relax(true);

        midTest = snapshot.compare(reactant);

        if (midTest){
            min = mid;
        } else {
            *product = snapshot;
            max = mid;
        }
    }

    return (min+max)/2 + 1;
}

double DynamicsSaddleSearch::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix DynamicsSaddleSearch::getEigenvector()
{
    return eigenvector;
}
