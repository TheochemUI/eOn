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
    std::vector<double> MDTimes;
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

        if ((i+1) % recordInterval == 0) {
            Matter *tmp = new Matter(parameters);    
            *tmp = *saddle;
            MDSnapshots.push_back(tmp);
            MDTimes.push_back((i+1)*parameters->mdTimeStepInput);
        }

        //if (parameters->writeMovies == true) {
        //    saddle->matter2con("dynamics", true);
        //}

        if ((i+1)%checkInterval == 0) {
            log("Checking, step %i\n", i+1);

            *product = *saddle;
            product->relax(false, false);

            if (!product->compare(reactant)) {
                log("Force calls total: %i\n", Potential::fcallsTotal);
                log("Found new state\n");
                int image = refineTransition(MDSnapshots, product);
                *saddle = *MDSnapshots[image];
                log("Found trasition at image %i\n", image);
                time = MDTimes[image];
                log("Transition Time %.2f fs\n", time);

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

                AtomMatrix mode;
                if (parameters->nebMaxIterations > 0) {
                    neb.compute();
                    neb.printImageData(true);
                    int extremumImage = -1; 
                    int j=0;
                    for (j=0;j<neb.numExtrema;j++) {
                        if (neb.extremumCurvature[j] < 0.0) { 
                            extremumImage = floor(neb.extremumPosition[j]);
                            log("chose image %i as extremum image\n", extremumImage);
                            break;
                        }
                    }
                    if (extremumImage != -1) {
                        *saddle = *neb.image[extremumImage];
                        double interpDistance = neb.extremumPosition[j] - (double)extremumImage;
                        log("interpDistance %f\n", interpDistance);
                        AtomMatrix bandDirection = saddle->pbc(neb.image[extremumImage+1]->getPositions() - 
                                                               neb.image[extremumImage]->getPositions());
                        saddle->setPositions(interpDistance * bandDirection + saddle->getPositions());
                        mode = saddle->pbc( neb.image[extremumImage+1]->getPositions() - saddle->getPositions());
                        mode.normalize();
                    }else{
                        log("no maxima found, using max energy non-endpoint image\n");
                        double maxEnergy = -INFINITY;
                        for (int image=1;image<=neb.images;image++) {
                            double U = neb.image[image]->getPotentialEnergy();
                            if (U > maxEnergy) {
                                maxEnergy = U;
                                *saddle = *neb.image[image];
                                mode = saddle->pbc(neb.image[image+1]->getPositions() - saddle->getPositions());
                                mode.normalize();
                            }
                        }
                    }
                }else{
                    neb.maxEnergyImage = neb.images/2 + 1;
                }



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
                MDTimes.clear();
                log("Force calls total: %i\n", Potential::fcallsTotal);
                return MinModeSaddleSearch::STATUS_GOOD; 
            }else{
                log("Still in original state\n");
                MDTimes.clear();
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

    log("refining transition time\n");

    while( (max-min) > 1 ) {
        mid = min + (max-min)/2;
        log("minimizing image %i\n", mid);
        Matter snapshot(parameters);
        snapshot = *MDSnapshots[mid];

        snapshot.relax(false);

        midTest = snapshot.compare(reactant);

        if (midTest){
            log("image %i minimizes to reactant\n", mid);
            min = mid;
        }else{
            log("image %i minimizes to product\n", mid);
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
