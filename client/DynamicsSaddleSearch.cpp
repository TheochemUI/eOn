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

    ifstream massFile("masses.dat");
    if (massFile.is_open()) {
        log("Found mass weights file\n");
        massFile.close();
        VectorXd masses = helper_functions::loadMasses("masses.dat", saddle->numberOfAtoms());
        saddle->setMasses(masses);
        log("Applied mass weights\n");
    }else{
        log("No mass weights file found\n");
        massFile.close();
    }

    Dynamics dyn(saddle, parameters);
    log("Initializing velocities from Maxwell-Boltzmann distribution\n");
    dyn.setTemperature(parameters->saddleDynamicsTemperature);
    dyn.setThermalVelocity();

    int dephaseSteps = int(floor(parameters->parrepDephaseTime/parameters->mdTimeStep+0.5));

    while (true) {

        log("Dephasing: %i steps\n", dephaseSteps);
        //always start from the initial configuration
        *saddle = *reactant;
        dyn.setThermalVelocity();

        // Dephase MD trajectory
        for (int step=1; step<=dephaseSteps; step++) {
            dyn.oneStep(step);
        }

        // Check to see if a transition occured
        Matter min(parameters);
        min = *saddle;
        min.relax();

        if (min.compare(reactant)) {
            log("Dephasing successful\n");
            break;
        }else{
            log("Transition occured during dephasing; Restarting\n");
            dephaseSteps /= 2;
            if (dephaseSteps < 1) dephaseSteps = 1;
        }
    }

    BondBoost bondBoost(saddle, parameters);
    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        log("Initializing Bond Boost\n");
        bondBoost.initialize();
    }

    int checkInterval = int(parameters->saddleDynamicsStateCheckInterval/parameters->mdTimeStep+0.5);
    int recordInterval = int(parameters->saddleDynamicsRecordInterval/parameters->mdTimeStep+0.5);

    if (parameters->writeMovies == true) {
        saddle->matter2con("dynamics", false);
    }

    for (int step=1;step<=parameters->mdSteps;step++) {
        dyn.oneStep(step);

        if (step % recordInterval == 0 && recordInterval != 0) {
            log("recording configuration at step %i time %.3f\n", step, step*parameters->mdTimeStep*parameters->timeUnit);
            Matter *tmp = new Matter(parameters);    
            *tmp = *saddle;
            MDSnapshots.push_back(tmp);
            MDTimes.push_back(step*parameters->mdTimeStep);
        }

        if (parameters->writeMovies == true) {
            saddle->matter2con("dynamics", true);
        }

        if (step%checkInterval == 0) {
            log("Minimizing trajectory, step %i\n", step);

            *product = *saddle;
            product->relax(false, false);

            if (!product->compare(reactant)) {
                log("Force calls total: %i\n", Potential::fcallsTotal);
                log("Found new state\n");
                int image = refineTransition(MDSnapshots, product);
                *saddle = *MDSnapshots[image];
                log("Found transition at snapshot image %i\n", image);
                for (int ii=0;ii<(int)MDTimes.size();ii++) log("MDTimes[%i] = %.3f\n", ii, MDTimes[ii]*parameters->timeUnit);
                // subtract off half the record interval in order to not introduce a systematic
                // bias towards longer times.
                time = MDTimes[image] - parameters->saddleDynamicsRecordInterval/2.0;
                log("Transition time %.2f fs\n", time*parameters->timeUnit);

                NudgedElasticBand neb(reactant, product, parameters);

                if (parameters->saddleDynamicsLinearInterpolation == false) {
                    log("Interpolating initial band through MD transition state\n");
                    AtomMatrix reactantToSaddle = saddle->pbc(saddle->getPositions()  - reactant->getPositions());
                    AtomMatrix saddleToProduct  = saddle->pbc(product->getPositions() - saddle->getPositions());
                    log("Initial band saved to neb_initial_band.con\n");
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
                }else{
                    log("Linear interpolation between minima used for initial band\n");
                    neb.image[0]->matter2con("neb_initial_band.con", false);
                    for(int j=1; j<=neb.images+1; j++){
                        neb.image[j]->matter2con("neb_initial_band.con", true);
                    }
                }

                AtomMatrix mode;
                if (parameters->nebMaxIterations > 0) {
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

                    neb.compute();
                    neb.printImageData(true);
                    int extremumImage = -1; 
                    int j;
                    for (j=0;j<neb.numExtrema;j++) {
//                        if (neb.extremumCurvature[j] < 0.0) { 
                        if (neb.extremumCurvature[j] < parameters->saddleDynamicsMaxInitCurvature) { 
                            extremumImage = (int)floor(neb.extremumPosition[j]);
                            *saddle = *neb.image[extremumImage];
                            double interpDistance = neb.extremumPosition[j] - (double)extremumImage;
                            AtomMatrix bandDirection = saddle->pbc(neb.image[extremumImage+1]->getPositions() - 
                                                                   neb.image[extremumImage]->getPositions());
                            saddle->setPositions(interpDistance * bandDirection + saddle->getPositions());
                            mode = saddle->pbc( neb.image[extremumImage+1]->getPositions() - saddle->getPositions());
                            mode.normalize();
                            minModeMethod->compute(saddle, mode);
                            double eigenvalue = minModeMethod->getEigenvalue();
                            log("extrema #%i has eigenvalue %.8f\n", j+1, eigenvalue);

                            if (eigenvalue < 0) {
                                log("chose image %i (extrema #%i) as extremum image\n", extremumImage, j+1);
                                break;
                            }else{
                                extremumImage=-1;
                            }
                        }
                    }

                    delete minModeMethod;

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
                        if (maxEnergy <= reactant->getPotentialEnergy()) {
                            log("warning: no barrier found\n");
                            return MinModeSaddleSearch::STATUS_BAD_NO_BARRIER;
                        }
                    }
                }else{
                    neb.maxEnergyImage = neb.images/2 + 1;
                }

                log("Initial saddle guess saved to saddle_initial_guess.con\n");
                saddle->matter2con("saddle_initial_guess.con");
                MinModeSaddleSearch search = MinModeSaddleSearch(saddle, 
                        mode, reactant->getPotentialEnergy(), parameters);
                int minModeStatus = search.run();

                if (minModeStatus != MinModeSaddleSearch::STATUS_GOOD) {
                    log("error in min mode saddle search\n");
                    return minModeStatus;
                }

                eigenvalue = search.getEigenvalue();
                eigenvector = search.getEigenvector();
                log("eigenvalue: %.3f\n", eigenvalue);

                double barrier = saddle->getPotentialEnergy()-reactant->getPotentialEnergy();
                log("found barrier of %.3f\n", barrier);
                for (unsigned int i=0;i<MDSnapshots.size();i++) {
                    delete MDSnapshots[i];
                }
                MDSnapshots.clear();
                MDTimes.clear();
                log("Force calls total: %i\n", Potential::fcallsTotal);
                return MinModeSaddleSearch::STATUS_GOOD; 
            }else{
                log("Still in original state\n");
                for (unsigned int i=0;i<MDSnapshots.size();i++) {
                    delete MDSnapshots[i];
                }
                MDTimes.clear();
                MDSnapshots.clear();
            }
        }
    }

    MDSnapshots.clear();
    time = parameters->mdSteps*parameters->mdTimeStep;
    return MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT;
}

int DynamicsSaddleSearch::refineTransition(std::vector<Matter*> MDSnapshots, Matter *product)
{
    int min, max, mid;
    bool midTest;
    min = 0;
    max = MDSnapshots.size() - 1;
    if (max == 0) {
        return 0;
    }

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

    return (min+max)/2;
}

double DynamicsSaddleSearch::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix DynamicsSaddleSearch::getEigenvector()
{
    return eigenvector;
}
