#include "NudgedElasticBandJob.h"
#include "ConjugateGradients.h"
#include "Potential.h"
#include "Log.h"

#include <stdio.h>
#include <string>

using namespace std;

NudgedElasticBandJob::NudgedElasticBandJob(Parameters *parametersPassed)
{
    parameters = parametersPassed;
    fCallsNEB = 0;
}

NudgedElasticBandJob::~NudgedElasticBandJob()
{}

std::vector<std::string> NudgedElasticBandJob::run(void)
{
    long status;
    int f1;

    string reactantFilename = helper_functions::getRelevantFile("reactant.con");
    string productFilename = helper_functions::getRelevantFile("product.con");

    string transitionStateFilename = helper_functions::getRelevantFile("ts.con");
    Matter *transitionState = NULL;
    bool tsInterpolate = false;
    FILE *fhTransitionState = fopen("ts.con", "r");
    if (fhTransitionState != NULL) {
        tsInterpolate = true;
        fclose(fhTransitionState);
        transitionState = new Matter(parameters);
        transitionState->con2matter(transitionStateFilename);
    }

    Matter *initial = new Matter(parameters);
    Matter *final = new Matter(parameters);

    initial->con2matter(reactantFilename);
    final->con2matter(productFilename);

    NudgedElasticBand *neb = new NudgedElasticBand(initial, final, parameters);

    if (tsInterpolate) {
        AtomMatrix reactantToTS = transitionState->pbc(transitionState->getPositions()  - initial->getPositions());
        AtomMatrix TSToProduct  = transitionState->pbc(final->getPositions() - transitionState->getPositions());
        for (int image=1;image<=neb->images;image++) {
            int mid = neb->images/2 + 1;
            if (image < mid) {
                double frac = ((double)image) / ((double)mid);
                neb->image[image]->setPositions(initial->getPositions() + frac * reactantToTS);
            }else if (image > mid) {
                double frac = (double)(image-mid) / (double)(neb->images - mid + 1);
                neb->image[image]->setPositions(transitionState->getPositions() + frac * TSToProduct);
            }else if (image == mid) {
                neb->image[image]->setPositions(transitionState->getPositions());
            }
        }
    }

    f1 = Potential::fcalls;
    status = neb->compute();
    fCallsNEB += Potential::fcalls - f1;

    if (status == NudgedElasticBand::STATUS_INIT) {
        status = NudgedElasticBand::STATUS_GOOD;
    }

    printEndState(status);
    saveData(status, neb);

    delete neb;
    delete initial;
    delete final;

    return returnFiles;
}

void NudgedElasticBandJob::saveData(int status, NudgedElasticBand *neb)
{
    FILE *fileResults, *fileNEB;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    
    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_neb\n", fCallsNEB);
    fprintf(fileResults, "%f energy_reference\n", neb->image[0]->getPotentialEnergy());
    fprintf(fileResults, "%li number_of_images\n", neb->images);
    for(long i=0; i<=neb->images+1; i++) {
        fprintf(fileResults, "%f image%li_energy\n", neb->image[i]->getPotentialEnergy()-neb->image[0]->getPotentialEnergy(), i);
        fprintf(fileResults, "%f image%li_force\n", neb->image[i]->getForces().norm(), i);
        fprintf(fileResults, "%f image%li_projected_force\n", neb->projectedForce[i]->norm(), i);
    }
    fprintf(fileResults, "%li number_of_extrema\n", neb->numExtrema);
    for(long i=0; i<neb->numExtrema; i++) {
        fprintf(fileResults, "%f extremum%li_position\n", neb->extremumPosition[i], i);
        fprintf(fileResults, "%f extremum%li_energy\n", neb->extremumEnergy[i], i);
    }

    fclose(fileResults);

    std::string nebFilename("neb.con");
    returnFiles.push_back(nebFilename);
    fileNEB = fopen(nebFilename.c_str(), "wb");
    for(long i=0; i<=neb->images+1; i++) {
        neb->image[i]->matter2con(fileNEB);
    }
    fclose(fileNEB);

    returnFiles.push_back("neb.dat");
    neb->printImageData(true);
}

void NudgedElasticBandJob::printEndState(int status)
{
    log("\nFinal state: ");
    if(status == NudgedElasticBand::STATUS_GOOD)
        log("Nudged elastic band, successful.\n");
    else if(status == NudgedElasticBand::STATUS_BAD_MAX_ITERATIONS)
        log("Nudged elastic band, too many iterations.\n");
    else
        log("Unknown status: %i!\n", status);
    return;
}

