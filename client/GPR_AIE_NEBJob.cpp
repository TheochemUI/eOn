#include "GPR_AIE_NEBJob.h"
#include "ConjugateGradients.h"
#include "Potential.h"
#include "Log.h"

#include <stdio.h>
#include <string>

GPR_AIE_NEBJob::GPR_AIE_NEBJob(Parameters *parametersPassed) : parameters{parametersPassed},
    fCallsNEB{0}, fCallsGPR{0}
{
}

GPR_AIE_NEBJob::~GPR_AIE_NEBJob()
{}

std::vector<std::string> GPR_AIE_NEBJob::run(void)
{
    long status;
    int f1;

    // TODO: Deal with tsInterpolate, see NudgedElasticBandJob
    string reactantFilename = helper_functions::getRelevantFile("reactant.con");
    string productFilename = helper_functions::getRelevantFile("product.con");

    Parameters eonp = *this->parameters;
    // Setup the Matter objects
    Matter reactant{&eonp}, product{&eonp};
    reactant.con2matter(reactantFilename);
    product.con2matter(productFilename);
    // Setup the observations
    auto imgArray = helper_functions::prepInitialPath(&eonp);
    // Setup GPR
    auto gpf = std::make_shared<GPRobj>(reactant, eonp);
    gpf->trainGPR(imgArray);

   std::vector<GPRMatter> vecgpr;
    for (auto& image : imgArray){
        GPRMatter tmpmatter(image, gpf);
        vecgpr.push_back(tmpmatter);
    }

    auto gpnebInit = GPRNEB(vecgpr, eonp);

    status = gpnebInit.compute();
    this->fCallsGPR += 1;

        // printEndState(status);
        // saveData(status, &gpnebInit);

        // return returnFiles;
    bool mustUpdate = gpnebInit.needsRetraining();

    f1 = Potential::fcalls;
    std::vector<Matter> matvec;
    std::vector<GPRMatter> gpmvec;
    if(mustUpdate){
        matvec = gpnebInit.getCurPath();
    }

    while(mustUpdate){
        gpf->retrainGPR(matvec);
        gpmvec.clear();
        for (auto& image : imgArray){
            GPRMatter tmpmatter(image, gpf);
            gpmvec.push_back(tmpmatter);
        }
        auto gpnebTwo = GPRNEB(gpmvec, eonp);
        gpnebTwo.compute();
        this->fCallsGPR += 1;
        mustUpdate = gpnebTwo.needsRetraining();
        matvec = gpnebTwo.getCurPath();
    };
    // If there is only one round and no updates, do not redo NEB
    if (this->fCallsGPR > 1){
        // Final round
        gpf->retrainGPR(matvec);
        gpmvec.clear();
        for (auto& image : imgArray){
            GPRMatter tmpmatter(image, gpf);
            gpmvec.push_back(tmpmatter);
        }
        auto gpnebFin = GPRNEB(gpmvec, eonp);
        gpnebFin.compute();
        this->fCallsGPR += 1;

        fCallsNEB += Potential::fcalls - f1;

        if (status == GPRNEB::STATUS_INIT) {
            status = GPRNEB::STATUS_GOOD;
        }

        printEndState(status);
        saveData(status, &gpnebFin);

        return returnFiles;
    } else {
        // TODO: Be more elegant when one round has occurred
        fCallsNEB += Potential::fcalls - f1;

        if (status == GPRNEB::STATUS_INIT) {
            status = GPRNEB::STATUS_GOOD;
        }

        printEndState(status);
        saveData(status, &gpnebInit);

        return returnFiles;
    }
}

void GPR_AIE_NEBJob::saveData(int status, GPRNEB *gpneb)
{
    FILE *fileResults, *fileNEB;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_neb\n", fCallsNEB);
    fprintf(fileResults, "%d gpr_created\n", fCallsGPR);
    auto pef_init = gpneb->imageArray.front().gpr_energy_forces();
    fprintf(fileResults, "%f energy_reference\n", std::get<double>(pef_init));
    fprintf(fileResults, "%li number_of_images\n", gpneb->nimages);
    for(size_t idx{0}; idx <= gpneb->nimages+1; idx++) {
        auto pef_cur = gpneb->imageArray[idx].gpr_energy_forces();
        fprintf(fileResults, "%f image%li_energy\n", std::get<double>(pef_cur)-std::get<double>(pef_init), idx);
        fprintf(fileResults, "%f image%li_force\n", std::get<AtomMatrix>(pef_cur).norm(), idx);
        fprintf(fileResults, "%f image%li_projected_force\n", gpneb->projectedForceArray[idx].norm(), idx);
    }
    fprintf(fileResults, "%li number_of_extrema\n", gpneb->numExtrema);
    for(size_t idx{0}; idx<gpneb->numExtrema; idx++) {
        fprintf(fileResults, "%f extremum%li_position\n", gpneb->extremumPositions[idx], idx);
        fprintf(fileResults, "%f extremum%li_energy\n", gpneb->extremumEnergies[idx], idx);
    }

    fclose(fileResults);

    std::string gpnebFilename("neb.con");
    returnFiles.push_back(gpnebFilename);
    fileNEB = fopen(gpnebFilename.c_str(), "wb");
    for(size_t idx{0}; idx <= gpneb->nimages+1; idx++) {
        gpneb->imageArray[idx].truePotMatter.matter2con(fileNEB);
    }
    fclose(fileNEB);

    returnFiles.push_back("neb.dat");
    gpneb->printImageData(true);
}

void GPR_AIE_NEBJob::printEndState(int status)
{
    log("\nFinal state: ");
    if(status == GPRNEB::STATUS_GOOD)
        log("GPR-AIE Nudged elastic band, successful.\n");
    else if(status == GPRNEB::STATUS_BAD_MAX_ITERATIONS)
        log("GPR-AIE Nudged elastic band, too many iterations.\n");
    else
        log("Unknown status: %i!\n", status);
    return;
}
