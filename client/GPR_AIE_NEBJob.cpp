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
    std::vector<Matter> matvec;
    std::vector<Matter> prevpath;
    double trueConvForce{999};
    bool stoppedEarly{false}, converged{false}, mustUpdate{true};
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
    while(!converged){
        // We are in the outer loop
        auto gpnebInit = GPRNEB(helper_functions::prepGPRMatterVec(imgArray, gpf), eonp);
        status = gpnebInit.compute();
        this->fCallsGPR += 1;
        trueConvForce = gpnebInit.getTrueConvForce();
        std::cout<<"\n Current trueConvForce is "<<trueConvForce<<" and needs to be "<<eonp.gprPotTol<<std::endl;
        converged = (trueConvForce<eonp.gprPotTol);
        mustUpdate = gpnebInit.needsRetraining(eonp.gprPotTol);
        f1 = Potential::fcalls;
        if(mustUpdate){
            matvec = gpnebInit.getCurPath();
            prevpath = gpnebInit.getCurPathFull();
        }
        while(mustUpdate){
            // Relaxation loop
            gpf->retrainGPR(matvec);
            auto gpnebTwo = GPRNEB(helper_functions::prepGPRMatterVec(imgArray, gpf), eonp);
            gpnebTwo.compute();
            this->fCallsGPR += 1;
            mustUpdate = gpnebTwo.needsRetraining(eonp.gprPotTol);
            trueConvForce = gpnebTwo.getTrueConvForce();
            std::cout<<"\n Current trueConvForce is "<<trueConvForce<<" and needs to be "<<eonp.gprPotTol<<std::endl;
            converged = (trueConvForce<eonp.gprPotTol);
            // TODO: Use the discount factor in parameters instead of 0.5
            stoppedEarly = gpnebTwo.stoppedEarly(prevpath, 0.5);
            if (stoppedEarly){
                std::cout<<"\nStopping relaxation phase due to early stopping\n";
                std::cout<<"Not using current path\n";
                break;
            }
            if (converged){
                std::cout<<"\nConverged within relaxation phase\n";
                if (status == GPRNEB::STATUS_INIT) {
                    status = GPRNEB::STATUS_GOOD;
                }
                printEndState(status);
                saveData(status, &gpnebTwo);
                return returnFiles;
            }
            if (mustUpdate){
                std::cout<<"\nHaven't bailed will update\n";
                std::cout<<"Adding current path to training data\n";
                matvec = gpnebTwo.getCurPath();
                prevpath = gpnebTwo.getCurPathFull();
            }
        }
        // Outer loop
        if (converged){
            std::cout<<"\nConverged within outer loop\n";
            if (status == GPRNEB::STATUS_INIT) {
                status = GPRNEB::STATUS_GOOD;
            }
            printEndState(status);
            saveData(status, &gpnebInit);
            return returnFiles;
        }
    }
    throw std::runtime_error("You shouldn't be here"s);
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
