#include "GPR_AIE_NEBJob.h"
#include "ConjugateGradients.h"
#include "Potential.h"
#include "Log.h"

#include <stdio.h>
#include <string>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/os.h>

GPR_AIE_NEBJob::GPR_AIE_NEBJob(Parameters *parametersPassed) : eonp{parametersPassed},
                                                               fCallsNEB{0}, fCallsGPR{1},
                                                               reactant{parametersPassed},
                                                               product{parametersPassed},
                                                               stoppedEarly{false},
                                                               converged{false},
                                                               mustUpdate{true},
                                                               isWithin{true},
                                                               useCI{false}

{
    reactantFilename = helper_functions::getRelevantFile("reactant.con");
    productFilename = helper_functions::getRelevantFile("product.con");
    // Setup the Matter objects
    reactant.con2matter(reactantFilename);
    product.con2matter(productFilename);
    gpf = std::make_shared<GPRobj>(reactant, *eonp);
    linearMatter = helper_functions::prepInitialPath(eonp);
    gpf->trainGPR(linearMatter);
    linearPath = helper_functions::prepGPRMatterVec(linearMatter, gpf);
    std::copy(linearMatter.begin()+1, linearMatter.end()-1,
              std::back_inserter(evaluatedIntermediates));
}

GPR_AIE_NEBJob::~GPR_AIE_NEBJob()
{}

void GPR_AIE_NEBJob::retrainGPR(std::vector<Matter>& newpath){
    this->gpf->retrainGPR(newpath);
    ++this->fCallsGPR;
    for (auto mat : newpath){ this->evaluatedIntermediates.push_back(mat);}
}

void GPR_AIE_NEBJob::runGPRNEB(GPRNEB& gprneb){
    int status;
    status = gprneb.compute(this->evaluatedIntermediates,
                            this->isWithin);
    ++this->fCallsNEB;
    saveData(0, &gprneb); // Force to STATUS_GOOD
    this->matvec = gprneb.getCurPath(); // Set current intermediates
    this->checkConvergence(gprneb.getConvergenceTrue()); // Set convergence
    if (!this->converged){
        this->mustUpdate = gprneb.needsRetraining(this->eonp->gprPotTol); // Set retraining
    }
}

void GPR_AIE_NEBJob::checkConvergence(std::pair<double, double> curTrueEnergy){
    if (fCallsGPR < 2 ) {return;}
    const auto [saddleConv, pathConv] = curTrueEnergy;
    const double saddleCriteria = this->eonp->spConvergence;
    const double pathCriteria = this->eonp->mepConvergence; // TODO: Let the path relaxation be a user parameter
    auto fstring = fmt::format("\n For GPR round {} and NEB round {}\n Current saddle force norm is {} and needs to be below {}\n Current path force convergence is {} and needs to be below {}\n",
                               fCallsGPR, fCallsNEB, saddleConv,
                               saddleCriteria, pathConv, pathCriteria);
    std::cout<<fstring;
    this->handleCI(pathConv < pathCriteria);
    if (saddleConv < saddleCriteria && pathConv < pathCriteria){
        this->converged = true;
        std::cout<<"\nCONVERGED\n";
    } else { this->converged = false; }
}

void GPR_AIE_NEBJob::handleCI(bool comparer){
        if (eonp->toggleCI){
            useCI = comparer;
            eonp->nebClimbingImageMethod = comparer;
            eonp->nebClimbingImageConvergedOnly = false;
        }
}

void GPR_AIE_NEBJob::runOuterLoop(bool retrain){
    if (retrain){
        this->retrainGPR(this->matvec);
    }
    auto gpnebInit = GPRNEB(this->linearPath, *eonp);
    this->runGPRNEB(gpnebInit);
}

void GPR_AIE_NEBJob::runRelaxationLoop(){
    this->retrainGPR(this->matvec);
    auto gpneb = GPRNEB(this->linearPath, *eonp);
    this->runGPRNEB(gpneb);
    if(!this->isWithin){
        this->runOuterLoop(true);
    }
}

std::vector<std::string> GPR_AIE_NEBJob::run(void)
{
    long status;
    double trueConvForce;

    // TODO: Deal with tsInterpolate, see NudgedElasticBandJob

    bool retrain{false};
    while(!this->converged){
        this->runOuterLoop(retrain);
        std::cout<<"\nFinished the outer loop";
        if (!this->converged){
            this->runRelaxationLoop();
            std::cout<<"\nFinished the relaxation loop";
            retrain = true;
        }
    }
    return returnFiles;
}

void GPR_AIE_NEBJob::saveData(int status, GPRNEB *gpneb)
{
    FILE *fileResults, *fileNEB;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%s potential_type\n", eonp->potential.c_str());
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

    auto gpnebFilename = fmt::format("neb_{:02d}.con", this->fCallsNEB);
    returnFiles.push_back(gpnebFilename);
    fileNEB = fopen(gpnebFilename.c_str(), "wb");
    for(size_t idx{0}; idx <= gpneb->nimages+1; idx++) {
        gpneb->imageArray[idx].truePotMatter.matter2con(fileNEB);
    }
    fclose(fileNEB);

    // TODO: Update return files with all names
    returnFiles.push_back("neb.dat");
    gpneb->printImageDataTrue(true, this->fCallsNEB);
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
