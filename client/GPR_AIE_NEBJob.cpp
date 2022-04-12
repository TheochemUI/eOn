#include "GPR_AIE_NEBJob.h"
#include "ConjugateGradients.h"
#include "Potential.h"
#include "Log.h"

#include <stdio.h>
#include <string>

GPR_AIE_NEBJob::GPR_AIE_NEBJob(Parameters *parametersPassed)
{
    parameters = parametersPassed;
    fCallsNEB = 0;
    gprfunc = std::make_unique<gpr::GaussianProcessRegression>();
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

    auto matterInit = std::make_unique<Matter>(this->parameters);
    auto matterFinal = std::make_unique<Matter>(this->parameters);

    matterInit->con2matter(reactantFilename);
    matterFinal->con2matter(productFilename);

    // Prepare data for GPR
    auto config_data = helper_functions::eon_matter_to_frozen_conf_info(matterInit.get(),
                                                                        this->parameters->gprPotActiveRadius);
    auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
    // Setup the runs
    auto imgArray = helper_functions::prepInitialPath(this->parameters);
    auto obspath = helper_functions::prepInitialObs(imgArray);
    // Setup GPR
    auto eondat = std::make_pair(*this->parameters, *matterInit);
    *this->gprfunc = helper_functions::initializeGPR(*this->gprfunc, atoms_config, obspath, eondat);
    this->gprfunc->setHyperparameters(obspath, atoms_config);
    this->gprfunc->optimize(obspath);
    // Prepare GPR potential
    GPRPotential gprpot{this->parameters};
    gprpot.registerGPRObject(this->gprfunc.get());
    // Set Potentials
    matterInit->setPotential(&gprpot);
    matterFinal->setPotential(&gprpot);

    auto nebInit = helper_functions::prepGPRNEBround(*this->gprfunc,
                                                *matterInit, *matterFinal,
                                                *this->parameters);
    nebInit->compute();

    bool mustUpdate = helper_functions::maybeUpdateObs(*nebInit, obspath);

    f1 = Potential::fcalls;

    while(mustUpdate){
        this->gprfunc->setHyperparameters(obspath, atoms_config, false);
        this->gprfunc->optimize(obspath);
        auto nebTwo = helper_functions::prepGPRNEBround(*this->gprfunc,
                                                        *matterInit, *matterFinal,
                                                        *this->parameters);
        nebTwo->compute();
        mustUpdate = helper_functions::maybeUpdateObs(*nebTwo, obspath);
    };
    // Final round
    auto nebFin = helper_functions::prepGPRNEBround(*this->gprfunc,
                                                    *matterInit, *matterFinal,
                                                    *this->parameters);
    status = nebFin->compute();

    fCallsNEB += Potential::fcalls - f1;

    if (status == NudgedElasticBand::STATUS_INIT) {
        status = NudgedElasticBand::STATUS_GOOD;
    }

    printEndState(status);
    saveData(status, nebFin.get());

    return returnFiles;
}

void GPR_AIE_NEBJob::saveData(int status, NudgedElasticBand *neb)
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

void GPR_AIE_NEBJob::printEndState(int status)
{
    log("\nFinal state: ");
    if(status == NudgedElasticBand::STATUS_GOOD)
        log("GPR-AIE Nudged elastic band, successful.\n");
    else if(status == NudgedElasticBand::STATUS_BAD_MAX_ITERATIONS)
        log("Nudged elastic band, too many iterations.\n");
    else
        log("Unknown status: %i!\n", status);
    return;
}
