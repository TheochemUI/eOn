//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <stdio.h>
#include <string>
#include <math.h>

#include "Dynamics.h"
#include "BasinHoppingJob.h"
#include "Potential.h"
#include "HelperFunctions.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "Log.h"

#ifdef BOINC
    #include <boinc/boinc_api.h>
    #include <boinc/diagnostics.h>
    #include <boinc/filesys.h>
#ifdef WIN32
    #include <boinc/boinc_win.h>
    #include <boinc/win_util.h>
#endif
#else
    #include "false_boinc.h"
#endif

using namespace std;
using namespace helper_functions;

BasinHoppingJob::BasinHoppingJob(Parameters *params)
{
    parameters = params;
    current = new Matter(parameters);
    trial = new Matter(parameters);
    fcalls = Potential::fcalls;
}

BasinHoppingJob::~BasinHoppingJob()
{
    delete current;
    delete trial;
}

std::vector<std::string> BasinHoppingJob::run(void)
{
    bool swapMove;
    double swap_accept = 0.0;
    jcount = 0;
    scount = 0;
    dcount = 0;
    int consecutive_rejected_trials = 0;
    double totalAccept = 0.0;
    Matter *minTrial = new Matter(parameters);
    Matter *swapTrial = new Matter(parameters);

//    current->con2matter("reactant_passed.con");
    string conFilename = getRelevantFile(parameters->conFilename);
    current->con2matter(conFilename);
    if(parameters->basinHoppingInitialMD == true){
        Dynamics dyn(current,parameters);
        dyn.setTemperature(parameters->basinHoppingInitialMDTemperature);
        dyn.run();
    }
    *trial = *current;
    *minTrial = *current;

    current->relax(true);

    double currentEnergy = current->getPotentialEnergy();
    double minimumEnergy = currentEnergy;

    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;
    int nsteps = parameters->basinHoppingSteps + parameters->basinHoppingQuenchingSteps;
    long totalfc;
    FILE * pFile;
    pFile = fopen("bh.dat","w");

    log("%4s %12s %12s %12s %4s\n", "mcs", "current", "trial", "global min", "fc");
    log("%4s %12s %12s %12s %4s\n", "---", "-------", "-----", "----------", "--");

    for (int step=0; step<nsteps; step++) {
        if(randomDouble(1.0)<parameters->basinHoppingSwapProbability && 
           step<parameters->basinHoppingSteps){
            *swapTrial = *current;
            randomSwap(swapTrial);
            swapMove = true;
            *minTrial = *swapTrial;
        }else{
            AtomMatrix displacement;
            displacement = displaceRandom();

            trial->setPositions(current->getPositions() + displacement);
            swapMove=false;
            *minTrial = *trial;
        }

        if (parameters->writeMovies == true) {
            trial->matter2con("trials", true);
        }

        Potential::fcalls = 0;
        minTrial->relax(true);
        int minfcalls = Potential::fcalls;

        double deltaE = minTrial->getPotentialEnergy()-currentEnergy;
        double p=0.0;
        if (step>=parameters->basinHoppingSteps) {
            if (deltaE <= 0.0) {
                    p = 1.0;
            }
        }else{
            if (deltaE <= 0.0) {
                p = 1.0;
            }else{
                p = exp(-deltaE / (parameters->temperature*8.617343e-5));
            }
        }

        if (randomDouble(1.0) < p) {
            if(parameters->basinHoppingSignificantStructure){
                *current = *minTrial;
            }else{
                *current = *trial;
            }
            if(swapMove){
                swap_accept = swap_accept + 1.0;
            }
            if(step<parameters->basinHoppingSteps) {
                totalAccept = totalAccept + 1.0;
            }

            currentEnergy = minTrial->getPotentialEnergy();
            if (currentEnergy < minimumEnergy) {
                minimumEnergy = currentEnergy;
                *minimumEnergyStructure = *minTrial;
            }
        }else{
            consecutive_rejected_trials++;
        }

        if (parameters->writeMovies == true) {
            minTrial->matter2con("movie", true);
        }

        totalfc = Potential::fcallsTotal;
        log("%4i %12.3f %12.3f %12.3f %4i\n",
               step+1, currentEnergy, minTrial->getPotentialEnergy(), minimumEnergy,
               minfcalls);
        fprintf(pFile, "%6i %9ld %12.4e %12.4e\n",step+1,totalfc,currentEnergy,
                minTrial->getPotentialEnergy());

        if(consecutive_rejected_trials == parameters->basinHoppingJumpMax && step<parameters->basinHoppingSteps){
            consecutive_rejected_trials = 0;
            AtomMatrix jump;
            for(int j=0; j<parameters->basinHoppingJumpSteps; j++){
                jcount++;
                jump = displaceRandom();
                current->setPositions(current->getPositions() + jump);
                if(parameters->basinHoppingSignificantStructure){
                    current->relax(true);
                }
                currentEnergy = current->getPotentialEnergy();
                if (currentEnergy < minimumEnergy) {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
        }
        boinc_fraction_done(((double)step+1.0)/(double)nsteps);
    }
    fclose(pFile);

    /* Save Results */

    FILE *fileResults, *fileProduct;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    if (parameters->writeMovies == true) {
        std::string movieFilename("movie.xyz");
        returnFiles.push_back(movieFilename);
    }

    fprintf(fileResults, "%d termination_reason\n", 0);
    fprintf(fileResults, "%.6f minimum_energy\n", minimumEnergy);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%.3f acceptance_ratio\n", totalAccept/parameters->basinHoppingSteps);
    if(parameters->basinHoppingSwapProbability>0){
        fprintf(fileResults, "%.3f swap_acceptance_ratio\n", swap_accept/double(scount));
    }
    fprintf(fileResults, "%ld total_normal_displacement_steps\n",dcount-jcount-parameters->basinHoppingQuenchingSteps);
    fprintf(fileResults, "%d total_jump_steps\n", jcount);
    fprintf(fileResults, "%d total_swap_steps\n", scount);
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    fclose(fileResults);

    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);
    fileProduct = fopen(productFilename.c_str(), "wb");
    minimumEnergyStructure->matter2con(fileProduct);
    fclose(fileProduct);

    std::string bhFilename("bh.dat");
    returnFiles.push_back(bhFilename);

    delete minTrial;
    delete minimumEnergyStructure;
    delete swapTrial;
    return returnFiles;
}

AtomMatrix BasinHoppingJob::displaceRandom()
{
    dcount++;
    // Create a random displacement.
    AtomMatrix displacement;
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();
    double md = parameters->basinHoppingMaxDisplacement;
    VectorXd distvec = calculateDistanceFromCenter(current);
    int num = trial->numberOfAtoms();
    int m = 0;
    if(parameters->basinHoppingSingleAtomDisplace) {
      m = randomInt(0, trial->numberOfAtoms()-1);
      num = m + 1;
    }

    for(int i = m; i < num; i++) {
        double dist = distvec(i);
        double mdp = 0.0;

        if(!trial->getFixed(i)) {
            if(parameters->basinHoppingMaxDisplacementAlgorithm=="standard") {
                mdp = md;
            }
            else if(parameters->basinHoppingMaxDisplacementAlgorithm=="linear") {
                double Cs = md/distvec.maxCoeff();
                mdp = Cs*dist;
            }
            else if(parameters->basinHoppingMaxDisplacementAlgorithm=="quadratic") {
                double Cq = md/(distvec.maxCoeff()*distvec.maxCoeff());
                mdp = Cq*dist*dist;
            }else{
                printf("Unknown max_displacement_algorithm\n");
                exit(1);
            }
            for(int j=0; j<3; j++) {
                if(parameters->basinHoppingDisplacementDistribution=="uniform") {
                    displacement(i, j) = randomDouble(2*mdp) - mdp;
                }
                else if(parameters->basinHoppingDisplacementDistribution=="gaussian") {
                    displacement(i,j) = gaussRandom(0.0, mdp);
                }
                else {
                    printf("Unknown displacement_distribution\n");
                    exit(1);
                }
            }
        }
    }
    return displacement;
}


void BasinHoppingJob::randomSwap(Matter *matter)
{
    scount++;
    vector<long> Elements;
    Elements=getElements(matter);

    long ela;
    long elb;
    long ia = randomInt(0, Elements.size()-1);
    ela = Elements.at(ia);
    Elements.erase(Elements.begin()+ia);

    long ib = randomInt(0, Elements.size()-1);
    elb = Elements.at(ib);

    int changera=0;
    int changerb=0;

    changera = randomInt(0, matter->numberOfAtoms()-1);
    while (matter->getAtomicNr(changera) != ela) {
        changera = randomInt(0, matter->numberOfAtoms()-1);
    }

    changerb = randomInt(0, matter->numberOfAtoms()-1);
    while (matter->getAtomicNr(changerb) != elb) {
        changerb = randomInt(0, matter->numberOfAtoms()-1);
    }

    double posax=matter->getPosition(changera, 0);
    double posay=matter->getPosition(changera, 1);
    double posaz=matter->getPosition(changera, 2);

    matter->setPosition(changera, 0, matter->getPosition(changerb, 0));
    matter->setPosition(changera, 1, matter->getPosition(changerb, 1));
    matter->setPosition(changera, 2, matter->getPosition(changerb, 2));

    matter->setPosition(changerb, 0, posax);
    matter->setPosition(changerb, 1, posay);
    matter->setPosition(changerb, 2, posaz);
}

vector<long> BasinHoppingJob::getElements(Matter *matter)
{
    int allElements[118] = {0};
    vector<long> Elements;

    for(long y=0; y<matter->numberOfAtoms(); y++) {
        if(!matter->getFixed(y)) {
            int index= matter->getAtomicNr(y);
            allElements[index] = 1;
        }
    }

    for(int i=0; i<118; i++) {
        if(allElements[i] != 0) {
            Elements.push_back(i);
        }
    }

    return Elements;
}

VectorXd BasinHoppingJob::calculateDistanceFromCenter(Matter *matter)
{
    AtomMatrix pos = matter->getPositions();
    double cenx = 0;
    double ceny = 0;
    double cenz = 0;
    int num = matter->numberOfAtoms();

    for(int k=0; k<num; k++) {
        cenx = cenx + pos(k,0);
        ceny = ceny + pos(k,1);
        cenz = cenz + pos(k,2);
    }
    cenx = cenx / (double)num;
    ceny = ceny / (double)num;
    cenz = cenz / (double)num;

    double cen[] = {cenx, ceny, cenz};

    VectorXd dist(num);
     
    for(int l = 0; l < num; l++)
    {
        double xd = pos(l,0)-cen[0];
        double yd = pos(l,1)-cen[1];
        double zd = pos(l,2)-cen[2];
        xd = xd*xd;
        yd = yd*yd;
        zd = zd*zd;
        dist(l) = sqrt(xd + yd + zd);
    }

    return dist;
}
