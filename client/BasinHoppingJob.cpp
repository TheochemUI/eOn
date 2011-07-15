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
#include "Constants.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "Potential.h"
#include "HelperFunctions.h"

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
}

BasinHoppingJob::~BasinHoppingJob()
{
    delete current;
    delete trial;
}

vector<long> BasinHoppingJob::getElements(Matter *matter)
{
    int allElements[118] = {0};
    vector<long> Elements;

    for(long y=0; y<matter->numberOfAtoms(); y++)
    {
        if(!matter->getFixed(y))
        {
            int index= matter->getAtomicNr(y);
            allElements[index]++;
        }
    }
    for(int i=0; i<118; i++)
    {
        if(allElements[i] > 0)
        {
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
    for(int k=0; k<matter->numberOfAtoms(); k++)
    {
        cenx = cenx + pos(k,0);
        ceny = ceny + pos(k,1);
        cenz = cenz + pos(k,2);
    }
    int num = matter->numberOfAtoms();
    cenx = cenx / num;
    ceny = ceny / num;
    cenz = cenz / num;

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

std::vector<std::string> BasinHoppingJob::run(void)
{
    bool swapMove;
    double swap_accept=0.0;
    jcount=0;
    scount=0;
    dcount=0;
    int jump_max_count=0;
    int jump_steps_count=0;
    double totalAccept=0.0;
    Matter *minTrial = new Matter(parameters);
    Matter *swapTrial = new Matter(parameters);

    current->con2matter("reactant_passed.con");
    if(parameters->basinHoppingMDFirst==true){
        Dynamics dyn(current,parameters);
        dyn.fullSteps(parameters->basinHoppingMDTemp);
    }
    *trial = *current;
    *minTrial = *current;


    Minimizer *minimizer = NULL; 
    Parameters minParameters = *parameters;
    minParameters.writeMovies = false;
    if (parameters->optMethod == "cg") {
        minimizer = new ConjugateGradients(current, &minParameters);
    }else if (parameters->optMethod == "qm"){
        minimizer = new Quickmin(current, &minParameters);
    }
    minimizer->setOutput(0);
    minimizer->fullRelax();
    delete minimizer;

    double currentEnergy = current->getPotentialEnergy();
    double minimumEnergy = currentEnergy;

    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;
    int nsteps = parameters->basinHoppingSteps + parameters->basinHoppingQuenchingSteps;
    long totalfc = 0;
    FILE * pFile;
    pFile = fopen("bh.dat","w");

    printf("%4s %12s %12s %12s %4s\n", "mcs", "current", "trial", "global min", "fc");
    printf("%4s %12s %12s %12s %4s\n", "---", "-------", "-----", "----------", "--");

    for (int step=0; step<nsteps; step++) {
        if(randomDouble(1.0)<parameters->basinHoppingSwapProbability && 
           step<parameters->basinHoppingSteps && 
           jump_max_count<parameters->basinHoppingJumpMax){
            *swapTrial = *current;
            randomSwap(swapTrial);
            swapMove=true;
            *minTrial = *swapTrial;
      	}else{
            AtomMatrix displacement;
            displacement = displaceRandom();

            trial->setPositions(current->getPositions() + displacement);
            swapMove=false;
            *minTrial = *trial;
        }

        if (parameters->writeMovies == true) {
            trial->matter2xyz("trials", true);
        }

        if (parameters->optMethod == "cg") {
            minimizer = new ConjugateGradients(minTrial, &minParameters);
        }else if (parameters->optMethod == "qm"){
            minimizer = new Quickmin(minTrial, &minParameters);
        }
        minimizer->setOutput(0);
        minimizer->fullRelax();
        double deltaE = minTrial->getPotentialEnergy()-currentEnergy;
        double p=0.0;
        if (step>=parameters->basinHoppingSteps) {
            if (deltaE < 0.0) {
                    p = 1.0;
            }
        }else if (jump_max_count>=parameters->basinHoppingJumpMax) {
            jump_steps_count++;
            jcount++;
            if (deltaE > 0.0) {
                p = 1.0;
            }
        }else{
            if (deltaE < 0.0) {
                p = 1.0;
            }else{
                p = exp(-deltaE / (parameters->temperature*8.617343e-5));
            }
        }

        if (randomDouble(1.0)<min(1.0, p)) {
            if(parameters->basinHoppingSignificantStructure){
                *current = *minTrial;
            }else{
                *current = *trial;
            }
            if(swapMove){
                swap_accept=swap_accept+1.0;
            }
            if(step<parameters->basinHoppingSteps) {
                totalAccept=totalAccept+1.0;
            }
       
            currentEnergy = minTrial->getPotentialEnergy();
            if (abs(deltaE)>parameters->structureComparisonEnergyDifference) {
                *current = *minTrial;
                if (currentEnergy < minimumEnergy) 
                {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
        }else{
            jump_max_count++;
        }

        if (jump_steps_count==parameters->basinHoppingJumpSteps) {
            jump_max_count=0;
            jump_steps_count=0;
        }

        if (parameters->writeMovies == true) {
            minTrial->matter2xyz("movie", true);
        }

        totalfc = totalfc + minimizer->totalForceCalls;
        //printf("mcs: %4i e_cur: %.3f e_trial: %.3f g_min: %.3f fc: %4ld\n",
        printf("%4i %12.3f %12.3f %12.3f %4ld\n",
               step+1, currentEnergy, minTrial->getPotentialEnergy(), minimumEnergy,
               minimizer->totalForceCalls);
        fprintf(pFile, "%6i %9ld %12.4e %12.4e\n",step+1,totalfc,currentEnergy,
                minTrial->getPotentialEnergy());

        boinc_fraction_done(((double)step+1.0)/(double)nsteps);
        delete minimizer;
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
    fprintf(fileResults, "%e minimum_energy\n", minimumEnergy);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%.3f acceptance_ratio\n", totalAccept/parameters->basinHoppingSteps);
    if(parameters->basinHoppingSwapProbability>0){
      fprintf(fileResults, "%.3f swap_acceptance_ratio\n", swap_accept/double(scount));
    }
    fprintf(fileResults, "%ld total_normal_displacement_steps\n",dcount-jcount-parameters->basinHoppingQuenchingSteps);
    fprintf(fileResults, "%d total_jump_steps\n", jcount);
    fprintf(fileResults, "%d total_swap_steps\n", scount);
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
      m = (long)random(trial->numberOfAtoms());
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
    int i=0;
    int j=0;
    while(j+i!=2) {
        for(long x=randomInt(0,matter->numberOfAtoms()-1); x<matter->numberOfAtoms(); x++){
            if(matter->getAtomicNr(x)==ela) {
                changera=x;
                i=1;
                break;
            }
        }
        for(long x=randomInt(0, matter->numberOfAtoms()-1); x<matter->numberOfAtoms(); x++){
            if(matter->getAtomicNr(x)==elb) {
                changerb=x;
                j=1;
                break;
            }
        }
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
