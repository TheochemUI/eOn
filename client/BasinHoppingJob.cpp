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

#include "BasinHoppingJob.h"
#include "Constants.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "false_boinc.h"
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

  // printf ("%.3f\n", cenx);
  //printf ("%.3f\n", ceny);
  //printf ("%.3f\n", cenz);
  double cen[] = {cenx, ceny, cenz};

  VectorXd dist(num);
  for(int i = 0; i < num; i++)
    {
       
     
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

    }
  return dist;
}

std::vector<std::string> BasinHoppingJob::run(void)
{
    printf("Hello\n");

    Matter *tmpMatter = new Matter(parameters);

    current->con2matter("reactant_passed.con");
    *trial = *current;
    *tmpMatter = *current;

    Minimizer *minimizer; 
    if (parameters->optMethod == "cg") {
        minimizer = new ConjugateGradients(tmpMatter, parameters);
    }else if (parameters->optMethod == "qm"){
        minimizer = new Quickmin(tmpMatter, parameters);
    }
    minimizer->setOutput(0);
    minimizer->fullRelax();
    delete minimizer;

    double currentEnergy = tmpMatter->getPotentialEnergy();
    double minimumEnergy = currentEnergy;

    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;
    int nsteps = parameters->basinHoppingSteps + parameters->basinHoppingQuenchingSteps;

    for (int step=0; step<nsteps; step++)
    {
        AtomMatrix displacement;
       
        displacement = displaceRandom();
 
        trial->setPositions(current->getPositions() + displacement);

        *tmpMatter = *trial;

        if (parameters->optMethod == "cg") {
            minimizer = new ConjugateGradients(tmpMatter, parameters);
        }else if (parameters->optMethod == "qm"){
            minimizer = new Quickmin(tmpMatter, parameters);
        }
        minimizer->setOutput(0);
        minimizer->fullRelax();

        double deltaE = tmpMatter->getPotentialEnergy()-currentEnergy;
        double p;
        if (step==parameters->basinHoppingSteps)
        {
            p = 1.0;
        }else{
            p = exp(-deltaE / (parameters->temperature*8.617343e-5));
        }

        if (randomDouble(1.0)<min(1.0, p)) 
        {
            *current = *trial;
            if(parameters->basinHoppingStayMinimized)
            {
                *current = *tmpMatter;
            }
            currentEnergy = tmpMatter->getPotentialEnergy();
            if (abs(deltaE)>parameters->structureComparisonEnergyDifference) 
            {
                *current = *tmpMatter;
                if (currentEnergy < minimumEnergy) 
                {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
            tmpMatter->matter2xyz("movie", true);
        }
        printf("step: %6i energy: %10.4f c_energy: %12.3e de: %10.2e min_fc: %ld\n",
               step+1, currentEnergy, current->getPotentialEnergy(),
               deltaE, minimizer->totalForceCalls);
        boinc_fraction_done((double)(step+1)/(double)parameters->basinHoppingSteps);
        delete minimizer;
    }

    /* Save Results */

    FILE *fileResults, *fileProduct;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", 0);
    fprintf(fileResults, "%e minimum_energy\n", minimumEnergy);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fclose(fileResults);

    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);
    fileProduct = fopen(productFilename.c_str(), "wb");
    minimumEnergyStructure->matter2con(fileProduct);
    fclose(fileProduct);

    delete tmpMatter;
    delete minimumEnergyStructure;
  
    return returnFiles;
}

AtomMatrix BasinHoppingJob::displaceRandom()
{
    // Create a random displacement.
    AtomMatrix displacement;
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();
    double md = parameters->basinHoppingMaxDisplacement;
    VectorXd distvec = calculateDistanceFromCenter(current);
    int num = trial->numberOfAtoms();
    int m = 0;
    if(parameters->basinHoppingSingleAtomDisplace)
    {
      m = (long)random(trial->numberOfAtoms());
      num = m + 1;
    }

    for(int i = m; i < num; i++)
    {
        double dist = distvec(i);
        double mdp = 0.0;

        if(!trial->getFixed(i))
        {
            if(parameters->basinHoppingMaxDisplacementAlgorithm=="standard")
            {
                mdp = md;
            }
            else if(parameters->basinHoppingMaxDisplacementAlgorithm=="linear")
            {
                double Cs = md/distvec.maxCoeff();
                mdp = Cs*dist;
            }
            else if(parameters->basinHoppingMaxDisplacementAlgorithm=="quadratic")
            {
                double Cq = md/(distvec.maxCoeff()*distvec.maxCoeff());
                mdp = Cq*dist*dist;
            }
            else
            {
                printf("Unknown max_displacement_algorithm\n");
                exit(1);
            }
            for(int j=0; j<3; j++)
            {
                displacement(i, j) = randomDouble(2*mdp) - mdp;
            }
        }
    }
    return displacement;
}
