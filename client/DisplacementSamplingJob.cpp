//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "DisplacementSamplingJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"
#include "HelperFunctions.h"
#include "EpiCenters.h"
#include "Dimer.h"
#include <math.h>

using namespace helper_functions;

DisplacementSamplingJob::DisplacementSamplingJob(Parameters *params)
{
    parameters = params;
}

DisplacementSamplingJob::~DisplacementSamplingJob(){ }

void DisplacementSamplingJob::run(int bundleNumber)
{
    // No bundling for this job, so bundleNumber is ignored.
    
    FILE *results = fopen("results.dat", "w");
    fprintf(results, "%15s   %15s   %15s\n", "cutoff", "magnitude", "fraction good");

    long   iterMax = parameters->displaceIterMax;
    double torqueConvergence = parameters->displaceTorqueConvergence;
    double maxCurvature = parameters->displaceMaxCurvature;
    double max_dE = parameters->displaceMaxDE;
    int nSamples = parameters->displaceNSamples;;

    parameters->dimerRotationsMax = 1;
    parameters->dimerRotationsMin = 1;
    parameters->dimerWindowMax = 0.0;
    parameters->dimerWindowMin = 0.0;

    std::stringstream cutoffStream(parameters->displaceCutoffs);
    while(true)
    {
        double cutoff;
        cutoffStream >> cutoff;
        if(!cutoffStream)
        {
            break;
        }
        std::stringstream magnitudeStream(parameters->displaceMagnitudes);
        while(true)
        {
            double magnitude;
            magnitudeStream >> magnitude;
            if(!magnitudeStream)
            {
                break;
            }
    
            double nGood = 0;

            for(int k = 0; k < nSamples; k++)
            {
                Matter *reactant = new Matter(parameters);
                reactant->con2matter("reactant_passed.con");
                double e0 = reactant->getPotentialEnergy();

                long epicenter = EpiCenters::minCoordinatedEpiCenter(reactant,parameters->neighborCutoff);
                
                // Create a random displacement.
                Matrix<double, Eigen::Dynamic, 3> displacement;        
                displacement.resize(reactant->numberOfAtoms(), 3);
                displacement.setZero();
                for(int i = 0; i < reactant->numberOfAtoms(); i++)
                {
                    if(reactant->distance(epicenter, i) < cutoff)
                    {
                        for(int j = 0; j < 3; j++)
                        {
                            if(!reactant->getFixed(i))
                            {
                                displacement(i, j) = randomDouble(1.0);
                            }
                        }
                    }
                }
                
                displacement.normalize(); 
                displacement *= gaussRandom(0.0, magnitude);
                
                reactant->setPositions(reactant->getPositions() + displacement);
                
                Dimer *d = new Dimer(reactant, parameters);
                
                d->startNewSearchAndCompute(reactant, displacement);
                d->moveAndCompute(reactant);
                int iterCount = 0;
                while(d->stats[0] > torqueConvergence && iterCount < iterMax)
                {
                    d->moveAndCompute(reactant);
                    iterCount++;
                }
                double e1 = reactant->getPotentialEnergy();
                if(e1 - e0 < max_dE && d->stats[1] < maxCurvature)
                {
                    nGood++;
                }

                delete d;
                delete reactant;
            }
            fprintf(results, "%15f   %15f   %15f\n", cutoff, magnitude, nGood / nSamples);
        }
    }
    fclose(results);
}


