#include "DisplacementSamplingJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"
#include "HelperFunctions.h"
#include "EpiCenters.h"
#include "Dimer.h"

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
    fprintf(results, "%15s   %15s   %15s   %15s   %15s\n", "cutoff", "magnitude", "curvature", "dE", "result");

    long   nSteps = parameters->displaceNSteps;
    double cutoff = parameters->displaceCutoff;
    double cutoffStep = parameters->displaceCutoffStep;
    double magnitude = parameters->displaceMagnitude;
    double magnitudeStep = parameters->displaceMagnitudeStep;
    long   iterMax = parameters->displaceIterMax;
    double torqueConvergence = parameters->displaceTorqueConvergence;
    double maxCurvature = parameters->displaceMaxCurvature;
    double max_dE = parameters->displaceMaxDE;

    parameters->dimerRotationsHigh = 1;
    parameters->dimerRotationsLow = 1;
    parameters->dimerWindowHigh = 0.0;
    parameters->dimerWindowLow = 0.0;

    for(int k = 0; k < nSteps; k++)
    {

        double newCutoff = guaRandom(cutoff, cutoffStep);
        double newMagnitude = guaRandom(magnitude, magnitudeStep);

        Matter *reactant = new Matter(parameters);
        reactant->con2matter("reactant_passed.con");
        double e0 = reactant->getPotentialEnergy();

        long epicenter = EpiCenters::minimalCoordinatedEpiCenter(reactant);
        
        
        // Create a random displacement.
        Matrix<double, Eigen::Dynamic, 3> displacement;        
        displacement.resize(reactant->numberOfAtoms(), 3);
        displacement.setZero();
        for(int i = 0; i < reactant->numberOfAtoms(); i++)
        {
            if(reactant->distance(epicenter, i) < newCutoff)
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
        displacement *= guaRandom(0.0, newMagnitude);
        
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
            cutoff = newCutoff;
            magnitude = newMagnitude;
            fprintf(results, "%15f   %15f   %15f   %15f   %15s\n", newCutoff, newMagnitude, d->stats[1], e1 - e0, "good");
        }
        else
        {
            fprintf(results, "%15f   %15f   %15f   %15f   %15s\n", newCutoff, newMagnitude, d->stats[1], e1 - e0, "bad");
        }

        delete d;
        delete reactant;
        
    }
    
    fclose(results);
}



























