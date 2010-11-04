#include "DimerDisplaceJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"
#include "HelperFunctions.h"
#include "EpiCenters.h"
#include "Dimer.h"

using namespace helper_functions;

DimerDisplaceJob::DimerDisplaceJob(Parameters *params)
{
    parameters = params;
}

DimerDisplaceJob::~DimerDisplaceJob(){ }

void DimerDisplaceJob::run(int bundleNumber)
{
    // No bundling for this job, so bundleNumber is ignored.
    
    // Load the displacement con file and get the position.
    
    FILE *results = fopen("results.dat", "w");
    fprintf(results, "%15s   %15s   %15s   %15s   %15s\n", "cutoff", "magnitude", "curvature", "dE", "result");

    int nSteps = 1024; //XXX Parameterize.
    double cutoff = 3.3; //XXX Parameterize.
    double cutoffStep = 0.25;//XXX Parameterize.
    double magnitude = 1.0; //XXX Parameterize.
    double magnitudeStep = 0.25;//XXX Parameterize.
    double iterMax = 32; //XXX Parameterize.
    double torqueConvergence = 0.01; //XXX Parameterize.
    double maxCurvature = -0.1;
    double max_dE = 10.0;

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



























