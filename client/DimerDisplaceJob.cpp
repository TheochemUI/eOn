#include "DimerDisplaceJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"
#include "HelperFunctions.h"
#include "EpiCenters.h"

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
    Matter *reactant = new Matter(parameters);
    reactant->con2matter("reactant_passed.con");

    long epicenter = EpiCenters::minimalCoordinatedEpiCenter(reactant);
    printf("\nepicenter: %d\n\n", epicenter);

    // Create a random displacement.
    Matrix<double, Eigen::Dynamic, 3> displacement;        
    displacement.resize(reactant->numberOfAtoms(), 3);
    displacement.setZero();
    for(int i = 0; i < reactant->numberOfAtoms(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(!reactant->getFixed(i))
            {
                displacement(i, j) = randomDouble(1.0);
            }
        }
    }
    displacement.normalize(); 


    
}
