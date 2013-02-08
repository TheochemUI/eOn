#include "BasinHoppingSaddleSearch.h"
#include "Log.h"
#include "NudgedElasticBand.h"
#include "MinModeSaddleSearch.h"
#include "LowestEigenmode.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"

BasinHoppingSaddleSearch::BasinHoppingSaddleSearch(Matter *reactantPassed, Matter *displacementPassed,
                                           Parameters *parametersPassed)
{
    reactant = new Matter(parameters);
    *reactant = *reactantPassed;
    parameters = parametersPassed;

    saddle = displacementPassed;

    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
}

BasinHoppingSaddleSearch::~BasinHoppingSaddleSearch()
{
    delete reactant;
    delete product;
}

int BasinHoppingSaddleSearch::run(void)
{
    //minimize "saddle"

    //*product = *saddle;
    //accept or reject based on boltzman exp(-de/(kB*parameters->temperature))

    // NEB reactant to minimized "saddle"

    // pick the maximum energy image along the band and do dimer
    
    //eigenvalue = ;
    //eigenvector = ;

}

double BasinHoppingSaddleSearch::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix BasinHoppingSaddleSearch::getEigenvector()
{
    return eigenvector;
}
