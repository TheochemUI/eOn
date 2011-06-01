//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "NudgedElasticBand.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include <cassert>

using namespace helper_functions;

const char NudgedElasticBand::OPT_QM[] = "qm";
const char NudgedElasticBand::OPT_CG[] = "cg";
const char NudgedElasticBand::OPT_LBFGS[] = "lbfgs";

NudgedElasticBand::NudgedElasticBand()
{
    return;
}

NudgedElasticBand::~NudgedElasticBand()
{
    clean();
    return;
}

NudgedElasticBand::NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed)
{
    initial = initialPassed;
    final = finalPassed;
    initialize(initialPassed, finalPassed, parametersPassed);
    return;
}

void NudgedElasticBand::clean(void)
{
    for(long i=0; i<images+2; i++) {
        delete neb[i];
    }
    return;
}

void NudgedElasticBand::initialize(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed)
{
    clean();
    initial = initialPassed;
    final = finalPassed;
    parameters = parametersPassed;
    images = parameters -> nebImages;
    Matter *neb[images+2];
    AtomMatrix tangent[images+2];

    for(long i=0; i<images+2; i++){
        neb[i] = new Matter(parameters);
        *neb[i] = *initial;
        tangent[i].resize(nAtoms,3);
    }
    *neb[images+1] = *final;  // final image
    nAtoms = initial->numberOfAtoms();
    assert(nAtoms == matterFinal->numberOfAtoms());

    AtomMatrix posInitial = neb[0]->getPositions();
    AtomMatrix posFinal = neb[images+1]->getPositions();
    AtomMatrix imageSep = neb[0]->pbc(posFinal-posInitial)/(images+1);
    for(long i=1; i<=images; i++) {
        neb[i]->setPositions(posInitial+imageSep*double(i));
    }
    // Make sure that the endpoints know their energy
    neb[0]->getPotentialEnergy();
    neb[images+1]->getPotentialEnergy();
    climbingImage = 0;
}

int NudgedElasticBand::compute(void)
{
    int status = 0;
    long iterations = 0;
    // optimizers
    Quickmin *qm[images+2];
    ConjugateGradients *cg[images+2];
    // temporary variables for cg
    AtomMatrix forcesStep;
    AtomMatrix posStep;
    AtomMatrix forces[images+2];
    AtomMatrix pos[images+2];
    
    updateForces();

    // Need to generalize this, but use QM for now
    if( parameters->nebOptMethod == OPT_QM )
    {
        for(long i=1; i<images+1; i++) {
            qm[i] = new Quickmin(neb[i], parameters); 
        }
    }
    else if ( parameters->nebOptMethod == OPT_CG ) 
    {
        for(long i=1; i<images+1; i++) {
            pos[i] = neb[i]->getPositions();
            forces[i] = neb[i]->getForces();
            cg[i] = new ConjugateGradients(neb[i], parameters, forces[i]);
        }
    }

    while( convergenceForce() > parameters->optConvergedForce && iterations < parameters->optMaxIterations )
    {

        if( parameters->nebOptMethod == OPT_QM ) // Quickmin
        {
            for(long i=1; i<images+1; i++) {
                qm[i]->oneStep();
            }
            updateForces();
        }
        else if( parameters->nebOptMethod == OPT_CG ) // Conjugate gradients
        {
            for(long i=1; i<images+1; i++)
            {
                posStep = cg[i]->makeInfinitesimalStepModifiedForces(pos[i]);
                neb[i]->setPositions(posStep);
            }
            updateForces();
            for(long i=1; i<images+1; i++)
            {
                forcesStep = neb[i]->getForces();
                pos[i] = cg[i]->getNewPosModifiedForces(pos[i], forces[i], forcesStep, parameters->optMaxMove);
                posStep = cg[i]->makeInfinitesimalStepModifiedForces(pos[i]);
                neb[i]->setPositions(posStep);
            }
            updateForces();
            for(long i=1; i<images+1; i++)
            {
                forces[i] = neb[i]->getForces();
            }
        }
        iterations++;
    }

    // Cleanup
    if( parameters->nebOptMethod == OPT_QM ) {
        for(long i=1; i<images+1; i++) {
            delete qm[i]; }
     }else if( parameters->nebOptMethod == OPT_CG ) {
        for(long i=1; i<images+1; i++) {
            delete cg[i]; }
    }

    return status;
}

// generate the force value which is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void)
{
    if( parameters->nebClimbingImageMethod && climbingImage!=0 ) {
        return neb[climbingImage]->getForces().norm();
    }
    double fmax = neb[1]->getForces().norm();
    for(long i=2; i<images+1; i++) {
        if( fmax < neb[i]->getForces().norm() ) {
            fmax = neb[i]->getForces().norm();
        }
    }
    return fmax;
}


// Update the forces, do the projections, and add spring forces
void NudgedElasticBand::updateForces(void)
{
    // variables for tangent
    double maxDiffEnergy, minDiffEnergy;
    double energyDiffPrev, energyDiffNext;
    double energy, energyPrev, energyNext;
    bool higherEnergyPrev, higherEnergyNext;

    // variables for climbing image
    double maxEnergy;
    long maxEnergyImage;

    // variables for force projections
    AtomMatrix force(3,nAtoms), forcePerp(3,nAtoms), forcePar(3,nAtoms);
    AtomMatrix forceSpringPar(3,nAtoms);
    AtomMatrix pos(3,nAtoms), posNext(3,nAtoms), posPrev(3,nAtoms);
    double distNext, distPrev;

    // update the forces on the images and find the highest energy image
    maxEnergy = neb[0]->getPotentialEnergy();
    maxEnergyImage = 0;
    for(long i=1; i<=images; i++) {
        neb[i]->getForces();
        if(neb[i]->getPotentialEnergy() > maxEnergy) {
            maxEnergy = neb[i]->getPotentialEnergy();
            maxEnergyImage = i;
        }
    }

    for(long i=1; i<=images; i++)
    {
        // set local variables
        force = neb[i]->getForces();
        pos = neb[i]->getPositions();
        posPrev = neb[i-1]->getPositions();
        posNext = neb[i+1]->getPositions();
        energy = neb[i]->getPotentialEnergy();
        energyPrev = neb[i-1]->getPotentialEnergy();
        energyNext = neb[i+1]->getPotentialEnergy();

        // determine the tangent
        if(parameters->nebOldTangent) {
            // old tangent
            tangent[i] = posNext - posPrev;
        }else{
            // new tangent
            higherEnergyPrev = energyPrev > energyNext;
            higherEnergyNext = energyNext > energyPrev;

            if(higherEnergyPrev != higherEnergyNext) {
                // we are not at an extremum
                if(higherEnergyPrev) {
                    tangent[i] = pos - posPrev;
                }else{
                    tangent[i] = posNext - pos;
                }
            }else{
                // we are at an extremum
                energyDiffPrev = energyPrev - energy;
                energyDiffNext = energyNext - energy;

                // calculate the energy difference to neighboring images
                minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
                maxDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));

                // use these energy differences to weight the tangent
                if(energyDiffPrev > energyDiffNext) {
                    tangent[i] = (posNext - pos) * minDiffEnergy;
                    tangent[i] += (pos - posPrev) * maxDiffEnergy;
                }else{ 
                    tangent[i] = (posNext - pos) * maxDiffEnergy;
                    tangent[i] += (pos - posPrev) * minDiffEnergy;
                }
            }
        }
        tangent[i].normalize();

        // project the forces and add springs
        force = neb[i]->getForces();

        if(parameters->nebClimbingImageMethod && i==maxEnergyImage)
        {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            neb[i]->setForces( force - 2.0 * (force.cwise() * tangent[i]).sum() * tangent[i]);
        }
        else  // all non-climbing images
        {
            // calculate the force perpendicular to the tangent
            forcePerp = force - (force.cwise() * tangent[i]).sum() * tangent[i];

            // calculate the spring force
            distPrev = (posPrev - pos).squaredNorm();
            distNext = (posNext - pos).squaredNorm();
            forceSpringPar = parameters->nebSpring * (distNext-distPrev) * tangent[i];

            // sum the spring and potential forces for the neb force
            neb[i]->setForces( forceSpringPar + forcePerp );
        }
    }
    return;
}

