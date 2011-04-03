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
#include <cassert>

using namespace helper_functions;

NEB::NEB(Matter const *matterInitial, Matter const *matterFinal, Parameters *params)
{
    parameters = params;
    images = parameters -> nebImages;
    Matter *neb[images+2];
    AtomMatrix tangent[images+2];

    for(long i=0; i<images+2; i++){
        neb[i] = new Matter(parameters);
        *neb[i] = *matterInitial;
        tangent[i].resize(nAtoms,3);
    }
    *neb[images+1] = *matterFinal;  // final image
    nAtoms = matterInitial->numberOfAtoms();
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
}

NEB::~NEB()
{
    for(long i=0; i<images+2; i++) {
        delete neb[i];
    }
}

void NEB::compute(void)
{
    return;
}

void NEB::updateForces(void)
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

        if(parameters->nebClimbingImage && i==maxEnergyImage)
        {
            // we are at the climbing image
            neb[i]->setForces( force - 2.0*(force.cwise()*tangent[i]).sum() * tangent[i]);
        }
        else  // all non-climbing images
        {
            // calculate the force perpendicular to the tangent
            forcePerp = force - (force.cwise()*tangent[i]).sum() * tangent[i];

            // calculate the spring force
            distPrev = (posPrev - pos).squaredNorm();
            distNext = (posNext - pos).squaredNorm();
            forceSpringPar = parameters->nebSpring*(distNext-distPrev) * tangent[i];

            // sum the spring and potential forces for the neb force
            neb[i]->setForces( forceSpringPar + forcePerp );
        }
    }
    return;
}

