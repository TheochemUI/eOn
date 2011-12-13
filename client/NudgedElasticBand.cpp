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

NudgedElasticBand::NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed)
{
    parameters = parametersPassed;
    images = parameters->nebImages;
    atoms = initialPassed->numberOfAtoms();
    image = new Matter *[images+2];
    tangent = new AtomMatrix *[images+2];
//    AtomMatrix tangent[images+2];

    cout <<"NEB: initialize\n";
    for(long i=0; i<images+2; i++){
        image[i] = new Matter(parameters);
        *image[i] = *initialPassed;
        tangent[i] = new AtomMatrix;
        tangent[i]->resize(atoms,3);
    }
    *image[images+1] = *finalPassed;  // final image
//    assert(atoms == matterFinal->numberOfAtoms());

    AtomMatrix posInitial = image[0]->getPositions();
    AtomMatrix posFinal = image[images+1]->getPositions();
    AtomMatrix imageSep = image[0]->pbc(posFinal-posInitial)/(images+1);
    for(long i=1; i<=images; i++) {
        image[i]->setPositions(posInitial+imageSep*double(i));
    }

    // Make sure that the endpoints know their energy
    image[0]->getPotentialEnergy();
    image[images+1]->getPotentialEnergy();
    climbingImage = 0;

    cout <<"energy 0: "<<image[0]->getPotentialEnergy()<<endl;
    cout <<"energy I: "<<image[images+1]->getPotentialEnergy()<<endl;
    return;
}

NudgedElasticBand::~NudgedElasticBand()
{
    clean();
    return;
}

void NudgedElasticBand::clean(void)
{
    for(long i=0; i<images+2; i++) {
        delete image[i];
        delete tangent[i];
    }
    delete [] image;
    delete [] tangent;
    return;
}

int NudgedElasticBand::compute(void)
{
    int status = 0;
    long iterations = 0;
    // optimizers
    //Quickmin *qm[images+2];
    //ConjugateGradients *cg[images+2];
    // temporary variables for cg
    AtomMatrix forcesStep;
    AtomMatrix posStep;
    AtomMatrix forces[images+2];
    AtomMatrix pos[images+2];

    cout <<"NEB: compute\n";
    cout <<"NEB: update forces, begin\n";
    updateForces();
    cout <<"NEB: update forces, end\n";

    // Need to generalize this, but use QM for now
    if( parameters->nebOptMethod == OPT_QM )
    {
        for(long i=1; i<images+1; i++) {
    //        qm[i] = new Quickmin(image[i], parameters); 
        }
    }
    else if ( parameters->nebOptMethod == OPT_CG ) 
    {
        for(long i=1; i<images+1; i++) {
            pos[i] = image[i]->getPositions();
            forces[i] = image[i]->getForces();
     //       cg[i] = new ConjugateGradients(image[i], parameters, forces[i]);
        }
    }

    while( convergenceForce() > parameters->optConvergedForce && iterations < parameters->optMaxIterations )
    {

        if( parameters->nebOptMethod == OPT_QM ) // Quickmin
        {
            for(long i=1; i<images+1; i++) {
      //          qm[i]->oneStep();
            }
            updateForces();
        }
        else if( parameters->nebOptMethod == OPT_CG ) // Conjugate gradients
        {
            for(long i=1; i<images+1; i++)
            {
       //         posStep = cg[i]->makeInfinitesimalStepModifiedForces(pos[i]);
                image[i]->setPositions(posStep);
            }
            updateForces();
            for(long i=1; i<images+1; i++)
            {
                forcesStep = image[i]->getForces();
       //         pos[i] = cg[i]->getNewPosModifiedForces(pos[i], forces[i], forcesStep, parameters->optMaxMove);
       //         posStep = cg[i]->makeInfinitesimalStepModifiedForces(pos[i]);
                image[i]->setPositions(posStep);
            }
            updateForces();
            for(long i=1; i<images+1; i++)
            {
                forces[i] = image[i]->getForces();
            }
        }
        iterations++;
    }

    // Cleanup
//    if( parameters->nebOptMethod == OPT_QM ) {
//        for(long i=1; i<images+1; i++) {
//            delete qm[i]; }
//     }else if( parameters->nebOptMethod == OPT_CG ) {
//        for(long i=1; i<images+1; i++) {
//            delete cg[i]; }
//    }

    return status;
}

// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void)
{
    if( parameters->nebClimbingImageMethod && climbingImage!=0 ) {
        return image[climbingImage]->getForces().norm();
    }
    double fmax = image[1]->getForces().norm();
    for(long i=2; i<images+1; i++) {
        if( fmax < image[i]->getForces().norm() ) {
            fmax = image[i]->getForces().norm();
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
    AtomMatrix force(3,atoms), forcePerp(3,atoms), forcePar(3,atoms);
    AtomMatrix forceSpringPar(3,atoms);
    AtomMatrix pos(3,atoms), posNext(3,atoms), posPrev(3,atoms);
    double distNext, distPrev;

    // update the forces on the images and find the highest energy image
    cout <<"before energy\n";
//    cout <<image[0]->getPositions();
    maxEnergy = image[0]->getPotentialEnergy();
//    cout <<"energy 0: "<<maxEnergy<<endl;
    maxEnergyImage = 0;
    for(long i=1; i<=images+1; i++) {
        image[i]->getForces();
//        cout <<"energy "<<i<<": "<<image[i]->getPotentialEnergy()<<endl;
        if(image[i]->getPotentialEnergy() > maxEnergy) {
            maxEnergy = image[i]->getPotentialEnergy();
            maxEnergyImage = i;
        }
    }
    cout <<"max energy image "<<maxEnergyImage<<endl;

    for(long i=1; i<=images; i++)
    {
        // set local variables
        force = image[i]->getForces();
        pos = image[i]->getPositions();
        posPrev = image[i-1]->getPositions();
        posNext = image[i+1]->getPositions();
        energy = image[i]->getPotentialEnergy();
        energyPrev = image[i-1]->getPotentialEnergy();
        energyNext = image[i+1]->getPotentialEnergy();

        // determine the tangent
        if(parameters->nebOldTangent) {
            // old tangent
            *tangent[i] = posNext - posPrev;
        }else{
            // new tangent
            higherEnergyPrev = energyPrev > energyNext;
            higherEnergyNext = energyNext > energyPrev;

            if(higherEnergyPrev != higherEnergyNext) {
                // we are not at an extremum
                if(higherEnergyPrev) {
                    *tangent[i] = pos - posPrev;
                }else{
                    *tangent[i] = posNext - pos;
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
                    *tangent[i] = (posNext - pos) * minDiffEnergy;
                    *tangent[i] += (pos - posPrev) * maxDiffEnergy;
                }else{
                    *tangent[i] = (posNext - pos) * maxDiffEnergy;
                    *tangent[i] += (pos - posPrev) * minDiffEnergy;
                }
            }
        }
        tangent[i]->normalize();

        // project the forces and add springs
        force = image[i]->getForces();

        if(parameters->nebClimbingImageMethod && i==maxEnergyImage)
        {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            image[i]->setForces( force - 2.0 * (force.cwise() * *tangent[i]).sum() * *tangent[i]);
        }
        else  // all non-climbing images
        {
            // calculate the force perpendicular to the tangent
            forcePerp = force - (force.cwise() * *tangent[i]).sum() * *tangent[i];

            // calculate the spring force
            distPrev = (posPrev - pos).squaredNorm();
            distNext = (posNext - pos).squaredNorm();
            forceSpringPar = parameters->nebSpring * (distNext-distPrev) * *tangent[i];

            // sum the spring and potential forces for the neb force
            image[i]->setForces( forceSpringPar + forcePerp );
        }
    }
    return;
}

class NEBObjectiveFunction : public ObjectiveFunction
{
    public:
        NEBObjectiveFunction(Matter *matterPassed,
                                Parameters *parametersPassed)
        {
            matter = matterPassed;
            parameters = parametersPassed;
        }
        ~NEBObjectiveFunction(void){};
        double getEnergy() { return matter->getPotentialEnergy(); }
        VectorXd getGradient(bool fdstep=false) { return -matter->getForcesFreeV(); }
        void setPositions(VectorXd x) { matter->setPositionsFreeV(x); }
        VectorXd getPositions() { return matter->getPositionsFreeV(); }
        int degreesOfFreedom() { return 3*matter->numberOfFreeAtoms(); }
        bool isConverged() { return getConvergence() < parameters->optConvergedForce; }
        double getConvergence() { return matter->maxForce(); }
    private:
        Matter *matter;
        Parameters *parameters;
};


