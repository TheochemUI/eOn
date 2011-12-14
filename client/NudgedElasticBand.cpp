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
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Log.h"

using namespace helper_functions;

class NEBObjectiveFunction : public ObjectiveFunction
{
    public:

        NEBObjectiveFunction(NudgedElasticBand *nebPassed, Parameters *parametersPassed)
        {
            neb = nebPassed;
            parameters = parametersPassed;
        }

        ~NEBObjectiveFunction(void){};

        VectorXd getGradient(bool fdstep=false)
        {
            VectorXd forceV;
            forceV.resize(3*neb->atoms*neb->images);
            if(neb->movedAfterForceCall) neb->updateForces();
            for(long i=1; i<=neb->images; i++){
               forceV.segment(3*neb->atoms*(i-1),3*neb->atoms) = VectorXd::Map(neb->image[i]->getForces().data(), 3*neb->atoms); 
            }
            return -forceV;
        }

        double getEnergy()
        {
            double Energy=0;
            for(long i=1; i<=neb->images; i++) {
                Energy += neb->image[i]->getPotentialEnergy();
            }
            return Energy;
        }

        void setPositions(VectorXd x)
        {
            neb->movedAfterForceCall = true;
            for(long i=1; i<=neb->images; i++) {
                neb->image[i]->setPositions(MatrixXd::Map(x.segment(3*neb->atoms*(i-1),3*neb->atoms).data(),neb->atoms,3));
            }
        }

        VectorXd getPositions()
        {
            VectorXd posV;
            posV.resize(3*neb->atoms*neb->images);
            for(long i=1; i<=neb->images; i++){
               posV.segment(3*neb->atoms*(i-1),3*neb->atoms) = VectorXd::Map(neb->image[i]->getPositions().data(), 3*neb->atoms);
            }
            return posV;
        }

        int degreesOfFreedom() { return 3*neb->images*neb->atoms; }

        bool isConverged() { return getConvergence() < parameters->optConvergedForce; }

        double getConvergence() { return neb->convergenceForce(); }

    private:
        NudgedElasticBand *neb;
        Parameters *parameters;
};


NudgedElasticBand::NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed)
{
    parameters = parametersPassed;
    images = parameters->nebImages;
    atoms = initialPassed->numberOfAtoms();
    image = new Matter *[images+2];
    tangent = new AtomMatrix *[images+2];

    cout <<"NEB: initialize\n";
    for(long i=0; i<=images+1; i++)
    {
        image[i] = new Matter(parameters);
        *image[i] = *initialPassed;
        tangent[i] = new AtomMatrix;
        tangent[i]->resize(atoms,3);
    }
    *image[images+1] = *finalPassed;  // final image

    AtomMatrix posInitial = image[0]->getPositions();
    AtomMatrix posFinal = image[images+1]->getPositions();
    AtomMatrix imageSep = image[0]->pbc(posFinal-posInitial)/(images+1);
    for(long i=1; i<=images; i++) {
        image[i]->setPositions(posInitial+imageSep*double(i));
    }

    movedAfterForceCall = true;

    // Make sure that the endpoints know their energy
    image[0]->getPotentialEnergy();
    image[images+1]->getPotentialEnergy();
    climbingImage = 0;

    return;
}

NudgedElasticBand::~NudgedElasticBand()
{
    clean();
    return;
}

void NudgedElasticBand::clean(void)
{
    for(long i=0; i<=images+1; i++) {
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
    long iteration = 0;

    log("Nudged elastic band calculation started.\n");

    updateForces();

    NEBObjectiveFunction objf(this, parameters);

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    cout <<endl;
    cout <<"iteration     force     cimage"<<endl;
    cout <<"------------------------------"<<endl;

    while (!objf.isConverged())
    {
        if(iteration) { // so that we print forces before taking an optimizer step
            if (iteration >= parameters->nebMaxIterations) {
                status = STATUS_BAD_MAX_ITERATIONS;
                break;
            }
            optimizer->step(parameters->optMaxMove);
        }
        iteration++;
        if( parameters->nebClimbingImageMethod ) {
            printf(" %7li  %10.4f     %4li\n",iteration,convergenceForce(),climbingImage);
        } else {
            printf(" %7li  %10.4f      - \n",iteration,convergenceForce());
        }
    }
    delete optimizer;
    return status;
}

// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void)
{
    if(movedAfterForceCall) updateForces();
    if( parameters->nebClimbingImageMethod && climbingImage!=0 ) {
        return image[climbingImage]->getForces().norm();
    }
    double fmax = image[1]->getForces().norm();
    for(long i=2; i<=images; i++) {
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
    maxEnergy = image[0]->getPotentialEnergy();
    maxEnergyImage = 0;
    for(long i=1; i<=images+1; i++) {
//        cout <<"image "<<i<<endl;
//        cout <<"pos: "<<endl<<image[i]->getPositions()<<endl;
//        cout <<"energy "<<image[i]->getPotentialEnergy()<<endl;
//        cout <<"forces: "<<endl<<image[i]->getForces()<<endl;
        image[i]->getForces();
        if(image[i]->getPotentialEnergy() > maxEnergy) {
            maxEnergy = image[i]->getPotentialEnergy();
            maxEnergyImage = i;
        }
    }
//    cout <<"max energy image "<<maxEnergyImage<<endl;

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
//            cout <<"tangent, image "<<i<<endl;
//            cout <<*tangent[i]<<endl;
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

            movedAfterForceCall = false;  // so that we don't repeat a force call

//            cout <<"image "<<i<<" distPrev "<<distPrev<<" distNext "<<distNext<<" force "<<image[i]->getForces().norm() <<endl;
        }
    }
    return;
}

