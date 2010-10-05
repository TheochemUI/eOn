#include "Mdynamics.h"

using namespace helper_functions;
using namespace constants;

Mdynamics::Mdynamics(Matter *matter,Parameters *parameters)
{
    matter_ = matter;    
    parameters_ = parameters;
    dtScale_ = 1.0; //in unit of 10fs
    
	nFreeCoord_ = 3*matter->numberOfFreeAtoms();
    tempListDouble_ = new double[nFreeCoord_];

};


Mdynamics::~Mdynamics()
{
    delete [] tempListDouble_;
    return;
};


void Mdynamics::oneStep()
{
	//Andersen(); //Wait to be implemented later;
    VerletStep1();
    VerletStep2();
    return;
};


void Mdynamics::VerletStep1()
{
    double *positions;
    double *velocities;
	double *accelerations;
	
    positions = new double[nFreeCoord_];
    velocities = new double[nFreeCoord_];
	accelerations= new double[nFreeCoord_];

	matter_->getFreePositions(positions);
    matter_->getFreeVelocities(velocities);
    matter_->getFreeAccelerations(accelerations);
     
	/* //test acc
    for (long int i = 0;i<matter->numberOfFreeAtoms();i++){
	     for (int j = 0; j < 3; j++){
     	    printf("acc for atom%ld, axis %d is %lf\n",i,j,accelerations[3*i+j]);    
		 }
    }
	*/
    multiplyScalar(tempListDouble_, accelerations, 0.5 * parameters_->mdTimeStep * dtScale_, nFreeCoord_);
    add(velocities, tempListDouble_, velocities, nFreeCoord_);
    matter_->setFreeVelocities(velocities);  // First update velocities

	/* //test vel
	for (long int i = 0;i<nFreeCoord_/3;i++){
		 for (int j = 0; j < 3; j++){
			 printf("Velcocity 2 for atom %ld,axis %d is %lf \n",i,j,velocity[3*i+j]);
		 }
	}
    */
    multiplyScalar(tempListDouble_, velocities, parameters_->mdTimeStep * dtScale_, nFreeCoord_);
    add(positions, tempListDouble_, positions, nFreeCoord_);
    matter_->setFreePositions(positions); // Update Positions
    delete [] positions;
    delete [] velocities;
	delete [] accelerations;
    return;
};

void Mdynamics::VerletStep2()
{
    double *velocities;
	double *accelerations;

    velocities = new double[nFreeCoord_];
	accelerations= new double[nFreeCoord_];

    matter_->getFreeVelocities(velocities);
    matter_->getFreeAccelerations(accelerations);	

    multiplyScalar(tempListDouble_, accelerations, 0.5 * parameters_->mdTimeStep * dtScale_, nFreeCoord_);
    add(velocities, tempListDouble_, velocities, nFreeCoord_);
    matter_->setFreeVelocities(velocities);// Second update Velocities
    delete [] velocities;
	delete [] accelerations;
    return;
};

void Mdynamics::fullSteps()
{
    bool stoped = false;
    long forceCallsTemp;
	long nsteps=0;
    forceCallsTemp = matter_->getForceCalls();  
	printf("test:steps=%ld\n",parameters_-> mdSteps);
    while(!stoped)
    {
        oneStep();
		nsteps++;
        printf("MDsteps %ld\n",nsteps);
		if (nsteps >= parameters_->mdSteps ){
	       stoped = true;
		}       
    }
    return;
};


