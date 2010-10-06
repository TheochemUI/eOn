#include "Dynamics.h"
#include <math.h>

using namespace helper_functions;

Dynamics::Dynamics(Matter *matter,Parameters *parameters)
{
    matter_ = matter;    
    parameters_ = parameters;
    dtScale_ = 1.0; //in unit of 10fs
    kb=1.0/11604.5; //Kb in unit of eV
	nFreeCoord_ = 3*matter->numberOfFreeAtoms();
    tempListDouble_ = new double[nFreeCoord_];

};


Dynamics::~Dynamics()
{
    delete [] tempListDouble_;
    return;
};


void Dynamics::oneStep()
{
	Andersen(); //Wait to be implemented later;
    VerletStep1();
    VerletStep2();
    return;
};


void Dynamics::VerletStep1()
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

void Dynamics::VerletStep2()
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

void Dynamics::fullSteps()
{
    bool stoped = false;
    long forceCallsTemp;
	long nsteps=0;
    double *freeVelocity;
    double EKin;
	double TKin,SumT=0.0,SumT2=0.0,AvgT,VarT;
   
    freeVelocity = new double[nFreeCoord_];
    forceCallsTemp = matter_->getForceCalls();  
	printf("test:steps=%ld\n",parameters_-> mdSteps);
    
	velocityScale();

    while(!stoped)
    {
        oneStep();
		nsteps++;
        
		matter_->getFreeVelocities(freeVelocity);
        EKin=matter_->kineticEnergy();
        TKin=(2*EKin/nFreeCoord_/kb); 
		SumT+=TKin;
		SumT2+=TKin*TKin;
        printf("MDsteps %ld Ekin = %lf Tkin = %lf \n",nsteps,EKin,TKin); 
		/*//test
		for (long int i = 0;i<matter_->numberOfFreeAtoms();i++){
	 			printf("%lf	%lf	%lf \n",freeVelocity[3*i],freeVelocity[3*i+1],freeVelocity[3*i+2]);
		}
        */
		if (nsteps >= parameters_->mdSteps ){
	       stoped = true;
		}       
    }
     
	AvgT=SumT/nsteps;
	VarT=SumT2/nsteps-AvgT*AvgT;
	printf("Tempeture : Average = %lf ; Variance = %lf ; Factor = %lf \n", AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord_/2);
	

    delete [] freeVelocity;
    return;
};

void Dynamics::Andersen()
{    
	 double temp,alpha,Tcol,Pcol;//temp,sigma,Tcol,Pcol should be got from parameter.dat.
	 double irand,v1,new_v,old_v;
	 double *mass;
	 double *freeVelocity;
	 long int nFreeAtoms;

	 temp = parameters_->mdTemperture; //unit K
	 alpha = parameters_->Andersen_Alpha; //collision strength
	 Tcol = parameters_->Andersen_Tcol; // Average time between collision, in unit of dt
	 Pcol = 1.0-exp(-1.0/Tcol);
     nFreeAtoms = matter_->numberOfFreeAtoms(); 

	 /* //test parameters
	 printf("temp=%lf\n",temp);
	 printf("alpha=%lf\n",alpha);
	 printf("Tcol=%lf\n",Tcol);
	 printf("Pcol=%lf\n",Pcol);
	 */

	 mass = new double[nFreeAtoms];
     freeVelocity = new double[nFreeCoord_];
     matter_->getFreeVelocities(freeVelocity);
	 matter_->getFreeMasses(mass);

	 for (long int i = 0;i<nFreeAtoms;i++){
		 for (int j = 0; j < 3; j++){
	 	    old_v = freeVelocity[3*i+j];
		    irand = randomDouble();
			//printf("irand=%lf\n",irand);
		    if( irand < Pcol){
			   v1 = sqrt(kb*temp/mass[i])*guaRandom(0.0,1.0);
			   new_v = sqrt(1-alpha*alpha)*old_v+alpha*v1;
		       freeVelocity[3*i+j] = new_v;
			}
		 }
	 }
     /* //TEST
	 for (long int i = 0;i<nFreeAtoms;i++){
		for (int j = 0; j < 3; j++){
		    printf("Andersen Velocity for atom %ld,axis %d is %lf \n",i,j,freeVelocity[3*i+j]);
		 }
	 }
     */
	 matter_->setFreeVelocities(freeVelocity);
     
     delete [] freeVelocity;
	 delete [] mass;
	 return;
}

void Dynamics::velocityScale(){
	double temp,new_v;
	double *mass;
	double *freeVelocity;
	long int nFreeAtoms;

	temp = parameters_->mdTemperture;
    nFreeAtoms = matter_->numberOfFreeAtoms(); 
	mass = new double[nFreeAtoms];
	freeVelocity = new double[nFreeCoord_];
	matter_->getFreeVelocities(freeVelocity);
	matter_->getFreeMasses(mass);

   	for (long int i = 0;i<nFreeAtoms;i++){
		for (int j = 0; j < 3; j++){
			new_v = sqrt(kb*temp/mass[i])*guaRandom(0.0,1.0);
			freeVelocity[3*i+j] = new_v;
		}
	}

	matter_->setFreeVelocities(freeVelocity);
     
	delete [] freeVelocity;
	delete [] mass;
	return;

}




