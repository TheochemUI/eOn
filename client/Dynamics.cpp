#include "Dynamics.h"
#include <math.h>

using namespace helper_functions;

Dynamics::Dynamics(Matter *matter_passed,Parameters *parameters_passed)
{
    matter = matter_passed;    
    parameters = parameters_passed;
    dtScale = 1.0; //in unit of 10fs
    kb=1.0/11604.5; //Kb in unit of eV
	nAtoms = matter->numberOfAtoms();

};


Dynamics::~Dynamics()
{
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
     Matrix<double, Eigen::Dynamic, 3> positions;
     Matrix<double, Eigen::Dynamic, 3> velocities;
     Matrix<double, Eigen::Dynamic, 3> accelerations;

     positions = matter->getPositions();
     velocities = matter->getVelocities();
     accelerations = matter->getAccelerations();
     
     velocities += accelerations * 0.5 * parameters->mdTimeStep * dtScale;
     matter->setVelocities(velocities);  // First update velocities

     positions += velocities * parameters->mdTimeStep * dtScale;
     matter->setPositions(positions); // Update Positions
};

void Dynamics::VerletStep2()
{
     Matrix<double, Eigen::Dynamic, 3> velocities(nAtoms,3);
     Matrix<double, Eigen::Dynamic, 3> accelerations(nAtoms,3);

     velocities = matter->getVelocities();
     accelerations = matter->getAccelerations();	

     velocities += accelerations * 0.5 * parameters->mdTimeStep * dtScale;
     matter->setVelocities(velocities);// Second update Velocities
};

void Dynamics::fullSteps()
{
     bool stoped = false;
     long forceCallsTemp;
     long nsteps=0;
     Matrix<double, Eigen::Dynamic, 3> velocity; 
     double EKin;
     double TKin,SumT=0.0,SumT2=0.0,AvgT,VarT;
  
     long nFreeCoord = matter->numberOfFreeAtoms()*3; 
     forceCallsTemp = matter->getForceCalls();  

     velocityScale();

     while(!stoped){
        oneStep();
	nsteps++;
     
	velocity = matter->getVelocities();
    EKin=matter->getKineticEnergy();
    TKin=(2*EKin/nFreeCoord/kb); 
	SumT+=TKin;
	SumT2+=TKin*TKin;
        //printf("MDsteps %ld Ekin = %lf Tkin = %lf \n",nsteps,EKin,TKin); 

	if (nsteps % 100 == 0){
	    matter->matter2xyz("movie", true);
	}

	if (nsteps >= parameters->mdSteps ){
	    stoped = true;
	}       
     }
     
     AvgT=SumT/nsteps;
     VarT=SumT2/nsteps-AvgT*AvgT;
     printf("Temperature : Average = %lf ; Variance = %lf ; Factor = %lf \n", AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord/2);
};

void Dynamics::Andersen(){
    
     double temp,alpha,Tcol,Pcol;//temp,sigma,Tcol,Pcol should be got from parameter.dat.
     double irand,v1,new_v,old_v;
     Matrix<double, Eigen::Dynamic, 1> mass;
     Matrix<double, Eigen::Dynamic, 3> velocity;

     temp = parameters->mdTemperature; //unit K
     alpha = parameters->Andersen_Alpha; //collision strength
     Tcol = parameters->Andersen_Tcol; // Average time between collision, in unit of dt
     Pcol = 1.0-exp(-1.0/Tcol);

     velocity = matter->getVelocities();
     mass = matter->getMasses();

     for (long int i = 0;i<nAtoms;i++)
     {
        if(!matter->getFixed(i))
        {
            for (int j = 0; j < 3; j++)
            {
                old_v = velocity(i,j);
                irand = randomDouble();
                if( irand < Pcol)
                {
                    v1 = sqrt(kb*temp/mass[i])*guaRandom(0.0,1.0);
                    new_v = sqrt(1-alpha*alpha)*old_v+alpha*v1;
                    velocity(i,j) = new_v;
                }
	        }
        }
     }
 
     matter->setVelocities(velocity);
     
     return;
}

void Dynamics::velocityScale(){
	
     double temp,new_v;
     Matrix<double, Eigen::Dynamic, 1> mass;
     Matrix<double, Eigen::Dynamic, 3> velocity;
     long int nFreeAtoms;

     temp = parameters->mdTemperature;
     nFreeAtoms = matter->numberOfFreeAtoms(); 
     velocity = matter->getVelocities();
     mass = matter->getMasses();

     for (long int i = 0;i<nFreeAtoms;i++)
     {
	    for (int j = 0; j < 3; j++)
        {
	       new_v = sqrt(kb*temp/mass[i])*guaRandom(0.0,1.0);
	       velocity(i,j) = new_v;
	    }
     }

     matter->setVelocities(velocity);
     
     return;
}

