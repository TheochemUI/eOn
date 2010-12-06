//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Dynamics.h"
#include <math.h>

using namespace helper_functions;

Dynamics::Dynamics(Matter *matter_passed,Parameters *parameters_passed)
{
    matter = matter_passed;    
    parameters = parameters_passed;
    dtScale = 1.0; // in unit of 10fs
    kb = 1.0/11604.5; // Kb in unit of eV
    nAtoms = matter->numberOfAtoms();
    init = true;
/*
    if(parameters->ThermoType == 1){
       printf("Dynamics Using Andersen Thermostat and Verlet Integration\n");
    }
    else if(parameters->ThermoType == 2){
       printf("Dynamics Using Nose-Hoover Thermostat and Verlet Integration\n");
    }
*/
}

Dynamics::~Dynamics()
{
    return;
}

void Dynamics::oneStep(double temperature)
{
    if(parameters->thermoType == ANDERSEN){
       andersenVerlet(temperature);
    }
    else if(parameters->thermoType == NOSE_HOVER){
       noseHoverVerlet(temperature);
    }
    return;
}

void Dynamics::andersenVerlet(double temperature)
{
     Matrix<double, Eigen::Dynamic, 3> positions;
     Matrix<double, Eigen::Dynamic, 3> velocities;
     Matrix<double, Eigen::Dynamic, 3> accelerations;
     
     positions = matter->getPositions();
     velocities = matter->getVelocities();
     accelerations = matter->getAccelerations();
     
     andersen(temperature);//andersen Thermostat;

     velocities += accelerations * 0.5 * parameters->mdTimeStep * dtScale;
     matter->setVelocities(velocities);  // first update velocities

     positions += velocities * parameters->mdTimeStep * dtScale;
     matter->setPositions(positions); // update positions

     accelerations = matter->getAccelerations();
     velocities += accelerations * 0.5 * parameters->mdTimeStep * dtScale;
     matter->setVelocities(velocities); // second update velocities
}


void Dynamics::noseHoverVerlet(double temperature){
    double smass = parameters->thermoNoseMass;
    Matrix<double, Eigen::Dynamic, 3> v;
    Matrix<double, Eigen::Dynamic, 3> p;
    Matrix<double, Eigen::Dynamic, 3> d2;
    Matrix<double, Eigen::Dynamic, 3> d3;
    Matrix<double, Eigen::Dynamic, 3> d2c;
    double s1,s2,s3,s4,fact,sqq,svel,sc;
    double ev2j = 1.60217733E-19, am2kg = 1.6605402E-27, bolkev = 8.6173857E-5;
    double kinE,error,tmp,eps,es;
    long nFree,i,j;
    bool stop;
    double ul=1.0E-10, ut = parameters->mdTimeStep*dtScale*1.0E-14;
//    printf("ut = %e\n",ut);

    fact = (am2kg/ev2j)*(ul/ut)*(ul/ut);
    nFree = nAtoms*3.0;
    eps = 0.0;
    es = 0.0;
    v = matter->getVelocities();
    p = matter->getPositions();
    d2c = matter->getAccelerations();

    if(smass > 0){
        sqq = smass * fact;
    }else{
        sqq = 1.0;
    }

    if(init == true){
        init = false;
        kinE = matter->getKineticEnergy();
        s1 = 1.0;
        s2 = 0.0;
        s3 = 0.0;
        if(smass > 0){
            s3 = (kinE-nFree*bolkev*temperature/2.0)*s1/sqq;
        }
        s4 = 0.0;
    }

    d2 = v;
    svel = s2;
    long Iter = 0;
    stop = false;
    error = 0.0;
    while(!stop){
        if(smass > 0){
            sc = svel/s1;
        }else{
            sc = 0.0;
        }

        kinE = matter->getKineticEnergy();
        for(i=0;i<nAtoms;i++){
            for(j=0;j<3;j++){
                error += (v(i,j)+d2c(i,j)-sc*0.5*d2(i,j)-d2(i,j))*(v(i,j)+d2c(i,j)-sc*0.5*d2(i,j)-d2(i,j));
                d2(i,j) = v(i,j)+d2c(i,j)-sc*0.5*d2(i,j);
            }
        }
 
        if(smass > 0){
            tmp = (s2+s1*((kinE-nFree*bolkev*temperature/2.0)/sqq+0.5*sc*sc)-svel);
            tmp = tmp * tmp;
            error += tmp;
            svel = s2+s1*((kinE-nFree*bolkev*temperature/2.0)/sqq+0.5*sc*sc);
        }

        Iter++;
        if(Iter >= 10){ stop = true; }
        else if(sqrt(error)<1E-10){ stop = true; }
    }
                 
//FINAL VERLET STEP
    if (smass > 0){
        sc = svel/s1;
        eps = 0.5*sqq*(svel/s1)*(svel/s1);
        es = nFree*bolkev*temperature*log(s1);
    }
    else{
        sc = 0.0;
        eps = 0.0;
        es = 0.0;
    }

    v += 2.0*d2c-sc*d2;
    p += v;

    if (smass > 0){
        s2 += s1*((2*kinE-bolkev*temperature*nFree)/sqq+sc*sc);
        s1 += s2;
    }

    matter->setVelocities(v);
    matter->setPositions(p);
    return;
}


void Dynamics::fullSteps(double temperature)
{
    bool stoped = false;
    long forceCallsTemp;
    long nsteps = 0;
    Matrix<double, Eigen::Dynamic, 3> velocity; 
    double kinE, kinT;
    double sumT=0.0, sumT2=0.0, avgT, varT;
    long nFreeCoord = matter->numberOfFreeAtoms()*3; 
    forceCallsTemp = matter->getForceCalls();  

    velocityScale(temperature);

    while(!stoped)
    {
        oneStep(temperature);
        nsteps++;

        velocity = matter->getVelocities();
        kinE = matter->getKineticEnergy();
        kinT = (2*kinE/nFreeCoord/kb); 
        sumT += kinT;
        sumT2 += kinT*kinT;
        //printf("MDsteps %ld kinE = %lf Tkin = %lf \n",nsteps,kinE,kinT); 
/*
        if (nsteps % 100 == 0){
            matter->matter2xyz("movie", true);
        }
*/
        if (nsteps >= parameters->mdSteps){
            stoped = true;
        }
    }

    avgT=sumT/nsteps;
    varT=sumT2/nsteps-avgT*avgT;
    printf("Temperature : Average = %lf ; Variance = %lf ; Factor = %lf \n", avgT,varT,varT/avgT/avgT*nFreeCoord/2);
}

void Dynamics::velocityScale(double temperature)
{
    double vNew, kinT, kinE;
    Matrix<double, Eigen::Dynamic, 1> mass;
    Matrix<double, Eigen::Dynamic, 3> velocity;
    long nFreeCoord = matter->numberOfFreeAtoms()*3;

    velocity = matter->getVelocities();
    mass = matter->getMasses();

    for (long int i = 0;i<nAtoms;i++)
    {
        if(!matter->getFixed(i))
        for (int j = 0; j < 3; j++)
        {
           // printf("mass[%ld] = %lf\n",i,mass(i));
           vNew = sqrt(kb*temperature/mass[i])*gaussRandom(0.0,1.0);
           velocity(i,j) = vNew;
        }
    }

    matter->setVelocities(velocity);  
    kinE = matter->getKineticEnergy();
    kinT = (2*kinE/nFreeCoord/kb);
//    printf("Tkin_1 = %lf\n",kinT);
    velocity = velocity*sqrt(temperature/kinT);
    matter->setVelocities(velocity);
//    kinE = matter->getKineticEnergy();
//    kinT = (2*EKin/nFreeCoord/kb);
//    printf("Tkin_2 = %lf\n",kinT);

    return;
}

void Dynamics::andersen(double temperature)
{
    double alpha, tCol, pCol; // sigma, tCol, pCol should be obtained from parameter.dat
    double irand, v1, vNew, vOld;
    Matrix<double, Eigen::Dynamic, 1> mass;
    Matrix<double, Eigen::Dynamic, 3> velocity;

    alpha = parameters->thermoAndersenAlpha; // collision strength
    tCol = parameters->thermoAndersenTcol; // average time between collision, in unit of dt
    pCol = 1.0-exp(-1.0/tCol);

    velocity = matter->getVelocities();
    mass = matter->getMasses();

    for (long int i = 0;i<nAtoms;i++)
    {
        if(!matter->getFixed(i))
        {
            for (int j = 0; j < 3; j++)
            {
                vOld = velocity(i,j);
                irand = randomDouble();
                if( irand < pCol)
                {
                    v1 = sqrt(kb*temperature/mass[i])*gaussRandom(0.0,1.0);
                    vNew = sqrt(1-alpha*alpha)*vOld+alpha*v1;
                    velocity(i,j) = vNew;
                }
            }
        }
    }
    matter->setVelocities(velocity);
    return;
}


  
