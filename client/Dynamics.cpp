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
    if(parameters->Thermostat == 1){
       printf("Dynamics Using Andersen Thermostat and Verlet Integration\n");
    }
    else if(parameters->Thermostat == 2){
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
    if(parameters->thermostat == ANDERSEN){
       andersen(temperature);
       andersenVerlet();
    }
    else if(parameters->thermostat == NOSE_HOVER){
       noseHoverVerlet(temperature);
    }
    return;
}

void Dynamics::andersenVerlet()
{
     Matrix<double, Eigen::Dynamic, 3> positions;
     Matrix<double, Eigen::Dynamic, 3> velocities;
     Matrix<double, Eigen::Dynamic, 3> accelerations;

     positions = matter->getPositions();
     velocities = matter->getVelocities();
     accelerations = matter->getAccelerations();

     velocities += accelerations * 0.5 * parameters->mdTimeStep * dtScale;
     matter->setVelocities(velocities);  // first update velocities

     positions += velocities * parameters->mdTimeStep * dtScale;
     matter->setPositions(positions); // update positions

    velocities = matter->getVelocities();
    accelerations = matter->getAccelerations();

    velocities += accelerations * 0.5 * parameters->mdTimeStep * dtScale;
    matter->setVelocities(velocities); // second update velocities
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

void Dynamics::noseHoverVerlet(double temperature){
    double smass = parameters->thermoNoseMass;
    Matrix<double, Eigen::Dynamic, 3> vel;
    Matrix<double, Eigen::Dynamic, 3> pos;
    Matrix<double, Eigen::Dynamic, 3> acc;
    double q1,q2,g1,g2,s,dt,dt2,dt4,dt8;
    double bolkev = 8.6173857E-5;
    double kinE,Temp;
    long nFree;
   
    dt = parameters->mdTimeStep;
    dt2 = 0.5 * dt;
    dt4 = 0.25 * dt;
    dt8 = 0.125 * dt;
    nFree = matter->numberOfFreeAtoms()*3;
    q1 = q2 = smass;
    g1 = 0.0;
    g2 = 0.0;
    Temp = bolkev*temperature; //imposed termperature

    vel = matter->getVelocities();
    pos = matter->getPositions();
    acc = matter->getAccelerations();

    if(init == true){
        init = false;
        vxi1 = 0.0;
        vxi2 = 0.0;
        xi1 = 0.0;
        xi2 = 0.0;
    }
    kinE =  matter->getKineticEnergy();

    g2 = (q1*vxi1*vxi1 - Temp);
    vxi2 += g2 * dt4;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2*kinE-nFree*Temp)/q1;
    vxi1 += g1 * dt4;
    vxi1 *= exp(-vxi2*dt8);
    xi1 += vxi1*dt2;
    xi2 += vxi2*dt2;
    s = exp(-vxi1*dt2);
    vel *= s;
    kinE *= s*s;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2*kinE-nFree*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *=  exp(-vxi2*dt8);
    g2 = (q1*vxi1*vxi1-Temp)/q2;
    vxi2 += g2*dt4;

    pos += vel * dt2;
    matter->setPositions(pos);


    acc = matter->getAccelerations();
    vel += acc * dt;
    pos += vel * dt2;
    kinE =  matter->getKineticEnergy();

    g2 = (q1*vxi1*vxi1 - Temp);
    vxi2 += g2 * dt4;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2*kinE-nFree*Temp)/q1;
    vxi1 += g1 * dt4;
    vxi1 *= exp(-vxi2*dt8);
    xi1 += vxi1*dt2;
    xi2 += vxi2*dt2;
    s = exp(-vxi1*dt2);
    vel *= s;
    kinE *= s*s;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2*kinE-nFree*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *=  exp(-vxi2*dt8);
    g2 = (q1*vxi1*vxi1-Temp)/q2;
    vxi2 += g2*dt4;

    matter->setPositions(pos);    
    matter->setVelocities(vel);
    return;
}

