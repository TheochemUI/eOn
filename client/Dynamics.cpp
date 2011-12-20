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
#include "Log.h"
#include <math.h>

using namespace helper_functions;

const char Dynamics::ANDERSEN[] = "andersen";
const char Dynamics::NOSE_HOOVER[] = "nose_hoover";
const char Dynamics::LANGEVIN[] = "langevin";
const char Dynamics::NONE[] = "none";

Dynamics::Dynamics(Matter *matter_passed, Parameters *parameters_passed)
{
    matter = matter_passed;
    parameters = parameters_passed;
    dt = parameters->mdTimeStep;
    nAtoms = matter->numberOfAtoms();
    nFreeCoords = matter->numberOfFreeAtoms()*3;
    temperature = parameters->temperature;
    kb = 8.6173857E-5; // Kb in unit of eV
    vxi1 = vxi2 = xi1 = xi2 = 0.0; // NoseHoover variables
}

Dynamics::~Dynamics()
{
    return;
}

void Dynamics::setTemperature(double temperature_passed)
{
    temperature = temperature_passed;
}

void Dynamics::oneStep(void)
{
    if(parameters->thermostat == ANDERSEN){
       andersenCollision();
       velocityVerlet();
    }
    else if(parameters->thermostat == NOSE_HOOVER){
       noseHooverVerlet();
    }
    else if(parameters->thermostat == LANGEVIN){
       langevinVerlet();
    }
    else if(parameters->thermostat == NONE){
       velocityVerlet();
    }
}

void Dynamics::velocityVerlet()
{
    AtomMatrix positions = matter->getPositions();
    AtomMatrix velocities = matter->getVelocities();
    AtomMatrix accelerationsInitial = matter->getAccelerations();

    positions += (dt*velocities) + (0.5*dt*dt*accelerationsInitial);
    matter->setPositions(positions);

    AtomMatrix accelerationsFinal = matter->getAccelerations();

    velocities += 0.5*dt*(accelerationsInitial + accelerationsFinal);
    matter->setVelocities(velocities);
}

void Dynamics::run()
{
    AtomMatrix velocity;
    double potE, kinE, kinT;
    double sumT = 0, sumT2 = 0, avgT, varT, stdT;

    setThermalVelocity();

    if(parameters->thermostat != NONE) {
        log("[Dynamics] Running NVT molecular dynamics at %8.2lf K for %10ld steps\n", temperature, parameters->mdSteps);
    }else{
        log("[Dynamics] Running NVE molecular dynamics for %10ld steps\n", parameters->mdSteps);
    }

    if (parameters->writeMovies == true) {
        matter->matter2con("dynamics", false);
    }

    log("[Dynamics] %8s %10s %12s %12s %10s\n", "Step", "KE", "PE", "TE", "kinT");

    for(long step=0; step<=parameters->mdSteps; step++)
    {
        oneStep();

        velocity = matter->getVelocities();
        kinE = matter->getKineticEnergy();
        potE = matter->getPotentialEnergy();
        kinT = (2.0*kinE/nFreeCoords/kb);
        sumT += kinT;
        sumT2 += kinT*kinT;

        log("[Dynamics] %8ld %10.4f %12.4f %12.4f %10.2f\n", step, kinE, potE, kinE+potE, kinT);

        if ( (parameters->writeMovies == true) && (step % parameters->writeMoviesInterval == 0) ) {
            matter->matter2con("dynamics", true);
        }
    }
    avgT = sumT/double(parameters->mdSteps);
    varT = sumT2/double(parameters->mdSteps) - avgT*avgT;
    stdT = sqrt(varT);
    log("[Dynamics] Temperature : Average = %.2lf ; StdDev = %.2lf ; Factor = %.2lf\n",
        avgT, stdT, varT/avgT/avgT*nFreeCoords/2.0);
}

void Dynamics::andersenCollision()
{
    double alpha, tCol, pCol;
    double vNew, vOld;
    Matrix<double, Eigen::Dynamic, 1> mass;
    AtomMatrix velocity;

    alpha = parameters->thermoAndersenAlpha; // collision strength
    tCol = parameters->thermoAndersenTcol; // average time between collisions, in unit of fs
    pCol = 1.0-exp(-dt/tCol);

    velocity = matter->getVelocities();
    mass = matter->getMasses();

    for (long i=0; i<nAtoms; i++)
    {
        if( (randomDouble() < pCol) && (!matter->getFixed(i)) )
        {
            for (int j=0; j<3; j++)
            {
                vOld = velocity(i,j);
                vNew = sqrt(kb*temperature/mass[i])*gaussRandom(0.0,1.0);
                velocity(i,j) = sqrt(1.0-alpha*alpha)*vOld + alpha*vNew;
            }
        }
    }
    matter->setVelocities(velocity);
}

void Dynamics::setThermalVelocity()
{
    AtomMatrix velocity = matter->getVelocities();
    Matrix<double, Eigen::Dynamic, 1> mass = matter->getMasses();

    for (long i=0; i<nAtoms; i++)
    {
        if(!matter->getFixed(i)){
            for (int j=0; j<3; j++)
            {
                velocity(i,j) = sqrt(kb*temperature/mass[i])*gaussRandom(0.0,1.0);
            }
        }
    }
    matter->setVelocities(velocity);
}

void Dynamics::rescaleVelocity()
{
    AtomMatrix velocity = matter->getVelocities();
    double kinE = matter->getKineticEnergy();
    double kinT = (2.0*kinE/nFreeCoords/kb);
    matter->setVelocities(velocity*sqrt(temperature/kinT));
}

void Dynamics::noseHooverVerlet()
{
    AtomMatrix vel, pos, acc;
    double q1, q2, g1, g2, s, dt2, dt4, dt8;
    double kinE, Temp;

    dt2 = 0.5*dt;
    dt4 = 0.25*dt;
    dt8 = 0.125*dt;
    q1 = q2 = parameters->thermoNoseMass;
    g1 = g2 = 0.0;
    Temp = kb*temperature; // imposed temperature

    vel = matter->getVelocities();
    pos = matter->getPositions();
    acc = matter->getAccelerations();

    kinE =  matter->getKineticEnergy();

    g2 = (q1*vxi1*vxi1-Temp);
    vxi2 += g2*dt4;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFreeCoords*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *= exp(-vxi2*dt8);
    xi1 += vxi1*dt2;
    xi2 += vxi2*dt2;
    s = exp(-vxi1*dt2);
    vel *= s;
    kinE *= s*s;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFreeCoords*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *=  exp(-vxi2*dt8);
    g2 = (q1*vxi1*vxi1-Temp)/q2;
    vxi2 += g2*dt4;

    pos += vel*dt2;
    matter->setPositions(pos);

    acc = matter->getAccelerations();
    vel += acc*dt;
    pos += vel*dt2;
    kinE =  matter->getKineticEnergy();

    g2 = (q1*vxi1*vxi1-Temp);
    vxi2 += g2*dt4;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFreeCoords*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *= exp(-vxi2*dt8);
    xi1 += vxi1*dt2;
    xi2 += vxi2*dt2;
    s = exp(-vxi1*dt2);
    vel *= s;
    kinE *= s*s;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFreeCoords*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *=  exp(-vxi2*dt8);
    g2 = (q1*vxi1*vxi1-Temp)/q2;
    vxi2 += g2*dt4;

    matter->setPositions(pos);
    matter->setVelocities(vel);
}

void Dynamics::langevinVerlet()
{
    Matrix<double, Eigen::Dynamic, 1> mass;
    AtomMatrix pos;
    AtomMatrix vel;
    AtomMatrix acc;
    AtomMatrix friction;
    AtomMatrix noise;
    double gamma;

    gamma = parameters->thermoLangvinFriction;
    pos = matter->getPositions();
    vel = matter->getVelocities();

    acc = matter->getAccelerations();
    noise = acc;
    mass = matter->getMasses();

    friction = - gamma * vel;
    for (long i=0; i<nAtoms; i++){
        if(!matter->getFixed(i)){
            for (int j=0; j<3; j++){
                noise(i,j) = sqrt(4.0*gamma*kb*temperature/dt/mass[i])*gaussRandom(0.0,1.0);
            }
        }
    }
    acc = acc + friction + noise;

    vel += acc*0.5*dt; // calculate velocites v(n+1/2)
    pos += vel*dt;
    matter->setPositions(pos); // update positions x(n+1)

    acc = matter->getAccelerations();
    friction = - gamma*vel;
    for (long i=0; i<nAtoms; i++){
        if(!matter->getFixed(i)){
            for (int j=0; j<3; j++){
                noise(i,j) = sqrt(4.0*gamma*kb*temperature/dt/mass[i])*gaussRandom(0.0,1.0);
            }
        }
    }
    acc = acc + friction + noise;
    vel += 0.5*dt*acc;
    matter->setVelocities(vel); // calculate velocities v(n+1)
}

