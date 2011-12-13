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
#include "Optimizer.h"
#include "Log.h"
#include <math.h>

using namespace helper_functions;

const char Dynamics::ANDERSEN[] = "andersen";
const char Dynamics::NOSE_HOOVER[] = "nose_hoover";
const char Dynamics::LANGEVIN[] = "langevin";
const char Dynamics::NVE[] = "nve";

Dynamics::Dynamics(Matter *matter_passed,Parameters *parameters_passed)
{
    matter = matter_passed;
    parameters = parameters_passed;
    dt = parameters->mdTimeStep;
    kb = 1.0/11604.5; // Kb in unit of eV
    nAtoms = matter->numberOfAtoms();
    init = true;
    min_fcalls = 0;
    md_fcalls = 0;
    rf_fcalls = 0;
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
    else if(parameters->thermostat == NOSE_HOOVER){
       noseHooverVerlet(temperature);
    }
    else if(parameters->thermostat == LANGEVIN){
       langevinVerlet(temperature);
    }
    else if(parameters->thermostat == NVE){
       velocityVerlet();
    }

    return;
}

void Dynamics::andersenVerlet()
{
    AtomMatrix positions;
    AtomMatrix velocities;
    AtomMatrix accelerations;

    positions = matter->getPositions();
    velocities = matter->getVelocities();
    accelerations = matter->getAccelerations();
    md_fcalls ++; // This is adding 2 force calls per md step; it should be one

    velocities += 0.5*dt*accelerations;
    matter->setVelocities(velocities); // first update velocities

    positions += dt*velocities;
    matter->setPositions(positions); // update positions

    velocities = matter->getVelocities();
    accelerations = matter->getAccelerations();
    md_fcalls ++;

    velocities += 0.5*dt*accelerations;
    matter->setVelocities(velocities); // second update velocities
}

void Dynamics::velocityVerlet()
{
    AtomMatrix positions;
    AtomMatrix velocities;
    AtomMatrix accelerations0;
    AtomMatrix accelerations1;

    positions = matter->getPositions();
    velocities = matter->getVelocities();
    accelerations0 = matter->getAccelerations();
    md_fcalls++;

    positions += (dt*velocities) + (0.5*dt*dt*accelerations0);
    matter->setPositions(positions);

    accelerations1 = matter->getAccelerations();
    md_fcalls++;

    velocities += 0.5*dt*(accelerations0 + accelerations1);
    matter->setVelocities(velocities);
}

void Dynamics::fullSteps(double temperature)
{
    bool stopped = false;
    long forceCallsTemp;
    long nsteps = 0;
    long movie_steps = parameters->writeMoviesSteps;
    AtomMatrix velocity;
    double kinE, kinT;
    double sumT = 0.0, sumT2 = 0.0, avgT, varT, stdT;
    long nFreeCoord = matter->numberOfFreeAtoms()*3;
    forceCallsTemp = matter->getForceCalls();


    if(parameters->thermostat != NVE) {
        log("[Dynamics] Running NVT molecular dynamics at %8.2lf K for %10ld steps\n", temperature, parameters->mdSteps);
        initialVel(temperature);
    }else{
        log("[Dynamics] Running NVE molecular dynamics for %10ld steps\n", parameters->mdSteps);
    }

    if (parameters->writeMovies == true) {
        matter->matter2con("dynamics", false);
    }

    log("[Dynamics] %8s %10s %12s %12s %10s\n", "Step", "KE", "PE", "TE", "kinT");

    while(!stopped) {
        oneStep(temperature);
        nsteps++;

        velocity = matter->getVelocities();
        kinE = matter->getKineticEnergy();
        double PE = matter->getPotentialEnergy();
        kinT = (2.0*kinE/nFreeCoord/kb);
        sumT += kinT;
        sumT2 += kinT*kinT;

        log("[Dynamics] %8ld %10.4f %12.4f %12.4f %10.2f\n", nsteps, kinE, PE, kinE+PE, kinT);

        if (parameters->writeMovies == true) {
            if (nsteps % movie_steps == 0){
                matter->matter2con("dynamics", true);
            }
        }

        if (nsteps >= parameters->mdSteps){
            stopped = true;
        }
    }

    avgT = sumT/nsteps;
    varT = sumT2/nsteps - avgT*avgT;
    stdT = sqrt(varT);
    log("[Dynamics] Temperature : Average = %.2lf ; StdDev = %.2lf ; Factor = %.2lf\n", avgT,stdT,varT/avgT/avgT*nFreeCoord/2.0);
}

void Dynamics::andersen(double temperature)
{
    double alpha, tCol, pCol; // sigma, tCol, pCol should be obtained from parameter.dat
    double irand, v1, vNew, vOld;
    Matrix<double, Eigen::Dynamic, 1> mass;
    AtomMatrix velocity;

    alpha = parameters->thermoAndersenAlpha; // collision strength
    tCol = parameters->thermoAndersenTcol; // average time between collisions, in unit of fs
    pCol = 1.0 - exp(-dt/tCol);

    velocity = matter->getVelocities();
    mass = matter->getMasses();

    for (long int i = 0; i < nAtoms; i++)
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
                    vNew = sqrt(1.0 - alpha*alpha)*vOld + alpha*v1;
                    velocity(i,j) = vNew;
                }
            }
        }
    }
    matter->setVelocities(velocity);
    return;
}

void Dynamics::initialVel(double temperature)
{
    double vNew, kinT, kinE;
    Matrix<double, Eigen::Dynamic, 1> mass;
    AtomMatrix velocity;
    long nFreeCoord = matter->numberOfFreeAtoms()*3;

    velocity = matter->getVelocities();
    mass = matter->getMasses();

    for (long int i = 0; i < nAtoms; i++)
    {
        if(!matter->getFixed(i))
        for (int j = 0; j < 3; j++)
        {
           vNew = sqrt(kb*temperature/mass[i])*gaussRandom(0.0,1.0);
           velocity(i,j) = vNew;
        }
    }

    matter->setVelocities(velocity);  
    kinE = matter->getKineticEnergy();
    kinT = (2.0*kinE/nFreeCoord/kb);
    velocity = velocity*sqrt(temperature/kinT);
    matter->setVelocities(velocity);

    return;
}

void Dynamics::velRescaling(double temperature)
{
    double  kinT, kinE;
    Matrix<double, Eigen::Dynamic, 1> mass;
    AtomMatrix velocity;
    long nFreeCoord = matter->numberOfFreeAtoms()*3;

    velocity = matter->getVelocities();
    mass = matter->getMasses();

    matter->setVelocities(velocity);
    kinE = matter->getKineticEnergy();
    kinT = (2.0*kinE/nFreeCoord/kb);
//    printf("Tkin_1 = %lf\n",kinT);
    velocity = velocity*sqrt(temperature/kinT);
    matter->setVelocities(velocity);

    return;
}

void Dynamics::noseHooverVerlet(double temperature){
    double smass = parameters->thermoNoseMass;
    AtomMatrix vel;
    AtomMatrix pos;
    AtomMatrix acc;
    double q1, q2, g1, g2, s, dt2, dt4, dt8;
    double bolkev = 8.6173857E-5;
    double kinE, Temp;
    long nFree;

    dt2 = 0.5*dt;
    dt4 = 0.25*dt;
    dt8 = 0.125*dt;
    nFree = matter->numberOfFreeAtoms()*3;
    q1 = q2 = smass;
    g1 = 0.0;
    g2 = 0.0;
    Temp = bolkev*temperature; //imposed temperature

    vel = matter->getVelocities();
    pos = matter->getPositions();
    acc = matter->getAccelerations();
    md_fcalls ++;

    if(init == true){
        init = false;
        vxi1 = 0.0;
        vxi2 = 0.0;
        xi1 = 0.0;
        xi2 = 0.0;
    }
    kinE =  matter->getKineticEnergy();

    g2 = (q1*vxi1*vxi1-Temp);
    vxi2 += g2*dt4;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFree*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *= exp(-vxi2*dt8);
    xi1 += vxi1*dt2;
    xi2 += vxi2*dt2;
    s = exp(-vxi1*dt2);
    vel *= s;
    kinE *= s*s;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFree*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *=  exp(-vxi2*dt8);
    g2 = (q1*vxi1*vxi1-Temp)/q2;
    vxi2 += g2*dt4;

    pos += vel*dt2;
    matter->setPositions(pos);

    acc = matter->getAccelerations();
    md_fcalls ++;
    vel += acc*dt;
    pos += vel*dt2;
    kinE =  matter->getKineticEnergy();

    g2 = (q1*vxi1*vxi1-Temp);
    vxi2 += g2*dt4;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFree*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *= exp(-vxi2*dt8);
    xi1 += vxi1*dt2;
    xi2 += vxi2*dt2;
    s = exp(-vxi1*dt2);
    vel *= s;
    kinE *= s*s;
    vxi1 *= exp(-vxi2*dt8);
    g1 = (2.0*kinE-nFree*Temp)/q1;
    vxi1 += g1*dt4;
    vxi1 *=  exp(-vxi2*dt8);
    g2 = (q1*vxi1*vxi1-Temp)/q2;
    vxi2 += g2*dt4;

    matter->setPositions(pos);
    matter->setVelocities(vel);
    return;
}


void Dynamics::langevinVerlet(double temperature)
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
    md_fcalls ++;
    noise = acc;
    mass = matter->getMasses();

    friction = - gamma * vel;
    for (long int i = 0; i < nAtoms; i++){
        if(!matter->getFixed(i)){
            for (int j = 0; j < 3; j++){
                noise(i,j) =  sqrt(4.0*gamma*kb*temperature/dt/mass[i])*gaussRandom(0.0,1.0);
            }
        }
    }
    acc = acc + friction + noise;

    vel += acc * 0.5 * dt; //first calculated velocites v(n+1/2)
    pos += vel * dt;
    matter->setPositions(pos); // update positions x(n+1)

    acc =  matter->getAccelerations();
    md_fcalls ++;
    friction = - gamma * vel;
    for (long int i = 0; i < nAtoms; i++){
        if(!matter->getFixed(i)){
            for (int j = 0; j < 3; j++){
                noise(i,j) = sqrt(4.0*gamma*kb*temperature/dt/mass[i])*gaussRandom(0.0,1.0);
            }
        }
    }

    acc = acc + friction + noise;
    vel += 0.5*dt*acc;
    matter->setVelocities(vel); // second update velocities v(n+1)
}

long Dynamics::getMDfcalls(){
    return md_fcalls;
}

long Dynamics::getMinfcalls(){
    return min_fcalls;
}

long Dynamics::getRefinefcalls(){
    return rf_fcalls;
}

bool Dynamics::checkState(Matter *matter, Matter *min1)
{
    Matter tmp(parameters);
    tmp = *matter;
    tmp.relax(true);
    min_fcalls += tmp.getForceCalls();

    if (tmp == *min1) {
        return false;
    }
    return true;
}

long Dynamics::refine(Matter *buff[],long length,Matter *min1)
{
    long a1, b1, test, refined , initial, final, diff, RefineAccuracy;
    long tmp_fcalls;
    bool ytest;

    RefineAccuracy = parameters->paraRepRefineAccuracy; 
    log("[Dynamics] Starting search for transition step with accuracy of %ld steps\n", RefineAccuracy);
    ytest = false;

    initial = 0;
    final = length - 1;
    a1 = initial;
    b1 = final;
    diff = final - initial;
    test = int((b1-a1)/2);

    tmp_fcalls = min_fcalls;
    min_fcalls = 0;
    while(diff > RefineAccuracy)
    {

     //   printf("a1 = %ld; b1= %ld; test= %ld; ytest= %d\n",a1,b1,test,ytest);
        test = a1 + int((b1-a1)/2);
        ytest = checkState(buff[test],min1);

        if (ytest == 0){
            a1 = test;
            b1 = b1;
        }
        else if (ytest == 1){
            a1 = a1;
            b1 = test;
        }
        else {
            log("Refine Step Failed ! \n");
            exit(1);
        }
        diff = abs(b1 - a1);
    }

    rf_fcalls = min_fcalls;
    min_fcalls = tmp_fcalls;

    refined = int((a1+b1)/2)+1;
    //printf("Refined mdsteps = %ld\n",nsteps_refined);
    return refined;
}

