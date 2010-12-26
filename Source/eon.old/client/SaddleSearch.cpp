//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "SaddleSearch.h"
#ifdef WITH_LANCZOS
    #include "Lanczos/lanczos_for_eon.hpp"
#endif
#include "ConjugateGradients.h"
#include "HelperFunctions.h"
#include "Dimer.h"
#include "EpiCenters.h"
#include "Constants.h"

#include <cstdlib>

using namespace helper_functions;

SaddlePoint::SaddlePoint(){
    lowestEigenmode = 0;
    forceCallsSaddlePointConcave = 0;
    forceCallsSaddlePointConvex = 0;
//    forceCallsMinimization = 0;
    eigenValue = 0;
    return;
}

SaddlePoint::~SaddlePoint(){
    clean();
    return;
}

SaddlePoint::SaddlePoint(Matter *initialPassed, Matter *saddlePassed, Parameters *parametersPassed){
    lowestEigenmode = 0;
    eigenValue = 0;
    initialize(initialPassed, saddlePassed, parametersPassed);
    return;
}

void SaddlePoint::clean(){
    if(lowestEigenmode != 0)
    {
        delete lowestEigenmode;
        lowestEigenmode = 0;
    }
    eigenValue = 0;
    return;
}

void SaddlePoint::initialize(Matter *initialPassed, Matter *saddlePassed, Parameters *parametersPassed)
{
    clean();
    initial = initialPassed;
    saddle = saddlePassed;
    parameters = parametersPassed;
    eigenMode.resize(saddlePassed->numberOfAtoms(), 3);
    eigenMode.setZero();
    if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
    {
        lowestEigenmode=new Dimer(saddle, parameters);
    }
    else if(parameters->saddleMinmodeMethod == MINMODE_LANCZOS)
    {
        #ifdef LANCZOS_FOR_EON_HPP
            lowestEigenmode = new Lanczos(saddle, parameters);
        #else
            std::cerr << "Lanczos not available. Compile client application with option LANCZOS\n";
            exit(EXIT_FAILURE);
        #endif
    }
    status = STATUS_INIT;
    eigenValue = 0;

    return;
}

void SaddlePoint::loadMode(string filename) {
    FILE *modeFile;

    modeFile = fopen(filename.c_str(), constants::READ.c_str());
    loadMode(modeFile);
    fclose(modeFile);
}

void SaddlePoint::loadMode(FILE *modeFile){
    long nall=0, nfree=0;
    fscanf(modeFile, "%ld %ld", &nall, &nfree);
    nFreeCoord = nfree;
    mode.resize(nall/3, 3);
    mode.setZero();
    for (int i=0; i < nall/3; i++) 
    {
        for(int j=0; j<3; j++)
        {
            fscanf(modeFile, "%lf", &mode(i,j));
        }
    }
}

void SaddlePoint::saveMode(FILE *modeFile)
{
    long const nAtoms = saddle->numberOfAtoms();
    fprintf(modeFile, "%ld %ld\n", nAtoms*3, nFreeCoord);
    for (long i=0; i < nAtoms; ++i) {
        if (saddle->getFixed(i)) {
            fprintf(modeFile, "0 0 0\n");
        }
        else {
            fprintf(modeFile, "%lf\t%lf \t%lf\n", eigenMode(i,0), eigenMode(i,1), eigenMode(i,2));
        }
    }
    return;
}

long SaddlePoint::locate(void){//(Matter *min1, Matter *min2) {
    double initialEnergy;
    eigenValue = 0;
    initialEnergy = saddle->getPotentialEnergy();

    fprintf(stdout, "  Saddle point search started.\n");

    // either an initial displacement is performed and the search is started
    // or a series of jumps is performed to reach a convex region 
    if (parameters->saddleDisplace) 
    {
        displaceAndSetMode(saddle);
    }
    
    lowestEigenmode->initialize(saddle, mode);
    lowestEigenmode->compute(saddle);
    eigenMode = lowestEigenmode->getEigenvector();
    eigenValue = lowestEigenmode->getEigenvalue();
    
    fprintf(stdout, "  Saddle point displaced.\n");

    if(status == STATUS_INIT)
    {
       searchForSaddlePoint(initialEnergy);
    }

//    if(status == STATUS_INIT)
//    {
//        fprintf(stdout, "    Saddle point determined.\n");
//        relaxFromSaddle(min1, min2);
//        fprintf(stdout, "    Minima determined.\n");
//    }
    return(status);
}

long SaddlePoint::getnFreeCoord() const
{
    return nFreeCoord;
}

Matrix<double, Eigen::Dynamic, 3> SaddlePoint::getEigenMode() 
{
    return eigenMode;
}

void SaddlePoint::displaceAndSetMode(Matter *matter)
{
    long nAtoms = saddle->numberOfAtoms();
    long j, indexEpiCenter = 0;
    double diffR;

    Matrix<double, Eigen::Dynamic, 3> initialDisplace(nAtoms, 3);
    initialDisplace.setZero(); 

    if(parameters->saddleDisplaceType == DISP_NOT_FCC_OR_HCP)
    {
        indexEpiCenter = EpiCenters::cnaEpiCenter(matter, parameters->neighborCutoff);
    }
    else if(parameters->saddleDisplaceType == DISP_LAST_ATOM)
    {
        indexEpiCenter = EpiCenters::lastAtom(matter);
    }
    else if(parameters->saddleDisplaceType == DISP_MIN_COORDINATED)
    {
        indexEpiCenter = EpiCenters::minCoordinatedEpiCenter(matter, parameters->neighborCutoff);
    }
    else
    {
        indexEpiCenter = EpiCenters::randomFreeAtomEpiCenter(matter);
    }
    printf("Chose atom %li as the epicenter.\n", indexEpiCenter);

    // To keep track of free coordinates when setting atoms that have moved
    // Create an array containing initialDisplace_ of atoms 
    // in the vicinity of the epicenter atom
    j = 0;
    for(int i = 0; i < nAtoms; i++)
    {
        if(matter->getFixed(i) == false)
        {
            diffR = matter->distance(i, indexEpiCenter);
            if(diffR < parameters->saddleDisplaceRadius)
            {
                initialDisplace(i,0) = 2 * randomDouble() - 1;
                initialDisplace(i,1) = 2 * randomDouble() - 1;
                initialDisplace(i,2) = 2 * randomDouble() - 1;
            }
        }
        j++;
    }
    initialDisplace.normalize();
    initialDisplace *= parameters->saddleDisplaceMagnitude;

    //XXX: There is probably a more idomatic way to do this with Eigen
    for(int i = 0; i < 3 * nAtoms; i++)
    {
        if(parameters->saddleMaxSingleDisplace < initialDisplace[i])
        {
            initialDisplace[i] = parameters->saddleMaxSingleDisplace;
        }
        else if(initialDisplace[i] < -parameters->saddleMaxSingleDisplace)
        {
            initialDisplace[i] = -parameters->saddleMaxSingleDisplace;
        }
    }
    // Adding the initialDisplace
    matter->setPositions(matter->getPositions() + initialDisplace);
    
    // Sets the initial mode for the SP search
    mode = initialDisplace;

    return;
}

Matrix<double,Eigen::Dynamic, 3> SaddlePoint::correctingForces(Matrix<double, Eigen::Dynamic, 3> force){
 
    Matrix<double, Eigen::Dynamic, 3> proj;
    proj = (force.cwise() * eigenMode).sum() * eigenMode.normalized();

    if (0 < eigenValue){
        if (parameters->saddlePerpForceRatio > 0.0) {
            // reverse force parallel to eigenvector, and reduce perpendicular force
            double const d=parameters->saddlePerpForceRatio;
            force = d*force - (1+d)*proj;
        }
        else {
            // follow eigenmode
            force = -proj;
        }
    }
    else{
        // reversing force parallel to eigenmode
        force += -2*proj;
    }
    return force;
}
/*
void SaddlePoint::relaxFromSaddle(Matter *min1, Matter *min2){

    Matrix<double, Eigen::Dynamic, 3> posSaddle = saddle->getPositions();

    Matrix<double, Eigen::Dynamic, 3> displacedPos;
    //----- Initialize end -----
    //std::cout<<"relaxFromSaddle\n";

    // Displace saddle point configuration along the lowest eigenmode and minimize
    *min1 = *saddle;
    //XXX: the distance displaced from the saddle should be a parameter
    displacedPos = posSaddle - eigenMode * 0.2;
    min1->setPositions(displacedPos);
    ConjugateGradients cgMin1(min1, parameters);
    cgMin1.fullRelax();
    forceCallsMinimization += cgMin1.totalForceCalls;

    *min2 = *saddle;
    displacedPos = posSaddle + eigenMode * 0.2;
    min2->setPositions(displacedPos);
    ConjugateGradients cgMin2(min2, parameters);  
    cgMin2.fullRelax();
    forceCallsMinimization += cgMin2.totalForceCalls;

    return;
}
*/
void SaddlePoint::addForceCallsSaddlePoint(long fcalls, double eigenvalue){
    if(0 < eigenvalue)
        forceCallsSaddlePointConcave += fcalls;
    else
        forceCallsSaddlePointConvex += fcalls;
    return;
}

void SaddlePoint::displaceInConcaveRegion()
{
    displaceAndSetMode(saddle);
    
    lowestEigenmode->initialize(saddle, mode);
    lowestEigenmode->compute(saddle);
    eigenMode = lowestEigenmode->getEigenvector();
    eigenValue = lowestEigenmode->getEigenvalue();

    return;
}

void SaddlePoint::searchForSaddlePoint(double initialEnergy)
{
    bool converged = false;
    long forceCallsSaddle;
    long iterations = 0;
    long concaveSeries = 0;
    double maxStep;
    double energySaddle;

    Matrix<double, Eigen::Dynamic, 3> forcesStep;
    Matrix<double, Eigen::Dynamic, 3> posStep;
    Matrix<double, Eigen::Dynamic, 3> forces;
    Matrix<double, Eigen::Dynamic, 3> pos;

    pos = saddle->getPositions();
    //----- Initialize end -----
    //std::cout<<"searchForSaddlePoint\n";
    forces = saddle->getForces();

    eigenValue = lowestEigenmode->getEigenvalue();
    eigenMode = lowestEigenmode->getEigenvector();
    forces = correctingForces(forces);
    ConjugateGradients cgSaddle(saddle, parameters, forces);
    #ifndef NDEBUG
        static int run;
        ostringstream climb;
        climb << "climb_" << run;
        initial->matter2xyz(climb.str(), false);
        saddle->matter2xyz(climb.str(), true);
        ++run;
        if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
        {
            printf("DIMER ---------------------------------------------------------------------------------------------\n");    
            printf("DIMER  %9s   %9s   %9s   %9s   %9s   %9s  %9s   %9s\n", "Step", "Force", "Torque", 
                   "Energy", "Curvature", "Angle", "Rotations", "Step Size");
            printf("DIMER ---------------------------------------------------------------------------------------------\n");    
        }
        else {
            printf("LANCZOS ---------------------------------------------------------------------------------------\n");    
            printf("LANCZOS  %9s  %9s  %9s  %9s\n", "Step", "Force", "Energy", "Step Size");
            printf("LANCZOS ---------------------------------------------------------------------------------------\n");    
        };
    #endif
    do
    {
        forceCallsSaddle = saddle->getForceCalls();        
        posStep = cgSaddle.makeInfinitesimalStepModifiedForces(pos);
        // Determining the optimal step
        saddle->setPositions(posStep);
        forcesStep = saddle->getForces();
        forcesStep = correctingForces(forcesStep);
        
        maxStep = parameters->saddleMaxStepSize;
        pos = cgSaddle.getNewPosModifiedForces(pos, forces, forcesStep, maxStep);
        #ifndef NDEBUG
        double stepSize = sqrt(((saddle->getPositions() - pos).cwise().square()).sum());
        #endif
        // The system (saddle) is moved to a new configuration
        saddle->setPositions(pos);
        forces = saddle->getForces();
        // The new lowest eigenvalue
        lowestEigenmode->compute(saddle);
        eigenValue = lowestEigenmode->getEigenvalue();
        eigenMode = lowestEigenmode->getEigenvector();
        // Updating the conjugated object to the new configuration
        forces = correctingForces(forces);
        cgSaddle.setForces(forces);
        if(eigenValue < 0)
        {
            converged = cgSaddle.isItConverged(parameters->optConvergedForce);
            concaveSeries = 0;
        }
        else
        {
            converged = false;
            concaveSeries = concaveSeries + 1;
        }
        forceCallsSaddle = saddle->getForceCalls()-forceCallsSaddle;        
        addForceCallsSaddlePoint(forceCallsSaddle, eigenValue);

        iterations++;
        #ifndef NDEBUG
            if(parameters->saddleMinmodeMethod == MINMODE_DIMER)        
            {
                double *stats = lowestEigenmode->stats;
                printf("DIMER  %9ld  % 9.3e  % 9.3e  % 10.3f  % 9.3e  % 9.3e  %9d  % 9.3e \n", iterations, 
                       sqrt((saddle->getForces().cwise().square()).sum()), stats[0], 
                       saddle->getPotentialEnergy(), stats[1], stats[2], (int)stats[3], stepSize);
            } else {
                printf("LANCZOS  %9ld  % 9.5f  % 9.5f  % 9.5f\n", iterations, 
                       sqrt((saddle->getForces().cwise().square()).sum()), saddle->getPotentialEnergy(), stepSize);
            }
            saddle->matter2xyz(climb.str(), true);
        #endif
        energySaddle = saddle->getPotentialEnergy();
    }while(!converged && 
           (iterations < parameters->saddleMaxIterations) && 
           (energySaddle-initialEnergy < parameters->saddleMaxEnergy));
    if(!converged)
    {
        if(parameters->saddleMaxIterations <= iterations) 
        {
            status = STATUS_BAD_MAX_ITERATIONS;
        }
        if(parameters->saddleMaxEnergy <= energySaddle-initialEnergy) 
        {
            status = STATUS_BAD_HIGH_BARRIER;
        }
    }
    return; 
}

LowestEigenmodeInterface const * SaddlePoint::getLowestEigenmode() const
{
      return lowestEigenmode;
}

Matrix<double, Eigen::Dynamic, 3> SaddlePoint::getSaddlePositions()
{
      return saddle->getPositions();	
}

