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
#include "ConjugateGradients.h"
#include "HelperFunctions.h"
#include "Lanczos.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "ExactMinMode.h"
#include "EpiCenters.h"
#include "Constants.h"

#include <cstdlib>

using namespace helper_functions;

const char SaddlePoint::MINMODE_DIMER[] =           "dimer";
const char SaddlePoint::MINMODE_LANCZOS[] =         "lanczos";
const char SaddlePoint::MINMODE_EXACT[] =           "exact";
const char SaddlePoint::DISP_LOAD[] =               "client_load";
const char SaddlePoint::DISP_NOT_FCC_OR_HCP[] =     "client_not_fcc_hcp_coordinated";
const char SaddlePoint::DISP_MIN_COORDINATED[] =    "client_least_coordinated";
const char SaddlePoint::DISP_LAST_ATOM[] =          "client_last_atom";
const char SaddlePoint::DISP_RANDOM[] =             "client_random";

SaddlePoint::SaddlePoint()
{
    lowestEigenmode = 0;
    forceCallsSaddlePointConcave = 0;
    forceCallsSaddlePointConvex = 0;
    eigenValue = 0;
    foundNegativeMode = false;
    return;
}

SaddlePoint::~SaddlePoint()
{
    clean();
    return;
}

SaddlePoint::SaddlePoint(Matter *initialPassed, Matter *saddlePassed, Parameters *parametersPassed)
{
    lowestEigenmode = 0;
    eigenValue = 0;
    initialize(initialPassed, saddlePassed, parametersPassed);
    return;
}

void SaddlePoint::clean()
{
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
    initialDisplacement = saddlePassed->getPositions() - initialPassed->getPositions();
    initialDisplacement /= initialDisplacement.norm();
    parameters = parametersPassed;
    eigenMode.resize(saddlePassed->numberOfAtoms(), 3);
    eigenMode.setZero();
    if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
    {
        if(parameters->dimerImproved)
        {
            lowestEigenmode = new ImprovedDimer(saddle, parameters);
        }
        else
        {
            lowestEigenmode = new Dimer(saddle, parameters);
        }
        
    }
    else if(parameters->saddleMinmodeMethod == MINMODE_LANCZOS)
    {
        lowestEigenmode = new Lanczos(saddle, parameters);
    }
    else if(parameters->saddleMinmodeMethod == MINMODE_EXACT)
    {
        lowestEigenmode = new ExactMinMode(saddle, parameters);
    }
    status = STATUS_INIT;
    eigenValue = 0;
    FILE *fp = fopen("saddlesearch.dat", "w");
    if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
    {
        fprintf(fp, "DIMER  %9s   %9s   %9s   %9s   %9s   %9s   %9s   %9s\n", 
                    "Step", "Step Size", "Energy", "Force", "Curvature", 
                    "Torque", "Angle", "Rotations");
    }
    else 
    {
        fprintf(fp, "LANCZOS  %9s  %9s  %9s  %9s  %9s\n", 
                    "Step", "Step Size", "Energy", "Force", "Curvature");
    }
    fclose(fp);
    return;
}

void SaddlePoint::loadMode(string filename)
{
    FILE *modeFile;
    modeFile = fopen(filename.c_str(), constants::READ.c_str());
    if (!modeFile) {
        cerr << "File " << filename << " was not found.\n";
        printf("Stop\n");
        exit(1);
    }
    loadMode(modeFile);
    fclose(modeFile);
}

void SaddlePoint::loadMode(FILE *modeFile)
{
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

long SaddlePoint::locate(void)
{
    double initialEnergy;
    eigenValue = 0;
    initialEnergy = saddle->getPotentialEnergy();

    fprintf(stdout, "  Saddle point search started.\n");

    // either an initial displacement is performed and the search is started
    // or a series of jumps is performed to reach a convex region

    // the displacement is done on the client
    if (parameters->saddleDisplaceType != DISP_LOAD) 
    {
        displaceAndSetMode(saddle);
    }

    fprintf(stdout, "  Saddle point displaced.\n");

    if(status == STATUS_INIT)
    {
       searchForSaddlePoint(initialEnergy);
    }

    return(status);
}

long SaddlePoint::getnFreeCoord() const
{
    return nFreeCoord;
}

AtomMatrix SaddlePoint::getEigenMode() 
{
    return eigenMode;
}

void SaddlePoint::displaceAndSetMode(Matter *matter)
{
    long nAtoms = saddle->numberOfAtoms();
    long j, indexEpiCenter = 0;
    double diffR;

    AtomMatrix initialDisplace(nAtoms, 3);
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
    else if(parameters->saddleDisplaceType == DISP_RANDOM)
    {
        indexEpiCenter = EpiCenters::randomFreeAtomEpiCenter(matter);
    }
    printf("Chose atom %li as the epicenter.\n", indexEpiCenter);

    initialDisplace(indexEpiCenter,0) = gaussRandom(0., parameters->saddleDisplaceMagnitude);
    initialDisplace(indexEpiCenter,1) = gaussRandom(0., parameters->saddleDisplaceMagnitude);
    initialDisplace(indexEpiCenter,2) = gaussRandom(0., parameters->saddleDisplaceMagnitude);

    // To keep track of free coordinates when setting atoms that have moved
    // Create an array containing initialDisplace_ of atoms in the vicinity of the epicenter atom
    j = 0;
    for(int i = 0; i < nAtoms; i++)
    {
        if(matter->getFixed(i) == false)
        {
            diffR = matter->distance(i, indexEpiCenter);
            if(diffR < parameters->saddleDisplaceRadius)
            {
                initialDisplace(i,0) = gaussRandom(0, parameters->saddleDisplaceMagnitude);
                initialDisplace(i,1) = gaussRandom(0, parameters->saddleDisplaceMagnitude);
                initialDisplace(i,2) = gaussRandom(0, parameters->saddleDisplaceMagnitude);
            }
        }
        j++;
    }

    //XXX: There is probably a more idiomatic way to do this with Eigen
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
    mode.normalize();

    return;
}

AtomMatrix SaddlePoint::projectedForce(AtomMatrix force){

    AtomMatrix proj;
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
    
    lowestEigenmode->compute(saddle,mode);
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

    AtomMatrix forcesStep;
    AtomMatrix posStep;
    AtomMatrix forces;
    AtomMatrix pos;

    if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
    {
        printf("DIMER  %9s   %9s   %9s   %9s   %9s   %9s   %9s   %9s\n", 
               "Step", "Step Size", "Energy", "Force", "Curvature", 
               "Torque", "Angle", "Rotations");
    }
    else 
    {
        printf("LANCZOS  %9s  %9s  %9s  %9s  %9s\n", 
               "Step", "Step Size", "Energy", "Force", "Curvature");
    }

    pos = saddle->getPositions();
    //----- Initialize end -----
    //std::cout<<"searchForSaddlePoint\n";
    forces = saddle->getForces();

    lowestEigenmode->compute(saddle, mode);
    eigenValue = lowestEigenmode->getEigenvalue();
    eigenMode = lowestEigenmode->getEigenvector();
    forces = projectedForce(forces);
    // GH: this should be generalized to other optimizers
    ConjugateGradients cgSaddle(saddle, parameters, forces);
    ostringstream climb;
    climb << "climb";
    if(parameters->writeMovies)
    {
        initial->matter2con(climb.str(), false);
        saddle->matter2con(climb.str(), true);
    }
    do
    {
        forceCallsSaddle = saddle->getForceCalls();
        if((eigenValue < -parameters->saddleNonzeroMode || foundNegativeMode) || !parameters->saddlePushDisplacement)
        {
            // Determine a CG step.
            posStep = cgSaddle.makeInfinitesimalStepModifiedForces(pos);
            saddle->setPositions(posStep);
            forcesStep = saddle->getForces();
            forcesStep = projectedForce(forcesStep);
            maxStep = parameters->saddleMaxStepSize;
            pos = cgSaddle.getNewPosModifiedForces(pos, forces, forcesStep, maxStep);
            foundNegativeMode = true;
        }
        else
        {
            pos = saddle->getPositions() + initialDisplacement * parameters->saddleConcaveStepSize;
        }
        // The system (saddle) is moved to a new configuration
        double stepSize = (saddle->pbc(saddle->getPositions() - pos )).norm();
        saddle->setPositions(pos);
        forces = saddle->getForces();
        // The new lowest eigenvalue
        lowestEigenmode->compute(saddle,eigenMode);
        eigenValue = lowestEigenmode->getEigenvalue();
        eigenMode = lowestEigenmode->getEigenvector();
        // Updating the conjugated object to the new configuration
        forces = projectedForce(forces);
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
        FILE *fp = fopen("saddlesearch.dat", "a");
        if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
        {
            fprintf(fp, "DIMER  %9ld  % 9.3e  % 9.3e  % 9.3e  % 9.3e  % 9.3e  % 9.3e   % 9d\n",
                        iterations, stepSize, saddle->getPotentialEnergy(),
                        saddle->getForces().norm(),
                        lowestEigenmode->getEigenvalue(),
                        lowestEigenmode->statsTorque,
                        lowestEigenmode->statsAngle,
                        (int)lowestEigenmode->statsRotations);
            printf("DIMER  %9ld  % 9.3e  % 9.3e  % 9.3e  % 9.3e  % 9.3e  % 9.3e   % 9d\n",
                   iterations, stepSize, saddle->getPotentialEnergy(),
                   saddle->getForces().norm(),
                   lowestEigenmode->getEigenvalue(),
                   lowestEigenmode->statsTorque,
                   lowestEigenmode->statsAngle,
                   (int)lowestEigenmode->statsRotations);
        }
        else 
        {
            fprintf(fp, "LANCZOS  %9ld  % 9.3f  % 9.3f  % 9.3f  % 9.3f\n", 
                        iterations, stepSize, saddle->getPotentialEnergy(),
                        saddle->getForces().norm(),
                        lowestEigenmode->getEigenvalue());
            printf("LANCZOS  %9ld  % 9.3f  % 9.3f  % 9.3f  % 9.3f\n", 
                   iterations, stepSize, saddle->getPotentialEnergy(),
                   saddle->getForces().norm(),
                   lowestEigenmode->getEigenvalue());
        }
        fclose(fp);
        if(parameters->writeMovies)
        {
            saddle->matter2con(climb.str(), true);
        }
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

AtomMatrix SaddlePoint::getSaddlePositions()
{
      return saddle->getPositions();
}

