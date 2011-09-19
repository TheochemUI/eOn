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
#include "Log.h"

#include <cstdlib>

using namespace helper_functions;

const char SaddleSearch::MINMODE_DIMER[] =           "dimer";
const char SaddleSearch::MINMODE_LANCZOS[] =         "lanczos";
const char SaddleSearch::MINMODE_EXACT[] =           "exact";
const char SaddleSearch::DISP_LOAD[] =               "load";
const char SaddleSearch::DISP_NOT_FCC_OR_HCP[] =     "not_fcc_hcp_coordinated";
const char SaddleSearch::DISP_MIN_COORDINATED[] =    "least_coordinated";
const char SaddleSearch::DISP_LAST_ATOM[] =          "last_atom";
const char SaddleSearch::DISP_RANDOM[] =             "random";

SaddleSearch::SaddleSearch()
{
    lowestEigenmode = 0;
    forceCallsSaddleSearchConcave = 0;
    forceCallsSaddleSearchConvex = 0;
    eigenValue = 0;
    return;
}

SaddleSearch::~SaddleSearch()
{
    clean();
    return;
}

SaddleSearch::SaddleSearch(Matter *initialPassed, Matter *saddlePassed, Parameters *parametersPassed)
{
    lowestEigenmode = 0;
    eigenValue = 0;
    initialize(initialPassed, saddlePassed, parametersPassed);
    return;
}

void SaddleSearch::clean()
{
    if(lowestEigenmode != 0)
    {
        delete lowestEigenmode;
        lowestEigenmode = 0;
    }
    eigenValue = 0;
    return;
}

void SaddleSearch::initialize(Matter *initialPassed, Matter *saddlePassed, Parameters *parametersPassed)
{
    clean();
    initial = initialPassed;
    saddle = saddlePassed;
    parameters = parametersPassed;
    if (parameters->saddleDisplaceType == DISP_LOAD)
    {
        initialDisplacement = saddlePassed->getPositions() - initialPassed->getPositions();
        initialDisplacement /= initialDisplacement.norm();
    }
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
    if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
    {
        log("[Dimer]  %9s   %9s   %16s   %9s   %9s   %9s   %9s   %9s\n", 
            "Step", "Step Size", "Energy", "Force", "Curvature", 
            "Torque", "Angle", "Rotations");
    }
    else 
    {
        log("[Lanczos]  %9s  %9s  %16s  %9s  %9s\n", 
            "Step", "Step Size", "Energy", "Force", "Curvature");
    }
    return;
}

void SaddleSearch::loadMode(string filename)
{
    FILE *modeFile;
    modeFile = fopen(filename.c_str(), constants::READ.c_str());
    if (!modeFile) {
        cerr << "File " << filename << " was not found.\n";
        log("Stop\n");
        exit(1);
    }
    loadMode(modeFile);
    fclose(modeFile);
}

void SaddleSearch::loadMode(FILE *modeFile)
{
    long const nAtoms = saddle->numberOfAtoms();
    mode.resize(nAtoms, 3);
    mode.setZero();
    for (int i=0; i < nAtoms; i++)
    {
        fscanf(modeFile, "%lf %lf %lf", &mode(i,0), &mode(i,1), &mode(i,2));
    }
}

void SaddleSearch::saveMode(FILE *modeFile)
{
    long const nAtoms = saddle->numberOfAtoms();
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

long SaddleSearch::locate(void)
{
    double initialEnergy;
    eigenValue = 0;
    initialEnergy = saddle->getPotentialEnergy();

    log("  Saddle point search started from reactant with energy %f eV.\n", initialEnergy);

    // either an initial displacement is performed and the search is started
    // or a series of jumps is performed to reach a convex region

    // the displacement is done on the client
    if (parameters->saddleDisplaceType != DISP_LOAD) 
    {
        displaceAndSetMode(saddle);
    }

    log("  Saddle point displaced.\n");

    if(status == STATUS_INIT)
    {
       searchForSaddlePoint(initialEnergy);
    }

    return(status);
}

AtomMatrix SaddleSearch::getEigenMode() 
{
    return eigenMode;
}

double SaddleSearch::getEigenValue()
{
    return eigenValue;
}

void SaddleSearch::displaceAndSetMode(Matter *matter)
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
    log("Chose atom %li as the epicenter.\n", indexEpiCenter);

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

AtomMatrix SaddleSearch::projectedForce(AtomMatrix force){

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


void SaddleSearch::addForceCallsSaddleSearch(long fcalls, double eigenvalue){
    if(0 < eigenvalue)
        forceCallsSaddleSearchConcave += fcalls;
    else
        forceCallsSaddleSearchConvex += fcalls;
    return;
}

void SaddleSearch::displaceInConcaveRegion()
{
    displaceAndSetMode(saddle);
    
    lowestEigenmode->compute(saddle,mode);
    eigenMode = lowestEigenmode->getEigenvector();
    eigenValue = lowestEigenmode->getEigenvalue();

    return;
}

void SaddleSearch::searchForSaddlePoint(double initialEnergy)
{
    ostringstream climb;
    climb << "climb";
    if(parameters->writeMovies)
    {
        if (parameters->checkpoint) {
            if (fopen("climb.con", "r") == NULL) {
                initial->matter2con(climb.str(), false);
            }
        }else{
            initial->matter2con(climb.str(), false);
        }
        saddle->matter2con(climb.str(), true);
    }

    bool converged = false;
    long forceCallsSaddle;
    long concaveSeries = 0;
    double maxStep;
    double energySaddle;
    iterations = 0;

    AtomMatrix forcesStep;
    AtomMatrix posStep;
    AtomMatrix forces;
    AtomMatrix pos;

    pos = saddle->getPositions();

    // GH: this should be generalized to other optimizers
    ConjugateGradients cgSaddle(saddle, parameters);

    eigenMode = mode;
    forces = saddle->getForces();

    do
    {
        forceCallsSaddle = saddle->getForceCalls();

        // The new lowest eigenvalue
        lowestEigenmode->compute(saddle, eigenMode);
        eigenValue = lowestEigenmode->getEigenvalue();
        eigenMode = lowestEigenmode->getEigenvector();

        // Updating the conjugated object to the new configuration
        forces = projectedForce(forces);
        cgSaddle.setForces(forces);

        // Determine a CG step.
        posStep = cgSaddle.makeInfinitesimalStepModifiedForces(pos);
        saddle->setPositions(posStep);
        forcesStep = saddle->getForces();
        forcesStep = projectedForce(forcesStep);
        maxStep = parameters->saddleMaxStepSize;
        pos = cgSaddle.getNewPosModifiedForces(pos, forces, forcesStep, maxStep);

        // The system (saddle) is moved to a new configuration
        double stepSize = (saddle->pbc(saddle->getPositions() - pos )).norm();
        saddle->setPositions(pos);
        forces = saddle->getForces();

        if(eigenValue < 0)
        {
            converged = cgSaddle.isItConverged(parameters->saddleConvergedForce);
            concaveSeries = 0;
        }
        else
        {
            converged = false;
            concaveSeries = concaveSeries + 1;
        }
        forceCallsSaddle = saddle->getForceCalls()-forceCallsSaddle;
        addForceCallsSaddleSearch(forceCallsSaddle, eigenValue);
        iterations++;
        if(parameters->saddleMinmodeMethod == MINMODE_DIMER)
        {
            log("[Dimer]  %9ld  % 9.3e   %16.4f  % 9.3e  % 9.3e  % 9.3e  % 9.3e   % 9d\n",
                        iterations, stepSize, saddle->getPotentialEnergy(),
                        saddle->getForces().norm(),
                        lowestEigenmode->getEigenvalue(),
                        lowestEigenmode->statsTorque,
                        lowestEigenmode->statsAngle,
                        (int)lowestEigenmode->statsRotations);
        }
        else 
        {
            log("[Lanczos]  %9ld  % 9.3f   %16.4f  % 9.3f  % 9.3f\n", 
                iterations, stepSize, saddle->getPotentialEnergy(),
                saddle->getForces().norm(),
                lowestEigenmode->getEigenvalue());
        }
        if(parameters->writeMovies)
        {
            saddle->matter2con(climb.str(), true);
        }
        energySaddle = saddle->getPotentialEnergy();
        if (parameters->checkpoint) {
            FILE *fileMode;
            fileMode = fopen("mode_checkpoint.dat", "wb");
            saveMode(fileMode);
            fclose(fileMode);

            FILE *fileSaddle;
            fileSaddle = fopen("displacement_checkpoint.con", "wb");
            saddle->matter2con(fileSaddle);
            fclose(fileSaddle);
        }
    }while(!converged && 
           (iterations < parameters->saddleMaxIterations) && 
           (energySaddle-initialEnergy < parameters->saddleMaxEnergy));
    if(!converged) {
        if (parameters->saddleMaxIterations <= iterations) {
            status = STATUS_BAD_MAX_ITERATIONS;
        }
        if (parameters->saddleMaxEnergy <= energySaddle-initialEnergy) {
            status = STATUS_BAD_HIGH_BARRIER;
        }
    }
    return;
}

LowestEigenmodeInterface const * SaddleSearch::getLowestEigenmode() const
{
      return lowestEigenmode;
}

AtomMatrix SaddleSearch::getSaddlePositions()
{
      return saddle->getPositions();
}

