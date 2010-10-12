/*
 *===============================================
 * EON SaddlePoint.cpp
 *===============================================
*/

#include "SaddlePoint.h"
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
    eigenMode.setZero();
    initialDisplacement.setZero();
    forceCallsSaddlePointConcave = 0;
    forceCallsSaddlePointConvex = 0;
    return;
}

SaddlePoint::~SaddlePoint(){
    clean();
    return;
}

SaddlePoint::SaddlePoint(Matter * initial_passed, Matter *saddlepassed, Parameters *parameterspassed){
    lowestEigenmode = 0;
    eigenMode.resize(saddle->numberOfAtoms(), 3);
    initialize(initial_passed, saddlepassed, parameterspassed);
    return;
}

void SaddlePoint::clean(){
    if(lowestEigenmode != 0)
    {
        delete lowestEigenmode;
        lowestEigenmode = 0;
    }
    return;
}

void SaddlePoint::initialize(Matter * initial_passed, Matter *saddlepassed, Parameters *parameterspassed)
{
    clean();
    initial=initial_passed;
    saddle = saddlepassed;
    parameters = parameterspassed;
    if(parameters->saddleLowestEigenmodeDetermination == minmodeDimer)
    {
        lowestEigenmode=new Dimer(saddle, parameters);
    }
    else if(parameters->saddleLowestEigenmodeDetermination == minmodeLanczos)
    {
        #ifdef LANCZOS_FOR_EON_HPP
            lowestEigenmode = new Lanczos(saddle, parameters);
        #else
            std::cerr << "Lanczos not available. Compile client application with option LANCZOS\n";
            exit(EXIT_FAILURE);
        #endif
    }
    nFreeCoord = 3 * saddle->numberOfFreeAtoms();
    status = statusInit;

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
    for (int i=0; i < nall; i++) 
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

long SaddlePoint::locate(Matter *min1, Matter *min2) {
    double initialEnergy;
    eigenValue = 0;
    initialEnergy = saddle->getPotentialEnergy();

    fprintf(stdout, "  Saddle point search started.\n");

    // either an initial displacement is performed and the search is started
    // or a series of jumps is performed to reach a convex region 
    if (parameters->saddleRefine) {
        lowestEigenmode->startNewSearchAndCompute(saddle, mode);
        eigenMode = lowestEigenmode->getEigenvector();
        eigenValue = lowestEigenmode->getEigenvalue();
    }else{
        if(parameters->saddleMaxJumpAttempts <= 0){
            displaceInConcaveRegion();
        }else{
            jumpToConvexRegion();
        }
    }
    fprintf(stdout, "  Saddle point displaced.\n");

    if(status == statusInit)
       searchForSaddlePoint(initialEnergy);
        
    if(status == statusInit){
        fprintf(stdout, "    Saddle point determined.\n");        
        relaxFromSaddle(min1, min2);
        fprintf(stdout, "    Minima determined.\n");
    }
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

void SaddlePoint::displaceState(Matter *matter)
{
    long nAtoms = saddle->numberOfAtoms();
    long j, indexEpiCenter = 0;
    double diffR;

    
    Matrix<double, Eigen::Dynamic, 3> initialDisplacement(nAtoms, 3);
    initialDisplacement.setZero(); 
    
    if(parameters->saddleTypePerturbation == dispNotFccOrHcp)
    {
        indexEpiCenter = EpiCenters::cnaEpiCenter(matter);            
    }
    else if(parameters->saddleTypePerturbation == dispLastAtom)
    {
        indexEpiCenter = EpiCenters::lastAtom(matter);
    }
    else if(parameters->saddleTypePerturbation == dispMinCoordinated)
    {
        indexEpiCenter = EpiCenters::minimalCoordinatedEpiCenter(matter);
    }
    else
    {
        indexEpiCenter = EpiCenters::randomFreeAtomEpiCenter(matter);
    }
    printf("Chose atom %li as the epicenter.\n", indexEpiCenter);
 
    // To keep track of free coordinates when setting atoms that have moved
    // Create an array containing initialDisplacement_ of atoms 
    // in the vicinity of the epicenter atom
    j = 0;
    for(int i = 0; i < nAtoms; i++)
    {
        if(matter->getFixed(i) == false)
        {
            diffR = matter->distance(i, indexEpiCenter);
            if(diffR < parameters->saddleWithinRadiusPerturbated)
            {
                initialDisplacement(i,0) = 2 * randomDouble() - 1;
                initialDisplacement(i,1) = 2 * randomDouble() - 1;
                initialDisplacement(i,2) = 2 * randomDouble() - 1;
            }
        }
        j++;
    }
    initialDisplacement.normalize();

    initialDisplacement *= parameters->saddleNormPerturbation;
 
    //XXX: There is probably a more idomatic way to do this with Eigen
    for(int i = 0; i < 3 * nAtoms; i++)
    {
        if(parameters->saddleMaxSinglePerturbation < initialDisplacement[i])
        {
            initialDisplacement[i] = parameters->saddleMaxSinglePerturbation;
        }
        else if(initialDisplacement[i] < -parameters->saddleMaxSinglePerturbation)
        {
            initialDisplacement[i] = -parameters->saddleMaxSinglePerturbation;
        }
    }
    // Adding the initialDisplacement
    matter->setPositions(matter->getPositions() + initialDisplacement);
 
    return;
}

Matrix<double,Eigen::Dynamic, 3> SaddlePoint::correctingForces(Matrix<double, Eigen::Dynamic, 3> force){
 
    Matrix<double, Eigen::Dynamic, 3> proj;
    proj = force.dot(eigenMode) * eigenMode.normalized();
 
    if (0 < eigenValue){
        if (parameters->saddlePerpendicularForceRatio > 0.0) {
            // reverse force parallel to eigenvector, and reduce perpendicular force
            double const d=parameters->saddlePerpendicularForceRatio;
            force = d*force - (1+d)*proj;
        }
        else {
            // Follow eigenmode
            force = -proj;
        }
    }
    else{
        // Reversing force parallel to eigenmode
        force += -2*proj;
    }
    return force;
}

void SaddlePoint::relaxFromSaddle(Matter *min1, Matter *min2){
 
    Matrix<double, Eigen::Dynamic, 3> posSaddle = saddle->getPositions();

    Matrix<double, Eigen::Dynamic, 3> displacedPos;
    //----- Initialize end -----
    //std::cout<<"relaxFromSaddle\n";
 
    // Displace saddle point configuration along the lowest eigenmode and minimize
    *min1 = *saddle;
    //XXX: the distance displced from the saddle should be a parameter
    displacedPos = posSaddle - eigenMode * 0.2;
    min1->setPositions(displacedPos);
    ConjugateGradients cgMin1(min1, parameters);
    cgMin1.fullRelax();
 
    *min2 = *saddle;
    displacedPos = posSaddle + eigenMode * 0.2;
    min2->setPositions(displacedPos);
    ConjugateGradients cgMin2(min2, parameters);  
    cgMin2.fullRelax();
 
    return;
}

void SaddlePoint::addForceCallsSaddlePoint(long fcalls, double eigenvalue){
    if(0 < eigenvalue)
        forceCallsSaddlePointConcave += fcalls;
    else
        forceCallsSaddlePointConvex += fcalls;
    return;
}

void SaddlePoint::jumpToConvexRegion(){
    long forceCallsSaddle;
    long iterations = 0;
 
    Matrix<double, Eigen::Dynamic, 3> pos = saddle->getPositions();
 
    forceCallsSaddle = saddle->getForceCalls();

    if(parameters->saddleTypePerturbation!=dispNone){
        do{
            saddle->setPositions(pos);
            displaceState(saddle);
            lowestEigenmode->startNewSearchAndCompute(saddle, initialDisplacement);
            eigenMode = lowestEigenmode->getEigenvector();
            eigenValue = lowestEigenmode->getEigenvalue();
            iterations++;
        }while((0 < eigenValue) && 
               (iterations < parameters->saddleMaxJumpAttempts));
    }
    if(0 < eigenValue)
        status = statusBadNoConvex;

    forceCallsSaddle = saddle->getForceCalls()-forceCallsSaddle;
    addForceCallsSaddlePoint(forceCallsSaddle, eigenValue);

    return;
}

void SaddlePoint::displaceInConcaveRegion()
{
    displaceState(saddle);
    lowestEigenmode->startNewSearchAndCompute(saddle, initialDisplacement);
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
 
    double *forcesStep;
    double *posStep;
    double *forces;
    double *pos;

    forcesStep = new double [nFreeCoord];
    posStep = new double [nFreeCoord];
    forces = new double [nFreeCoord];
    pos = new double[nFreeCoord];

    saddle->getFreePositions(pos);
    //----- Initialize end -----
    //std::cout<<"searchForSaddlePoint\n";
    saddle->getFreeForces(forces);
 
    eigenValue = lowestEigenmode->returnLowestEigenmode(eigenMode);
    forces = correctingForces(forces);
    ConjugateGradients cgSaddle(saddle, parameters, forces);
#ifndef NDEBUG
    saddle->matter2xyz("climb", false);
#endif
    do
    {
        forceCallsSaddle = saddle->getForceCalls();        
        cgSaddle.makeInfinitesimalStepModifiedForces(posStep, pos);
        // Determining the optimal step
        saddle->setFreePositions(posStep);
        saddle->getFreeForces(forcesStep);
        forcesStep = correctingForces(forcesStep);
/*        if(0 < eigenValue)
            maxStep = parameters->getMaxStepSizeConcave_SP();
        else
            maxStep = parameters->getMaxStepSizeConvex_SP();*/
        maxStep = parameters->saddleMaxStepSize;
        cgSaddle.getNewPosModifiedForces(pos, forces, forcesStep, maxStep);
        // The system (saddle) is moved to a new configuration
        saddle->setFreePositions(pos);
        saddle->getFreeForces(forces);
        // The new lowest eigenvalue
        lowestEigenmode_->moveAndCompute(saddle);
        eigenValue_ = lowestEigenmode_->returnLowestEigenmode(eigenMode_);
        // Updating the conjugated object to the new configuration
        forces = correctingForces(forces);
        cgSaddle.setFreeAtomForcesModifiedForces(forces);
        if(eigenValue_ < 0)
        {
            converged = cgSaddle.isItConverged(parameters->saddleConverged);
            concaveSeries = 0;
        }
        else
        {
            converged = false;
            concaveSeries = concaveSeries + 1;
        }
        forceCallsSaddle = saddle->getForceCalls()-forceCallsSaddle;        
        addForceCallsSaddlePoint(forceCallsSaddle, eigenValue_);

        iterations++;
        #ifndef NDEBUG
            printf("climb = %ld, max force = %f\n", iterations, saddle->maxForce());
            saddle->matter2xyz("climb", true);
        #endif
        energySaddle = saddle->potentialEnergy();
    }while(!converged && 
           (iterations < parameters->saddleMaxIterations) && 
           (energySaddle-initialEnergy < parameters->saddleMaxEnergy));
    if(!converged){
        if(parameters->saddleMaxIterations <= iterations) {
            status = statusBadMaxIterations;
        }
        if(parameters->saddleMaxEnergy <= energySaddle-initialEnergy) {
            status = statusBadHighBarrier;
        }
    }
    delete [] forcesStep;
    delete [] posStep;
    delete [] forces;
    delete [] pos;
    return; 
}

LowestEigenmodeInterface const * SaddlePoint::getLowestEigenmode() const
{
      return lowestEigenmode;
}
