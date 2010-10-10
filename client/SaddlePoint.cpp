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
    lowestEigenmode_ = 0;
    eigenMode_ = 0;
    initialDisplacement_ = 0;
    forceCallsSaddlePointConcave_ = 0;
    forceCallsSaddlePointConvex_ = 0;
    return;
}

SaddlePoint::~SaddlePoint(){
    clean();
    return;
}

SaddlePoint::SaddlePoint(Matter * initial, Matter *saddle, Parameters *parameters){
    lowestEigenmode_ = 0;
    eigenMode_ = 0;
    initialize(initial, saddle, parameters);
    return;
}

void SaddlePoint::clean(){
    // saddle_ should not be deleted, copy of a pointer that was passed in
    // delete saddle_
    if(lowestEigenmode_ != 0){
        delete lowestEigenmode_;
        lowestEigenmode_ = 0;
    }
    if(eigenMode_ != 0){
        delete eigenMode_;
        eigenMode_ = 0;
    }
    return;
}

void SaddlePoint::initialize(Matter * initial, Matter *saddle, Parameters *parameters)
{
    clean();
    initial_=initial;
    saddle_ = saddle;
    parameters_ = parameters;
    if(parameters_->saddleLowestEigenmodeDetermination == minmodeDimer)
    {
        lowestEigenmode_=new Dimer(saddle_, parameters_);
    }
    else if(parameters_->saddleLowestEigenmodeDetermination == minmodeLanczos)
    {
        #ifdef LANCZOS_FOR_EON_HPP
            lowestEigenmode_ = new Lanczos(saddle_, parameters_);
        #else
            std::cerr << "Lanczos not available. Compile client application with option LANCZOS\n";
            exit(EXIT_FAILURE);
        #endif
    }
    nFreeCoord_ = 3 * saddle->numberOfFreeAtoms();
    eigenMode_ = new double[nFreeCoord_];
    status_ = statusInit;

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
    mode = new double[nall];
    for (int i=0, j=0; i < nall; ++i) {
        fscanf(modeFile, "%lf", &mode[j]);
        ++j;
    }
}

void SaddlePoint::saveMode(FILE *modeFile)
{
    long const nAtoms = saddle_->numberOfAtoms();
    fprintf(modeFile, "%ld %ld\n", nAtoms*3, nFreeCoord_);
    for (long i=0, j=0; i < nAtoms; ++i) {
        if (saddle_->getFixed(i)) {
            fprintf(modeFile, "0 0 0\n");
        }
        else {
            fprintf(modeFile, "%lf\t%lf \t%lf\n", eigenMode_[j], eigenMode_[j+1], eigenMode_[j+2]);
            j+=3;
        }
    }
    return;
}

long SaddlePoint::locate(Matter *min1, Matter *min2) {
    double initialEnergy;
    eigenValue_ = 0;
    initialEnergy = saddle_->potentialEnergy();

    fprintf(stdout, "  Saddle point search started.\n");

    // either an initial displacement is performed and the search is started
    // or a series of jumps is performed to reach a convex region 
    if (parameters_->saddleRefine) {
        lowestEigenmode_->startNewSearchAndCompute(saddle_, mode);
        eigenValue_ = lowestEigenmode_->returnLowestEigenmode(eigenMode_);
    }else{
        if(parameters_->saddleMaxJumpAttempts <= 0){
            displaceInConcaveRegion();
        }else{
            jumpToConvexRegion();
        }
    }
    fprintf(stdout, "  Saddle point displaced.\n");

    if(status_ == statusInit)
       searchForSaddlePoint(initialEnergy);
        
    if(status_ == statusInit){
        fprintf(stdout, "    Saddle point determined.\n");        
        relaxFromSaddle(min1, min2);
        fprintf(stdout, "    Minima determined.\n");
    }
    return(status_);
}

long SaddlePoint::getnFreeCoord() const
{
    return nFreeCoord_;
}

double const *const SaddlePoint::getEigenMode() const
{
    return eigenMode_;
}

void SaddlePoint::displaceState(Matter *matter)
{
    long nAtoms = saddle_->numberOfAtoms();
    long j, indexEpiCenter = 0;
    double diffR;

    double *pos;
    //RT: commented following line; replaced with initialDisplacement_
    //double *displacement;
    
    pos = new double[3 * nAtoms];
    initialDisplacement_ = new double[3 * nAtoms];
    
    for(int i=0; i<nAtoms; i++)
    {
        initialDisplacement_[3 * i + 0] = 0;
        initialDisplacement_[3 * i + 1] = 0;
        initialDisplacement_[3 * i + 2] = 0;
    }
    
    if(parameters_->saddleTypePerturbation == dispNotFccOrHcp)
    {
        indexEpiCenter = EpiCenters::cnaEpiCenter(matter);            
    }
    else if(parameters_->saddleTypePerturbation == dispLastAtom)
    {
        indexEpiCenter = EpiCenters::lastAtom(matter);
    }
    else if(parameters_->saddleTypePerturbation == dispMinCoordinated)
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
            if(diffR < parameters_->saddleWithinRadiusPerturbated)
            {
                initialDisplacement_[3 * i + 0] = 2 * randomDouble() - 1;
                initialDisplacement_[3 * i + 1] = 2 * randomDouble() - 1;
                initialDisplacement_[3 * i + 2] = 2 * randomDouble() - 1;
            }
        }
        j++;
    }
    normalize(initialDisplacement_, 3 * nAtoms);
    multiplyScalar(initialDisplacement_, initialDisplacement_, parameters_->saddleNormPerturbation, 3 * nAtoms);
 
    for(int i = 0; i < 3 * nAtoms; i++)
    {
        if(parameters_->saddleMaxSinglePerturbation < initialDisplacement_[i])
        {
            initialDisplacement_[i] = parameters_->saddleMaxSinglePerturbation;
        }
        else if(initialDisplacement_[i] < -parameters_->saddleMaxSinglePerturbation)
        {
            initialDisplacement_[i] = -parameters_->saddleMaxSinglePerturbation;
        }
    }
    // Adding the initialDisplacement_
    matter->getFreePositions(pos);
    add(pos, pos, initialDisplacement_, 3*nAtoms);
    matter->setFreePositions(pos);
 
    delete [] pos;
    //RT: commented following line.
    //delete [] displacement;
    return;
}

void SaddlePoint::correctingForces(double *force){
 
    double *tempDoubleList;
 
    tempDoubleList = new double[nFreeCoord_];
    //----- Initialize end -----
    //std::cout<<"correctingForces\n";
 
    makeProjection(tempDoubleList, force, eigenMode_, nFreeCoord_);
 
    if (0 < eigenValue_){
        if (parameters_->saddlePerpendicularForceRatio > 0.0) {
            // reverse force parallel to eigenvector, and reduce perpendicular force
            double const d=parameters_->saddlePerpendicularForceRatio;
            multiplyScalar(force, force, d, nFreeCoord_);
            multiplyScalar(tempDoubleList, tempDoubleList, -1-d, nFreeCoord_);
            add(force, force, tempDoubleList, nFreeCoord_);
        }
        else {
            // Follow eigenmode
            multiplyScalar(force, tempDoubleList, -1, nFreeCoord_);
        }
    }
    else{
        // Reversing force parallel to eigenmode
        multiplyScalar(tempDoubleList, tempDoubleList, -2, nFreeCoord_);
        add(force, force, tempDoubleList, nFreeCoord_);
    }
    delete [] tempDoubleList;
    return;
}

void SaddlePoint::relaxFromSaddle(Matter *min1, Matter *min2){
    double *posSaddle;
    double *displacedPos;
 
    posSaddle = new double[nFreeCoord_];
    displacedPos = new double[nFreeCoord_];
 
    saddle_->getFreePositions(posSaddle);
    //----- Initialize end -----
    //std::cout<<"relaxFromSaddle\n";
 
    // Displace saddle point configuration along the lowest eigenmode and minimize
    *min1 = *saddle_;
    //XXX: the distance displced from the saddle should be a parameter
    multiplyScalar(displacedPos, eigenMode_, 0.2, nFreeCoord_);
    // NOTE using subtract
    subtract(displacedPos, posSaddle, displacedPos, nFreeCoord_);
    min1->setFreePositions(displacedPos);
    ConjugateGradients cgMin1(min1, parameters_);
    cgMin1.fullRelax();
 
    *min2 = *saddle_;
    multiplyScalar(displacedPos, eigenMode_, 0.2, nFreeCoord_);
    // NOTE using add
    add(displacedPos, posSaddle, displacedPos, nFreeCoord_);
    min2->setFreePositions(displacedPos);
    ConjugateGradients cgMin2(min2, parameters_);  
    cgMin2.fullRelax();
 
    delete [] posSaddle;
    delete [] displacedPos;
    return;
}

void SaddlePoint::addForceCallsSaddlePoint(long fcalls, double eigenvalue){
    if(0 < eigenvalue)
        forceCallsSaddlePointConcave_ += fcalls;
    else
        forceCallsSaddlePointConvex_ += fcalls;
    return;
}

void SaddlePoint::jumpToConvexRegion(){
    long forceCallsSaddle;
    long iterations = 0;
 
    double *pos;
    pos = new double[nFreeCoord_];
 
    saddle_->getFreePositions(pos);
    forceCallsSaddle = saddle_->getForceCalls();
    //----- Initialize end -----
    //std::cout<<"jumpToConvexRegion\n";

    if(parameters_->saddleTypePerturbation!=dispNone){
        do{
            saddle_->setFreePositions(pos);
            displaceState(saddle_);
            lowestEigenmode_->startNewSearchAndCompute(saddle_, initialDisplacement_);
            eigenValue_ = lowestEigenmode_->returnLowestEigenmode(eigenMode_);
            iterations = iterations+1;
        }while((0 < eigenValue_) && 
               (iterations < parameters_->saddleMaxJumpAttempts));
    }
    if(0 < eigenValue_)
        status_ = statusBadNoConvex;

    forceCallsSaddle = saddle_->getForceCalls()-forceCallsSaddle;
    addForceCallsSaddlePoint(forceCallsSaddle, eigenValue_);

    delete [] pos;
    return;
}

void SaddlePoint::displaceInConcaveRegion()
{
    displaceState(saddle_);
    lowestEigenmode_->startNewSearchAndCompute(saddle_, initialDisplacement_);
    eigenValue_ = lowestEigenmode_->returnLowestEigenmode(eigenMode_);
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

    forcesStep = new double [nFreeCoord_];
    posStep = new double [nFreeCoord_];
    forces = new double [nFreeCoord_];
    pos = new double[nFreeCoord_];

    saddle_->getFreePositions(pos);
    //----- Initialize end -----
    //std::cout<<"searchForSaddlePoint\n";
    saddle_->getFreeForces(forces);
 
    eigenValue_ = lowestEigenmode_->returnLowestEigenmode(eigenMode_);
    correctingForces(forces);
    ConjugateGradients cgSaddle(saddle_, parameters_, forces);
#ifndef NDEBUG
    saddle_->matter2xyz("climb", false);
#endif
    do
    {
        forceCallsSaddle = saddle_->getForceCalls();        
        cgSaddle.makeInfinitesimalStepModifiedForces(posStep, pos);
        // Determining the optimal step
        saddle_->setFreePositions(posStep);
        saddle_->getFreeForces(forcesStep);
        correctingForces(forcesStep);
/*        if(0 < eigenValue_)
            maxStep = parameters_->getMaxStepSizeConcave_SP();
        else
            maxStep = parameters_->getMaxStepSizeConvex_SP();*/
        maxStep = parameters_->saddleMaxStepSize;
        cgSaddle.getNewPosModifiedForces(pos, forces, forcesStep, maxStep);
        // The system (saddle_) is moved to a new configuration
        saddle_->setFreePositions(pos);
        saddle_->getFreeForces(forces);
        // The new lowest eigenvalue
        lowestEigenmode_->moveAndCompute(saddle_);
        eigenValue_ = lowestEigenmode_->returnLowestEigenmode(eigenMode_);
        // Updating the conjugated object to the new configuration
        correctingForces(forces);
        cgSaddle.setFreeAtomForcesModifiedForces(forces);
        if(eigenValue_ < 0)
        {
            converged = cgSaddle.isItConverged(parameters_->saddleConverged);
            concaveSeries = 0;
        }
        else
        {
            converged = false;
            concaveSeries = concaveSeries + 1;
        }
        forceCallsSaddle = saddle_->getForceCalls()-forceCallsSaddle;        
        addForceCallsSaddlePoint(forceCallsSaddle, eigenValue_);

        iterations++;
        #ifndef NDEBUG
            printf("climb = %ld, max force = %f\n", iterations, saddle_->maxForce());
            saddle_->matter2xyz("climb", true);
        #endif
        energySaddle = saddle_->potentialEnergy();
    }while(!converged && 
           (iterations < parameters_->saddleMaxIterations) && 
           (energySaddle-initialEnergy < parameters_->saddleMaxEnergy));
    if(!converged){
        if(parameters_->saddleMaxIterations <= iterations) {
            status_ = statusBadMaxIterations;
        }
        if(parameters_->saddleMaxEnergy <= energySaddle-initialEnergy) {
            status_ = statusBadHighBarrier;
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
      return lowestEigenmode_;
}
