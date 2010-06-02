/*
 *===============================================
 *  Created by Andreas Pedersen on 10/5/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *-----------------------------------------------
 *  Heavily inspired of codes by:
 *      Graeme Henkelmann
 *      Andri Arnaldsson
 *      Roar Olsen
 *===============================================
 */
#include "Dimer.h"

using namespace helper_functions;
using namespace constants;

Dimer::Dimer(Matter const *matter, Parameters *parameters){
    long nAllCoord;

    matterInitial_ = new Matter();
    matterDimer_ = new Matter();
    *matterInitial_ = *matter;
    *matterDimer_ = *matter;

    nAllCoord = 3*matter->numberOfAtoms();
    nFreeCoord_ = 3*matter->numberOfFreeAtoms();    
    
    // There should be space for both free and frozen atoms
    tempListDouble_ = new double[nAllCoord];
    
    directionNorm_ = new double[nFreeCoord_];
    rotationalPlaneNorm_ = new double[nFreeCoord_];
    
    parameters_ = parameters;
}

Dimer::~Dimer(){
    delete [] directionNorm_;
    delete [] rotationalPlaneNorm_;
    delete [] tempListDouble_;

    delete matterInitial_;
    delete matterDimer_;
}

void Dimer::moveAndCompute(Matter const *matter){
    *matterInitial_ = *matter;
//    determineLowestEigenmode(parameters_->getMaxRotationIterations_Dimer());
    estimateLowestEigenmode(parameters_->getMaxRotationIterations_Dimer());
    return;
}

void Dimer::startNewSearchAndCompute(Matter const *matter){
    *matterInitial_ = *matter;

    // create a initial direction for the dimer
    for(int i=0; i<nFreeCoord_; i++){
        directionNorm_[i] = random()-0.5;
        rotationalPlaneNorm_[i] = 0;
    }
    normalize(directionNorm_, nFreeCoord_);
    
    // The initial search should be allowed to use extra iterations to converge
//    determineLowestEigenmode(parameters_->getMaxRotationIterations_Dimer()+
//                             getDimerExtraRotationsInitial());
    estimateLowestEigenmode(parameters_->getMaxRotationIterations_Dimer()+
                            getDimerExtraRotationsInitial());

}

void Dimer::estimateLowestEigenmode(long rotationsToPerform){
//bool Dimer::determineLowestEigenmode(long maxRotations){
//    bool convergedBool = false;
    
    long rotations = 0;
    long forceCallsInitial;
    long forceCallsDimer;
    
    double forceDimer1AlongRotationalPlaneNorm;
    double forceDimer2AlongRotationalPlaneNorm;
    double curvature, rotationalForceChange, forceDimer, rotationAngle;
    double *rotationalForce, *rotationalForceOld, *rotationalPlaneNormOld;
    double lengthRotationalForceOld;

    rotationalForceChange = forceDimer = rotationAngle = curvature = 0;
    forceDimer1AlongRotationalPlaneNorm = 0;
    forceDimer2AlongRotationalPlaneNorm = 0;

    rotationalForce = new double[nFreeCoord_];
    rotationalForceOld = new double[nFreeCoord_];
    rotationalPlaneNormOld = new double[nFreeCoord_];

    // Values to ensure that gamma equals zero,
    // (a<b/2) is false in the part where the rotational plane is determined
    for(int i=0;i<nFreeCoord_;i++){
        rotationalForce[i] = 0;
        rotationalForceOld[i] = 0;
    } 
    lengthRotationalForceOld = 0;
    //----- Initialize end -----
    //std::cout<<"determineLowestEigenmode\n";
    //std::cout<<"estimateLowestEigenmode\n";

    forceCallsInitial = matterInitial_->getForceCalls();
    forceCallsDimer = matterDimer_->getForceCalls();
    
    while(rotations<rotationsToPerform){

        // First dimer
        curvature = calcRotationalForce(rotationalForce);        
        
//        // The calculation has converged 
//        // and the minimum number of rotations has been performed
//        if((length(rotationalForce,nFreeCoord_) < 
//            parameters_->getRotationalForceConverged_Dimer())
//           and (getDimerMinRotations() < rotations)){
//            convergedBool = true;
//            break;
//        }
        // The new rotational plane is determined
        determineRotationalPlane(rotationalForce, 
                                 rotationalForceOld, 
                                 rotationalPlaneNormOld,
                                 &lengthRotationalForceOld);
                
        // Rotational force along the rotational planes normal
        forceDimer1AlongRotationalPlaneNorm = dot(rotationalForce, 
                                                  rotationalPlaneNorm_, 
                                                  nFreeCoord_);

        rotateDimerAndNormalizeAndOrthogonalize(getDimerRotationAngle());
        
        // Second dimer
        curvature = calcRotationalForce(rotationalForce);
        
        forceDimer2AlongRotationalPlaneNorm = dot(rotationalForce, 
                                                  rotationalPlaneNorm_, 
                                                  nFreeCoord_);
        
        rotationalForceChange = ((forceDimer1AlongRotationalPlaneNorm-
                                 forceDimer2AlongRotationalPlaneNorm)/
                                 getDimerRotationAngle());
        
        forceDimer = (forceDimer1AlongRotationalPlaneNorm+
                      forceDimer2AlongRotationalPlaneNorm)/2;
        
        rotationAngle = (atan(2*forceDimer/rotationalForceChange)/2-
                         getDimerRotationAngle()/2);
        
        if(rotationalForceChange < 0)
            // M_PI is from the math.h library
            rotationAngle = rotationAngle+M_PI/2;
            
        rotateDimerAndNormalizeAndOrthogonalize(rotationAngle);
        
        copyRightIntoLeft(rotationalPlaneNormOld, rotationalPlaneNorm_, 
                          nFreeCoord_);
        rotations++;
    }
    forceCallsInitial = matterInitial_->getForceCalls()-forceCallsInitial;
    forceCallsDimer = matterDimer_->getForceCalls()-forceCallsDimer;
    parameters_->addForceCallsSaddlePoint(forceCallsInitial+forceCallsDimer);
    
    eigenvalue_ = curvature;
    delete [] rotationalForce;
    delete [] rotationalForceOld;
    delete [] rotationalPlaneNormOld;
//    return convergedBool;
    return;
}

double Dimer::returnLowestEigenmode(double *result){
    for(int i=0;i<nFreeCoord_;i++)
        result[i]=directionNorm_[i];
    return eigenvalue_;
}

double Dimer::calcRotationalForce(double *rotationalForce){
    
    double *posInitial;
    double *posDimer;
    double *forceA;
    double *forceB;
    double *forceInitial;
    double projectedForceA, projectedForceB;
    posInitial = new double[nFreeCoord_];
    posDimer = new double[nFreeCoord_];
    forceInitial = new double[nFreeCoord_];
    forceA = new double[nFreeCoord_];
    forceB = new double[nFreeCoord_];
    
    matterInitial_->getFreePositions(posInitial);    
    //----- Initialize end -----
    //std::cout<<"calcRotationalForce\n";

    // Displacing the one of the dimer configurations
    multiplyScalar(tempListDouble_,directionNorm_,getDimerSize(),nFreeCoord_);
    add(posDimer, posInitial, tempListDouble_, nFreeCoord_);

    // Obtaining the force for picture A
    matterDimer_->setFreePositions(posDimer);
    matterDimer_->getFreeForces(forceA);
    
    // Use "forward differencing" together with the central averaging formula
    // to obtain the force for configuration B
    matterInitial_->getFreeForces(forceInitial);
    multiplyScalar(tempListDouble_, forceInitial, 2, nFreeCoord_);
    subtract(forceB, tempListDouble_, forceA, nFreeCoord_);
    
    projectedForceA = dot(directionNorm_, forceA, nFreeCoord_);
    projectedForceB = dot(directionNorm_, forceB, nFreeCoord_);

    // Remove force component parallel to dimer
    makeOrthogonal(forceA, forceA, directionNorm_, nFreeCoord_);
    makeOrthogonal(forceB, forceB, directionNorm_, nFreeCoord_);
    
    // Determine difference in force orthogonal to dimer
    subtract(tempListDouble_, forceA, forceB, nFreeCoord_);
    divideScalar(rotationalForce, tempListDouble_, getDimerSize(), nFreeCoord_);
    
    delete [] posInitial;    
    delete [] posDimer;
    delete [] forceInitial;
    delete [] forceA;
    delete [] forceB;

    // Based on difference in force parallel to dimer    
    return (projectedForceB-projectedForceA)/(2*getDimerSize());
}

void Dimer::determineRotationalPlane(double *rotationalForce, 
                                     double *rotationalForceOld, 
                                     double *rotationalPlaneNormOld,
                                     double *lengthRotationalForceOld){
    double a, b, gamma = 0;
    double *rotationalPlane;
    rotationalPlane = new double[nFreeCoord_];
    //----- Initialize end -----
    //std::cout<<"determineRotationalPlane\n";

    
    a = fabs(dot(rotationalForce, rotationalForceOld, nFreeCoord_));
    b = dot(rotationalForceOld, rotationalForceOld, nFreeCoord_);
    if(a<0.5*b){
        subtract(tempListDouble_, 
                 rotationalForce, 
                 rotationalForceOld, 
                 nFreeCoord_);
        
        //Polak-Ribiere way to determine how much to mix in of old direction
        gamma = dot(rotationalForce, tempListDouble_, nFreeCoord_)/b;  
    }
    else
        gamma = 0;
    
    // The new rotational plane
    // Based on the current rotational force and the rotational plane force 
    // from the former iteration    
    multiplyScalar(tempListDouble_, rotationalPlaneNormOld, 
                   *(lengthRotationalForceOld)*gamma, nFreeCoord_);    
    
    add(rotationalPlane, rotationalForce, tempListDouble_, nFreeCoord_);    
    
    // The planes normal is normalized,
    // made orthogonal to the dimers direction and renormalized
    copyRightIntoLeft(rotationalPlaneNorm_, rotationalPlane, nFreeCoord_);
    
    normalize(rotationalPlaneNorm_, nFreeCoord_);
    
    makeOrthogonal(rotationalPlaneNorm_, rotationalPlaneNorm_,
                   directionNorm_, nFreeCoord_);
    
    normalize(rotationalPlaneNorm_, nFreeCoord_);

    copyRightIntoLeft(rotationalForceOld, rotationalForce, nFreeCoord_);
    *lengthRotationalForceOld = length(rotationalPlane, nFreeCoord_);
    
    delete [] rotationalPlane;
    return;
}

void Dimer::rotateDimerAndNormalizeAndOrthogonalize(double rotationAngle){
    double temp1, temp2, cosAngle, sinAngle;
    //----- Initialize end -----
    //std::cout<<"rotateDimerAndNormalizeAndOrthogonalize\n";

    cosAngle = cos(rotationAngle);
    sinAngle = sin(rotationAngle);
    
    for(int i=0; i<nFreeCoord_; i++){
        temp1 = directionNorm_[i]*cosAngle + rotationalPlaneNorm_[i]*sinAngle;
        temp2 = rotationalPlaneNorm_[i]*cosAngle - directionNorm_[i]*sinAngle;
        directionNorm_[i] = temp1;
        rotationalPlaneNorm_[i] = temp2;
    }
    
    normalize(directionNorm_, nFreeCoord_);
    normalize(rotationalPlaneNorm_, nFreeCoord_);
    
    // Remove component from rotationalPlaneNorm_ parallel to directionNorm_
    makeOrthogonal(rotationalPlaneNorm_, rotationalPlaneNorm_, 
                   directionNorm_, nFreeCoord_);
    normalize(rotationalPlaneNorm_, nFreeCoord_);

    return;
}
