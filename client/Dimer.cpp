/*
 *===============================================
 *  EON Dimer.cpp
 *-----------------------------------------------
 *  Todo:
 *      GH: add improved dimer method of Andreas Heyden
 *          - add 45 degree rotation angle
 *          - add rotation angle instead of rotational force criteria
 *      GH: add LBFGS optimizer
 *      GH: make torqueLimitHigh/Low and torqueMin/MaxRotations variables
 *===============================================
 */

#include "Dimer.h"

using namespace helper_functions;

Dimer::Dimer(Matter const *matter, Parameters *parameters)
{
    long nAllCoord;
    parameters_     = parameters;
    matterInitial_  = new Matter(parameters_);
    matterDimer_    = new Matter(parameters_);
    *matterInitial_ = *matter;
    *matterDimer_   = *matter;
    nAllCoord   = 3 * matter->numberOfAtoms();
    nFreeCoord_ = 3 * matter->numberOfFreeAtoms();    
    tempListDouble_ = new double[nAllCoord]; // There should be space for both free and frozen atoms.
    directionNorm_ = new double[nFreeCoord_];
    rotationalPlaneNorm_ = new double[nFreeCoord_];

    totalForceCalls = 0;
}


Dimer::~Dimer()
{
    delete [] directionNorm_;
    delete [] rotationalPlaneNorm_;
    delete [] tempListDouble_;
    delete matterInitial_;
    delete matterDimer_;
}


void Dimer::moveAndCompute(Matter const *matter)
{
    *matterInitial_ = *matter;
    estimateLowestEigenmode(parameters_->dimerRotations);
    return;
}


void Dimer::startNewSearchAndCompute(Matter const *matter, double *displacement)
{
    *matterInitial_ = *matter;
    long nAtoms = matter->numberOfAtoms();
    long index = 0;
        
    // Create an initial direction for the dimer
    for(int i = 0; i < nAtoms; i++)
    {
        if(!matter->getFixed(i))
        {
            directionNorm_[index + 0] = displacement[i * 3 + 0];
            directionNorm_[index + 1] = displacement[i * 3 + 1];
            directionNorm_[index + 2] = displacement[i * 3 + 2];
            rotationalPlaneNorm_[index + 0] = 0;
            rotationalPlaneNorm_[index + 1] = 0;
            rotationalPlaneNorm_[index + 2] = 0;
            index += 3;
        }
    }
    normalize(directionNorm_, nFreeCoord_);
    
    // RT: following part of the algorithm has been replaced with the torque window.
    // The initial search should be allowed to use extra iterations to obtain
    // better estimate
    //estimateLowestEigenmode(parameters_->getRotationsNewSearch_Dimer());
}


void Dimer::estimateLowestEigenmode(long rotationsToPerform)
{
    long rotations = 0;
    long forceCallsInitial;
    long forceCallsDimer;
    double forceDimer1AlongRotationalPlaneNorm;
    double forceDimer2AlongRotationalPlaneNorm;
    double curvature, rotationalForceChange, forceDimer, rotationAngle;
    double *rotationalForce, *rotationalForceOld, *rotationalPlaneNormOld;
    double lengthRotationalForceOld;
    double torqueMagnitude = 0.0;
    bool doneRotating = false;
    
    rotationalForceChange = forceDimer = rotationAngle = curvature = 0;
    forceDimer1AlongRotationalPlaneNorm = 0;
    forceDimer2AlongRotationalPlaneNorm = 0;
    rotationalForce = new double[nFreeCoord_];
    rotationalForceOld = new double[nFreeCoord_];
    rotationalPlaneNormOld = new double[nFreeCoord_];

    // Values to ensure that gamma equals zero,
    // (a<b/2) is false in the part where the rotational plane is determined
    for(int i = 0;i < nFreeCoord_; i++)
    {
        rotationalForce[i] = 0;
        rotationalForceOld[i] = 0;
        rotationalPlaneNormOld[i] = 0;
    } 
    lengthRotationalForceOld = 0;
    forceCallsInitial = matterInitial_->getForceCalls();
    forceCallsDimer = matterDimer_->getForceCalls();
    // Uses two force calls per rotation
    while(!doneRotating){
        // First dimer
        curvature = calcRotationalForce(rotationalForce);        
        
        // The new rotational plane is determined
        determineRotationalPlane(rotationalForce, 
                                 rotationalForceOld, 
                                 rotationalPlaneNormOld,
                                 &lengthRotationalForceOld);
                                 
        // Calculate the magnitude of the torque on the dimer.
        double sum = 0.0;
        for(int i = 0; i < nFreeCoord_; i++)
        {
            sum += rotationalForce[i] * rotationalForce[i];
        }
        torqueMagnitude = sqrt(sum);

        double torqueLimitHigh = 1.0;
        double torqueLimitLow = 0.1;
        int torqueMaxRotations = 32;
        int torqueMinRotations = 1;
  
        assert(( std::isnormal(torqueMagnitude) ));
        if(torqueMagnitude > torqueLimitHigh && rotations >= torqueMaxRotations)
        {
            doneRotating = true;
        }
        else if(torqueMagnitude < torqueLimitHigh && torqueMagnitude >= torqueLimitLow && rotations >= torqueMinRotations)
        {
            doneRotating = true;
        }
        else if(torqueMagnitude < torqueLimitLow)
        {
            doneRotating = true;
        }
                
        // Rotational force along the rotational planes normal
        forceDimer1AlongRotationalPlaneNorm = dot(rotationalForce, 
                                                  rotationalPlaneNorm_, 
                                                  nFreeCoord_);

        rotateDimerAndNormalizeAndOrthogonalize(parameters_->dimerRotationAngle);
        
        if(!doneRotating)
        {
            // Second dimer
            curvature = calcRotationalForce(rotationalForce);
            
            forceDimer2AlongRotationalPlaneNorm = dot(rotationalForce, rotationalPlaneNorm_, nFreeCoord_);
            
            rotationalForceChange = ((forceDimer1AlongRotationalPlaneNorm - forceDimer2AlongRotationalPlaneNorm) / 
                                     parameters_->dimerRotationAngle);
            
            forceDimer = (forceDimer1AlongRotationalPlaneNorm + forceDimer2AlongRotationalPlaneNorm) / 2;
            
            rotationAngle = (atan(2 * forceDimer / rotationalForceChange) / 2 - parameters_->dimerRotationAngle / 2);
            
            if(rotationalForceChange < 0)
            {
                rotationAngle = rotationAngle + M_PI / 2;
            }
                
            rotateDimerAndNormalizeAndOrthogonalize(rotationAngle);
            
            copyRightIntoLeft(rotationalPlaneNormOld, rotationalPlaneNorm_, nFreeCoord_);        
    
            rotations++;
        }
    }    
    eigenvalue_ = curvature;

    forceCallsInitial = matterInitial_->getForceCalls()-forceCallsInitial;
    forceCallsDimer = matterDimer_->getForceCalls()-forceCallsDimer;

    totalForceCalls += forceCallsInitial+forceCallsDimer;

    delete [] rotationalForce;
    delete [] rotationalForceOld;
    delete [] rotationalPlaneNormOld;
    return;
}


double Dimer::returnLowestEigenmode(double *result){
    for(int i=0; i<nFreeCoord_;i++)
        result[i]=directionNorm_[i];
    return eigenvalue_;
}

void Dimer::setEigenvector(long size, double const eigenvector[])
{
    assert(size == nFreeCoord_);
    for (int i=0; i < size; ++i)
        directionNorm_[i]=eigenvector[i];
    eigenvalue_=0.0;
}

double const * Dimer::getEigenvector(long & size)
{
      size=nFreeCoord_;
      return directionNorm_;
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

    // Displacing the one of the dimer configurations
    multiplyScalar(tempListDouble_,directionNorm_, parameters_->dimerSeparation,nFreeCoord_);
    add(posDimer, posInitial, tempListDouble_, nFreeCoord_);

    // Obtaining the force for configuration A
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
    divideScalar(rotationalForce, tempListDouble_, parameters_->dimerSeparation, nFreeCoord_);
    
    delete [] posInitial;    
    delete [] posDimer;
    delete [] forceInitial;
    delete [] forceA;
    delete [] forceB;

    // Based on difference in force parallel to dimer    
    return (projectedForceB-projectedForceA)/(2*parameters_->dimerSeparation);
}


void Dimer::determineRotationalPlane(double *rotationalForce, 
                                     double *rotationalForceOld, 
                                     double *rotationalPlaneNormOld,
                                     double *lengthRotationalForceOld){
    double a, b, gamma = 0;
    double *rotationalPlane;
    rotationalPlane = new double[nFreeCoord_];
    
    a = fabs(dot(rotationalForce, rotationalForceOld, nFreeCoord_));
    b = dot(rotationalForceOld, rotationalForceOld, nFreeCoord_);
    if(a<0.5*b)
    {
        subtract(tempListDouble_, rotationalForce, rotationalForceOld, nFreeCoord_);
        
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


void Dimer::rotateDimerAndNormalizeAndOrthogonalize(double rotationAngle)
{
    double temp1, temp2, cosAngle, sinAngle;

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
