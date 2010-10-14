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

Dimer::Dimer(Matter const *matter, Parameters *params)
{
    long nAllCoord;
    parameters     = params;
    matterInitial  = new Matter(parameters);
    matterDimer    = new Matter(parameters);
    *matterInitial = *matter;
    *matterDimer   = *matter;
    nAllCoord   = 3 * matter->numberOfAtoms();
    nFreeCoord = 3 * matter->numberOfFreeAtoms();  
    nAtoms = matter->numberOfAtoms();
    
    directionNorm.resize(nAtoms, 3);
    rotationalPlaneNorm.resize(nAtoms, 3);
    directionNorm.setZero();
    rotationalPlaneNorm.setZero();
    totalForceCalls = 0;
}


Dimer::~Dimer()
{
    delete matterInitial;
    delete matterDimer;
}


void Dimer::moveAndCompute(Matter const *matter)
{
    *matterInitial = *matter;
    estimateLowestEigenmode(parameters->dimerRotations);
    return;
}


void Dimer::startNewSearchAndCompute(Matter const *matter, Matrix<double, Eigen::Dynamic, 3> displacement)
{
    *matterInitial = *matter;
    
    rotationalPlaneNorm.setZero();
    // Create an initial direction for the dimer

    directionNorm = displacement.cwise() * matter->getFree();
    directionNorm.normalize();
    
    // RT: following part of the algorithm has been replaced with the torque window.
    // The initial search should be allowed to use extra iterations to obtain
    // better estimate
    //estimateLowestEigenmode(parameters->getRotationsNewSearch_Dimer());
}


void Dimer::estimateLowestEigenmode(long rotationsToPerform)
{
    long rotations = 0;
    long forceCallsInitial;
    long forceCallsDimer;
    double forceDimer1AlongRotationalPlaneNorm;
    double forceDimer2AlongRotationalPlaneNorm;
    double curvature, rotationalForceChange, forceDimer, rotationAngle;
    double lengthRotationalForceOld;
    double torqueMagnitude = 0.0;
    bool doneRotating = false;
    
    rotationalForceChange = forceDimer = rotationAngle = curvature = 0;
    forceDimer1AlongRotationalPlaneNorm = 0;
    forceDimer2AlongRotationalPlaneNorm = 0;
    Matrix<double, Eigen::Dynamic, 3> rotationalForce(nAtoms,3);
    Matrix<double, Eigen::Dynamic, 3> rotationalForceOld(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> rotationalPlaneNormOld(nAtoms, 3);
    rotationalForce.setZero();
    rotationalForceOld.setZero();
    rotationalPlaneNormOld.setZero();
    

    lengthRotationalForceOld = 0;
    forceCallsInitial = matterInitial->getForceCalls();
    forceCallsDimer = matterDimer->getForceCalls();
    // Uses two force calls per rotation
    while(!doneRotating){
        // First dimer
        //cout<<"Before "<<rotationalForce.norm()<<endl;
        curvature = calcRotationalForce(rotationalForce);  
        //cout<<"After "<<rotationalForce.norm()<<endl;

        
        // The new rotational plane is determined
        
        //cout<<"rfo Before "<<rotationalForceOld.norm()<<endl;
        determineRotationalPlane(rotationalForce, 
                                 rotationalForceOld, 
                                 rotationalPlaneNormOld,
                                 &lengthRotationalForceOld);
                                 
        //cout<<"rfo After "<<rotationalForceOld.norm()<<endl;
        // Calculate the magnitude of the torque on the dimer.
        torqueMagnitude = rotationalForce.squaredNorm();

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
        forceDimer1AlongRotationalPlaneNorm = (rotationalForce.cwise()*rotationalPlaneNorm).sum();

        rotateDimerAndNormalizeAndOrthogonalize(parameters->dimerRotationAngle);
        
        if(!doneRotating)
        {
            // Second dimer
            curvature = calcRotationalForce(rotationalForce);
            
            forceDimer2AlongRotationalPlaneNorm = (rotationalForce.cwise()*rotationalPlaneNorm).sum();
            
            rotationalForceChange = ((forceDimer1AlongRotationalPlaneNorm - forceDimer2AlongRotationalPlaneNorm) / 
                                     parameters->dimerRotationAngle);
            
            forceDimer = (forceDimer1AlongRotationalPlaneNorm + forceDimer2AlongRotationalPlaneNorm) / 2;
            
            rotationAngle = (atan(2 * forceDimer / rotationalForceChange) / 2 - parameters->dimerRotationAngle / 2);
            
            if(rotationalForceChange < 0)
            {
                rotationAngle = rotationAngle + M_PI / 2;
            }
                
            rotateDimerAndNormalizeAndOrthogonalize(rotationAngle);
            
            rotationalPlaneNormOld = rotationalPlaneNorm; //XXX: Is this copying correctly???        
    
            rotations++;
        }
    }    
    eigenvalue = curvature;

    forceCallsInitial = matterInitial->getForceCalls()-forceCallsInitial;
    forceCallsDimer = matterDimer->getForceCalls()-forceCallsDimer;

    totalForceCalls += forceCallsInitial+forceCallsDimer;

    return;
}


double Dimer::getEigenvalue(){
    return eigenvalue;
}

void Dimer::setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector)
{
    directionNorm = eigenvector;
    eigenvalue=0.0;
}

Matrix<double, Eigen::Dynamic, 3> Dimer::getEigenvector()
{
      return directionNorm;
}

double Dimer::calcRotationalForce(Matrix<double, Eigen::Dynamic, 3> &rotationalForce){
    
    double projectedForceA, projectedForceB;
    Matrix<double, Eigen::Dynamic, 3> posInitial(nAtoms,3);
    Matrix<double, Eigen::Dynamic, 3> posDimer(nAtoms,3);
    Matrix<double, Eigen::Dynamic, 3> forceInitial(nAtoms,3);
    Matrix<double, Eigen::Dynamic, 3> forceA(nAtoms,3);
    Matrix<double, Eigen::Dynamic, 3> forceB(nAtoms,3);
   
    posInitial = matterInitial->getPositions();    

    // Displacing the one of the dimer configurations
    posDimer = posInitial + directionNorm*parameters->dimerSeparation;
    // Obtaining the force for configuration A
    matterDimer->setPositions(posDimer);
    forceA = matterDimer->getForces();
    
    // Use "forward differencing" together with the central averaging formula
    // to obtain the force for configuration B
    forceInitial = matterInitial->getForces();
    forceB = 2*forceInitial - forceA;
    
    //cout<<"force norm: "<<forceA.norm()<<" "<<forceB.norm()<<endl;
    projectedForceA = (directionNorm.cwise() * forceA).sum();
    projectedForceB = (directionNorm.cwise() * forceB).sum();
    //cout<<"projected force: "<<projectedForceA<<" "<<projectedForceB<<endl;
    // Remove force component parallel to dimer
    forceA = makeOrthogonal(forceA, directionNorm);
    forceB = makeOrthogonal(forceB, directionNorm);
    
    // Determine difference in force orthogonal to dimer
    rotationalForce = (forceA - forceB)/parameters->dimerSeparation;
    //XXX: Side effect

    

    // Based on difference in force parallel to dimer    
    return (projectedForceB-projectedForceA)/(2*parameters->dimerSeparation);
}


void Dimer::determineRotationalPlane(Matrix<double, Eigen::Dynamic, 3> rotationalForce, 
                                     Matrix<double, Eigen::Dynamic, 3> &rotationalForceOld, 
                                     Matrix<double, Eigen::Dynamic, 3> rotationalPlaneNormOld,
                                     double* lengthRotationalForceOld){
    double a, b, gamma = 0;
    
    a = fabs((rotationalForce.cwise() * rotationalForceOld).sum());
    b = rotationalForceOld.squaredNorm();
    //cout<<"a "<<a<<" b "<<b<<endl;
    if(a<0.5*b)
    {
        gamma = (rotationalForce.cwise() * (rotationalForce - rotationalForceOld)).sum()/b;
    }
    else
        gamma = 0;
   //cout<<"gamma: "<<gamma<<endl; 
    // The new rotational plane
    // Based on the current rotational force and the rotational plane force 
    // from the former iteration    
    rotationalPlaneNorm = rotationalForce + rotationalPlaneNormOld * (*(lengthRotationalForceOld)) * gamma;
    
    // The planes normal is normalized,
    // made orthogonal to the dimers direction and renormalized
    
    *lengthRotationalForceOld = rotationalPlaneNorm.norm();
    //cout<<"lRFO "<<*lengthRotationalForceOld<<endl;
    rotationalPlaneNorm.normalize(); 
    rotationalPlaneNorm = makeOrthogonal(rotationalPlaneNorm, directionNorm); 
    rotationalPlaneNorm.normalize(); 

    rotationalForceOld = rotationalForce;
    
    return;
}


void Dimer::rotateDimerAndNormalizeAndOrthogonalize(double rotationAngle)
{
    double cosAngle, sinAngle;

    cosAngle = cos(rotationAngle);
    sinAngle = sin(rotationAngle);
   
    directionNorm = directionNorm * cosAngle + rotationalPlaneNorm * sinAngle;
    rotationalPlaneNorm = rotationalPlaneNorm * cosAngle - directionNorm * sinAngle; 
    directionNorm.normalize();
    rotationalPlaneNorm.normalize(); 
    // Remove component from rotationalPlaneNorm parallel to directionNorm
    rotationalPlaneNorm = makeOrthogonal(rotationalPlaneNorm, directionNorm);
    rotationalPlaneNorm.normalize();
    return;
}
