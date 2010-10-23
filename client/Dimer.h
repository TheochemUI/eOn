/*
 *===============================================
 *  EON Dimer.h
 *===============================================
 */
#ifndef DIMER_H
#define DIMER_H

#include <math.h>
#include <cmath>
#include <cassert>
#include "debug.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "LowestEigenmodeInterface.h"

#define PI 3.141592653589793

/** Use the dimer method to estimate the lowest eigenvector of the Hessian and its corresponding eigenvalue (curvature). In this version a 'fictive' dimer is used, only one dimer configuration is calculated explicitly whereas the other configuration is estimated with forward differencing relative to the center of the dimer*/
class Dimer : public LowestEigenmodeInterface{
public:
    /** Constructor where object is initialized.
    @param[in]      *matter                    Pointer to Matter object. Only used to initialize the dimer object. The pointer is not stored.
    @param[in]      *parameters                Pointer to the Parameter object containing the runtime parameters.*/    
    Dimer(Matter const *matter, Parameters *parameters);
    
    ~Dimer();///< Destructor.
        
    void startNewSearchAndCompute(Matter const *matter, Matrix<double, Eigen::Dynamic, 3>);///< Is computational heavy. The matter object pointed at is copied into a local copy defining the center of the dimer.
    void moveAndCompute(Matter const *matter);///< Is computational heavy! The matter object pointed at is copied into a local copy defining the center of the dimer. As the function is intended to be used for moving the dimer up the potential energy surface, matter should be somehow alike the argument used in the former call to Dimer::moveAndCompute or Dimer::startNewSearchAndCompute.
    double getEigenvalue();///< The latest determined lowest eigenmode is returned in result, the array should have a length corrsponding to the number free coordinates.
    /** Set initial direction manually.*/
    void setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector);
        /// Return eigenvector.
    Matrix<double, Eigen::Dynamic, 3>  getEigenvector();
    long totalForceCalls;
private:
    Matter *matterInitial;///< Used for the center of the dimer.
    Matter *matterDimer;///< Used for the one of the configurations defining the dimer.
    Matrix<double, Eigen::Dynamic, 3> directionNorm;///< Direction descring how matterDimer_ is displaced relative to matterInitial_.
    Matrix<double, Eigen::Dynamic, 3> rotationalPlaneNorm;///< Used to store the normal of plane in which the dimer is ratated.
    double eigenvalue;///< An estimate for the lowest eigenvale, its most important feature is wheter it is negative (concave region) or positive (convex region).
    long nFreeCoord;///< Number of free coordinates.
    int nAtoms;
    Parameters *parameters;///< Pointer to the runtime parameters. Note that the structure is used to store how many force calls that was used in the locating the eigenmodes.
    void estimateLowestEigenmode();///< Is computational heavy! Try to obtain a converged result for the eigenmode within the number of rotations passed in \a maxRotations.

    /** The rotational plane that is going to be used is determined with the Conjugate Gradient method.
    @param[in]      *rotationalForce           Pointer to array of length nFreeCoord_. The rotational forces just calculated.
    @param[in,out]  *rotationalForceOld        Pointer to array of length nFreeCoord_. The rotational forces from the former iteration. The rotational forces just calculated has been copied into the array, to be used in the next iteration.
    @param[in]      *rotationalPlaneNormOld    Pointer to array of length nFreeCoord_. The rotational plane from the former iteration.      
    @param[in,out]  *lengthRotationalForceOld  Pointer to value. The length of the rotational forces from the former iteration. The length of the rotational forces just calcuated, to be used in the next iteration*/
    void determineRotationalPlane(Matrix<double, Eigen::Dynamic, 3> rotationalForce, 
                                  Matrix<double, Eigen::Dynamic, 3> &rotationalForceOld, 
                                  Matrix<double, Eigen::Dynamic, 3> rotationalPlaneNormOld,
                                  double *lengthRotationalForceOld);
    
    void rotateDimerAndNormalizeAndOrthogonalize(double rotationAngle);///< Rotating the dimer rotationAngle (radians), the rotational plane is also modified.
    double calcRotationalForce(Matrix<double, Eigen::Dynamic, 3> &forceDiffOrthognalToDimer);///< Determine the rotational force of the dimer. forceDiffOrthognalToDimer should be the difference in forces for the 'two' configurations defining the dimer where the forces along the dimer directionNorm_ has been projected out.
};
#endif
