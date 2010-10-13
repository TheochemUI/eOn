/*
 *===============================================
 *  EON ConjugateGradients.h
 *===============================================
 */
#ifndef CG_H
#define CG_H

#include "Parameters.h"
#include "Matter.h"
#include "MinimizersInterface.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

/** Functionality relying on the conjugate gradients algorithm. The object is capable of minimizing a Matter object or modified forces being passed in.*/
class ConjugateGradients : public MinimizersInterface{

public:
    /** Constructor to be used when a structure is minimized.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    ConjugateGradients(Matter *matter, Parameters *parameters);

    /** Constructor to be used when modified forces are used.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.
    @param[in]   *forces        Double array containing the forces acting on the matter.*/
    ConjugateGradients(Matter *matter, Parameters *parameters, Matrix<double, Eigen::Dynamic, 3> forces);
  
    ~ConjugateGradients();///< Destructor.

    void oneStep(); // do one iteration
    void fullRelax(); // relax the Matter object at the pointer that was passed with the constructor
    bool isItConverged(double convergeCriterion);
    //----Functions used when forces are modified (saddle point search)----
    /** Performs and infinitesimal step along the search direction.
    @param[out]  *posStep  Double array containing the modified positions of the atoms.
    @param[in]   *pos      Double array containing the original positions of the atoms.*/
    Matrix<double, Eigen::Dynamic, 3> makeInfinitesimalStepModifiedForces(Matrix<double, Eigen::Dynamic, 3> pos);

    /** Performs and infinitesimal step along the search direction.
    @param[out]  *pos              Double array containing the calculated new positions of the atoms.
    @param[in]   *forceBeforeStep  Double array, the forces before ConjugateGradients::makeInfinitesimalStepSaddleSearch was called.
    @param[in]   *forceAfterStep   Double array, the forces returned by ConjugateGradients::makeInfinitesimalStepSaddleSearch.
    @param[in]   maxStep           Double the maximal accepted step. The maximal value of norm of the displacement.*/
    Matrix<double, Eigen::Dynamic, 3> getNewPosModifiedForces(Matrix<double, Eigen::Dynamic, 3> pos, Matrix<double, Eigen::Dynamic, 3> forceBeforeStep, Matrix<double, Eigen::Dynamic, 3> forceAfterStep, double maxStep);
 
    void setForces(Matrix<double, Eigen::Dynamic, 3> forces); // enables the use of modified forces
 
private:
    long nFreeCoord_; // number of free coordinates

    Matter *matter_; // pointer to atom object outside the scope of the class
    Parameters *parameters_; // pointer to a structure outside the scope of the class containing runtime parameters

    /** Double arrrays with size equal to three times the number of atoms. */
    double *tempListDouble_;
    double *direction_;
    double *directionOld_;
    double *directionNorm_;
    double *force_;
    double *forceOld_;

    /** To initialize the object to perform minimization of Matter object.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    void initialize(Matter *matter, Parameters *parameters);

    void determineSearchDirection(); // determine the search direction according to Polak-Ribiere

    /** Determine the size step to be performed.
    @param[in]   *forceBeforeStep  Double array, the forces before ConjugateGradients::makeInfinitesimalStepModifiedForces was called.
    @param[in]   *forceAfterStep   Double array, the forces returned by ConjugateGradients::makeInfinitesimalStepModifiedForces.
    @param[in]   maxStep           Double the maximal accepted step. The maximal value of norm of the displacement.*/
    double stepSize(double *forceBeforeStep, double *forceAfterStep, double maxStep);

    double sign(double value){ return(-(value<0)*2+1);}; // determine the sign of value, \return {double being either 1 or -1}.
};

#endif
