//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef CG_H
#define CG_H

#include "Eigen.h"
#include "Matter.h"
#include "Minimizer.h"
#include "Parameters.h"
#include "HelperFunctions.h"

/* Functionality relying on the conjugate gradients algorithm.
 * The object is capable of minimizing a Matter object or modified forces being passed in.*/
class ConjugateGradients : public Minimizer{

public:
    /** Constructor to be used when a structure is minimized.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    ConjugateGradients(Matter *matter, Parameters *parameters);

    /** Constructor to be used when modified forces are used.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.
    @param[in]   *forces        Double array containing the forces acting on the matter.*/
    ConjugateGradients(Matter *matter, Parameters *parameters, AtomMatrix forces);

    ~ConjugateGradients();///< Destructor.

    void oneStep(); // do one iteration
    long fullRelax(); // relax the Matter object at the pointer that was passed with the constructor
    bool isItConverged(double convergeCriterion);
    void setOutput(int level);
    //----Functions used when forces are modified (saddle point search)----
    /** Performs and infinitesimal step along the search direction.
    @param[out]  *posStep  Double array containing the modified positions of the atoms.
    @param[in]   *pos      Double array containing the original positions of the atoms.*/
    AtomMatrix makeInfinitesimalStepModifiedForces(AtomMatrix pos);

    /** Performs and infinitesimal step along the search direction.
    @param[out]  *pos              Double array containing the calculated new positions of the atoms.
    @param[in]   *forceBeforeStep  Double array, the forces before ConjugateGradients::makeInfinitesimalStepSaddleSearch was called.
    @param[in]   *forceAfterStep   Double array, the forces returned by ConjugateGradients::makeInfinitesimalStepSaddleSearch.
    @param[in]   maxStep           Double the maximal accepted step. The maximal value of norm of the displacement.*/
    AtomMatrix getNewPosModifiedForces(AtomMatrix pos, AtomMatrix forceBeforeStep, AtomMatrix forceAfterStep, double maxStep);
 
    void setForces(AtomMatrix forces); // enables the use of modified forces
 
private:
    long nAtoms; // number of free coordinates
    int outputLevel;

    Matter *matter; // pointer to atom object outside the scope of the class
    Parameters *parameters; // pointer to a structure outside the scope of the class containing runtime parameters

    /** Double arrrays with size equal to three times the number of atoms. */
    AtomMatrix direction;
    AtomMatrix directionOld;
    AtomMatrix directionNorm;
    AtomMatrix force;
    AtomMatrix forceOld;

    /** To initialize the object to perform minimization of Matter object.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    void initialize(Matter *matter, Parameters *parameters);

    void determineSearchDirection(); // determine the search direction according to Polak-Ribiere

    /** Determine the size step to be performed.
    @param[in]   *forceBeforeStep  Double array, the forces before ConjugateGradients::makeInfinitesimalStepModifiedForces was called.
    @param[in]   *forceAfterStep   Double array, the forces returned by ConjugateGradients::makeInfinitesimalStepModifiedForces.
    @param[in]   maxStep           Double the maximal accepted step. The maximal value of norm of the displacement.*/
    AtomMatrix getStep(AtomMatrix forceBeforeStep, AtomMatrix forceAfterStep, double maxStep);

    double sign(double value){ return(-(value<0)*2+1);}; // determine the sign of value, \return {double being either 1 or -1}.
};

#endif
