//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
/*
 *===============================================
 *  EON Quickmin.h
 *===============================================
 */
#ifndef QUICKMIN_H
#define QUICKMIN_H

#include "MinimizersInterface.h"
#include "Matter.h"
//#include "HelperFunctions.h"
//#include "Constants.h"
//#include "Parameters.h"

/** Functionality relying on the conjugate gradients algorithm. The object is capable of minimizing an Matter object or modified forces being passed in.*/
class Quickmin : public MinimizersInterface{

public:
    /** Constructor to be used when a structure is minimized.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    Quickmin(Matter *matter, Parameters *parameters);

    ~Quickmin();///< Destructor.

    void oneStep();///< Do one iteration.
    void oneStepPart1(Matrix<double, Eigen::Dynamic, 3> forces);
    void oneStepPart2(Matrix<double, Eigen::Dynamic, 3> forces);
    
    void fullRelax();///< Relax the Matter object corresponding to the pointer that was passed with the constructor.
    bool isItConverged(double convergeCriterion);///< Determine if the norm of the force vector is bellow the \a convergeCriterion.
    
private:
    long nAtoms;///< Number of free coordinates.

    Matter *matter;///< Pointer to atom object \b outside the scope of the class.    
    Parameters *parameters;///< Pointer to a structure outside the scope of the class containing runtime parameters. 

    Matrix<double, Eigen::Dynamic, 3> forces;///< Double array, its size equals the number of \b free atoms times 3.
    double dtScale;

    Matrix<double, Eigen::Dynamic, 3> velocity; 
};

#endif
