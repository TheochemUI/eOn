/*
 *===============================================
 *  Mdynamics.h
 *-----------------------------------------------
 *  Created by Liang Zhang on 4/17/07.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef MDYNAMICS_H
#define MDYNAMICS_H

#include "MinimizersInterface.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Constants.h"
#include "Parameters.h"

/** Functionality relying on the conjugate gradients algorithm. The object is capable of minimizing an Matter object or modified forces being passed in.*/
class Mdynamics {

public:
    /** Constructor to be used when a structure is minimized.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    Mdynamics(Matter *matter, Parameters *parameters);

    ~Mdynamics();///< Destructor.

 	void oneStep();	
	void VerletStep1();
	void VerletStep2();
	void fullSteps();
	void Andersen();
private:
    long nFreeCoord_;///< Number of free coordinates.

    Matter *matter_;///< Pointer to atom object \b outside the scope of the class.    
    Parameters *parameters_;///< Pointer to a structure outside the scope of the class containing runtime parameters. 

    double *tempListDouble_;///< Double array, its size equals the number of atoms times 3.
    double dtScale_;
};

#endif
