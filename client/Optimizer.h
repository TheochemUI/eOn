//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Parameters.h"
#include "ObjectiveFunction.h"
#include "Eigen.h"

/** @defgroup Optimizers
 *
 * \brief ClientEON methods for optimizing atomic structures
 *
 * This page provides links to all the available optimizers that can be run by the
 * ClientEON, as well as documentation on the optimizer class. 
 *
 */ 

/**
 * @file
 * @ingroup Optimizers
 *
 * \brief The optimizer class is used to serve as an abstract class for all optimizers,
 *  as well as to call an optimizer at runtime based off of the passed in paramters.
 *
 * The set of optimizers are methods for optimizing atomic structures, solving unconstrained 
 * energy minimization. Only a certain set of job types that ClientEON runs can take advantage
 * of the numeric optimizer (SEE OVERVIEW) and are documented in their own files accordingly. 
 * 
 */

/**
 * Decleration of the optimizer class
 */

class Optimizer
{
    public:
	//! optimizer deconstructor
        virtual ~Optimizer(){};
	//! Template for stepping the optimizer, returns convergence
        virtual bool step(double maxMove) = 0;
	//! Template for runnning the optimizer; uses a series of steps, checking for convergence each time
        virtual bool run(int maxIterations, double maxMove) = 0;
        //! Grabs the correct optimizer as specified by the parameters
        /*!
 	 * \param *objf an ref ObjectiveFunction that tells the optimizer how to run
 	 * \param *parameters defined by the config.init file
 	 */	 
	static Optimizer *getOptimizer(ObjectiveFunction *objf, Parameters *parameters);
};

#endif
