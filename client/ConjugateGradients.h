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
#include "Optimizer.h"
#include "Parameters.h"
#include "HelperFunctions.h"

/**
 * @file
 * @ingroup Optimizers
 *
 * \brief Direct optimization for energy minimization
 *
 * The conjugate gradient method is an algorithm for the numerical solution of particular
 * stystems of linear equations, namely symetric,positive-definite ones.
 *
 */

/**
 * Decleration of the Conjugate Gradients optimizer
 */

class ConjugateGradients : public Optimizer
{
    public:
	//! Conjugate Gradients optimizer constructor
	/*!
 	 * \param *objf an ref ObjectiveFunction that tells the optimizer how to run
 	 * \param *parameters defined by the config.init file
 	 */	 
        ConjugateGradients(ObjectiveFunction *objf, Parameters *parameters);
	//! Conjugant Gradient deconstructor
        ~ConjugateGradients();

	//! Calls the next step in the algorithm
	/**
 	 * Either calls the single_step or line_search method depending on the parameters
 	 * \return whether or not the algorithm has converged
 	 */
        bool step(double maxMove);
	//! Runs the conjugate gradient 
	/**
 	 * \todo method should also return an error code and message if the algorithm errors out
 	 * \return algorithm convergence
	 */
        bool run(int maxIterations, double maxMove);
	//! Gets the direction of the next step
        VectorXd getStep();

    private:
	//! An objective function relating a certain job method to the conjugate gradient optimizer
        ObjectiveFunction *objf;
	//! Parameters set by the config.init file
        Parameters *parameters;

	//! Current step direction of the conjugate gradient
        VectorXd direction;
	//! Algorithms previous step direction
        VectorXd directionOld;
	//! Normalised version of the current direction vector
        VectorXd directionNorm;
	//! Current force vector
        VectorXd force;
	//! Previous force vector
        VectorXd forceOld;
    
        //! Counts the number of descrete steps untill algorithm convergence
	int cg_i;
    
	//! Steps the conjugate gradient
	/**
 	 *  Checks for convergence based on the change in displacement or direction	
         */
	bool single_step(double maxMove);
	//! Steps the conjugate gradient
	/**
 	* Checks for convergence based on the ratio of the projected force along a line
 	* to the norm of the total force 
 	*/ 
        bool line_search(double maxMove);
    
};

#endif
