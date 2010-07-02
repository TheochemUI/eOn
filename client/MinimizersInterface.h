/*
 *===============================================
 *  MinimizersInterface.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 4/17/07.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef MINIMIZERS_INTERFACE_H
#define MINIMIZERS_INTERFACE_H

/** Defining the interface for the potential.*/
class MinimizersInterface{
    
public:
    virtual ~MinimizersInterface(){};
    void virtual oneStep() = 0;///< Do one iteration.
    void virtual fullRelax() = 0;///< Relax the Matter object corresponding to the pointer that was passed.
    bool virtual isItConverged(double convergeCriterion) = 0;///< Determine if the norm of the force vector is bellow the convergeCriterion.
};
#endif
