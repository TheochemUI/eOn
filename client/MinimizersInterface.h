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
 *  EON MinimizersInterface.h
 *===============================================
 */
#ifndef MINIMIZERS_INTERFACE_H
#define MINIMIZERS_INTERFACE_H

/* Define the interface for the potential */
class MinimizersInterface{
    
public:
    virtual ~MinimizersInterface(){};
    void virtual oneStep() = 0; // Do one iteration
    void virtual fullRelax() = 0; // Relax the Matter object corresponding to the pointer that was passed.
    bool virtual isItConverged(double convergeCriterion) = 0; // Determine if the norm of the force vector is below the convergeCriterion
    long totalForceCalls;
};
#endif
