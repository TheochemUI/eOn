//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LOWEST_EIGENMODE_INTERFACE_H
#define LOWEST_EIGENMODE_INTERFACE_H

#include "Matter.h"

#include "Parameters.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN //It hurts every time I type this

/* Define the interface for the lowest eigenvalue determination algorithm */
class LowestEigenmodeInterface{
public:
    virtual ~LowestEigenmodeInterface() {}
    void virtual startNewSearchAndCompute(Matter const *matter, Matrix<double, Eigen::Dynamic, 3>) = 0; 
    void virtual moveAndCompute(Matter const *matter) = 0;  
    double virtual getEigenvalue() = 0;
    double *stats;
    
    /// Return eigenvector.
    virtual Matrix<double, Eigen::Dynamic, 3> getEigenvector()=0;
        /** Set initial direction manually.*/
    virtual void setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector)=0;
};
#endif
