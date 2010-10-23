/*
 *===============================================
 *  EON LowestEigenmodeInterface.h
 *===============================================
 */

#ifndef LOWEST_EIGENMODE_INTERFACE_H
#define LOWEST_EIGENMODE_INTERFACE_H

#include "Matter.h"

#include "Parameters.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN //It hurts every time I type this

/* Define the interface for the lowest eigenvalue determination algorithm */
class LowestEigenmodeInterface{
public:
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
