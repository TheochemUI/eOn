/*
 *===============================================
 *  LowestEigenmodeInterface.h
 *  eon2
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 1/3/07.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef LOWEST_EIGENMODE_INTERFACE_H
#define LOWEST_EIGENMODE_INTERFACE_H

#include "Matter.h"

#include "Parameters.h"

/** Defining the interface for the lowest eigenvalue determination algorithm.*/
class LowestEigenmodeInterface{
public:
    virtual ~LowestEigenmodeInterface(){};
    void virtual startNewSearchAndCompute(Matter const *matter, double *displacement) = 0; 
    void virtual moveAndCompute(Matter const *matter) = 0;  
    double virtual returnLowestEigenmode(double *result) = 0;
    /// Return eigenvector.
    virtual double const * getEigenvector(long & size) const
    {size=0; return 0;}
        /** Set initial direction manually.*/
    virtual void setEigenvector(long size, double const eigenvector[]) {}
};
#endif
