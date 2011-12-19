//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef SADDLESEARCH_H
#define SADDLESEARCH_H

#include "Matter.h"
#include "LowestEigenmode.h"
#include "Eigen.h"
#include "Optimizer.h"

#include <string>

class SaddleSearch
{

public:

    enum{
        // DONT CHANGE THE ORDER OF THIS LIST
        STATUS_GOOD, //0
        STATUS_INIT, //1
        STATUS_BAD_NO_CONVEX, //2
        STATUS_BAD_HIGH_ENERGY, //3
        STATUS_BAD_MAX_CONCAVE_ITERATIONS, //4
        STATUS_BAD_MAX_ITERATIONS, //5
        STATUS_BAD_NOT_CONNECTED, //6
        STATUS_BAD_PREFACTOR, //7
        STATUS_BAD_HIGH_BARRIER, //8
        STATUS_BAD_MINIMA, //9
        STATUS_FAILED_PREFACTOR, //10
        STATUS_POTENTIAL_FAILED, //11
        STATUS_NONNEGATIVE_ABORT, //12
    };

    // Methods for finding the minimum mode
    static const char MINMODE_DIMER[];
    static const char MINMODE_LANCZOS[];
    static const char MINMODE_EXACT[];
 
    SaddleSearch(Matter *matterPassed, AtomMatrix modePassed, double reactantEnergy, Parameters *parametersPassed);
    ~SaddleSearch();
    AtomMatrix getEigenvector();
    double getEigenvalue();

    int run();

    int iteration;
//    double eigenValue; // estimate for the lowest eigenvalue
//    AtomMatrix eigenMode; // lowest eigenmode
    int status;

private:
    AtomMatrix mode;
    Matter *matter;
    Parameters *parameters;
    LowestEigenmode *minModeMethod;
    double reactantEnergy;
};

#endif
