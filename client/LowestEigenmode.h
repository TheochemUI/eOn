#ifndef LOWESTEIGENMODE_H
#define LOWESTEIGENMODE_H

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"

/* Define the interface for the lowest eigenvalue determination algorithm */
class LowestEigenmode
{

public:
    // stats information
    long totalForceCalls;
    double statsTorque;
    double statsCurvature;
    double statsAngle;
    long statsRotations;
    static const char MINMODE_DIMER[];
    static const char MINMODE_LANCZOS[];

    virtual ~LowestEigenmode() {}

    //void virtual initialize(Matter const *matter, AtomMatrix displacement) = 0;
    virtual void compute(Matter *matter, AtomMatrix direction) = 0;

    virtual double getEigenvalue() = 0;
    virtual AtomMatrix getEigenvector() = 0;
};

#endif
