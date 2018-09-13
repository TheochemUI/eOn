#ifndef MIN_MODE_SADDLE_SEARCH_H
#define MIN_MODE_SADDLE_SEARCH_H

#include "SaddleSearchMethod.h"
#include "Matter.h"
#include "LowestEigenmode.h"
#include "Eigen.h"
#include "Optimizer.h"

#include <string>

class MinModeSaddleSearch : public SaddleSearchMethod
{

public:

    enum{
        // DO NOT CHANGE THE ORDER OF THIS LIST
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
        STATUS_NONLOCAL_ABORT, //13
        STATUS_NEGATIVE_BARRIER, //14
        STATUS_BAD_MD_TRAJECTORY_TOO_SHORT, //15
        STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE, //16
        STATUS_BAD_NO_BARRIER, //17
        STATUS_ZEROMODE_ABORT, //18
        STATUS_OPTIMIZER_ERROR //19
    };

    MinModeSaddleSearch(Matter *matterPassed, AtomMatrix modePassed,
                 double reactantEnergy, Parameters *parametersPassed);
    ~MinModeSaddleSearch();
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
