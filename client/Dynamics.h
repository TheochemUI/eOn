//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "Optimizer.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

#include "Eigen.h"

class Dynamics {

public:

    static const char ANDERSEN[];
    static const char NOSE_HOOVER[];
    static const char LANGEVIN[];
    static const char NONE[];

    Dynamics(Matter *matter, Parameters *parameters);

    ~Dynamics();

    void oneStep(double temperature);
    void andersenVerlet();
    void velocityVerlet();
    void fullSteps(double temperature); 
    void andersen(double temperature);
    void initialVel(double temperature);
    void velRescaling(double temperature);
    void noseHooverVerlet(double temperature);
    void langevinVerlet(double temperature);
    long getMDfcalls();
    long getMinfcalls();
    long getRefinefcalls();
    bool checkState(Matter *matter,Matter *min1);
    long refine(Matter *buff[],long length,Matter *min1);

//  AtomMatrix langevinAcc(double temperature);

private:
    long nAtoms;

    Matter *matter;
    Parameters *parameters;

    long min_fcalls;
    long md_fcalls;
    long rf_fcalls;
    double dt;
    double kb;
    bool init;
    double vxi1,vxi2,xi1,xi2;
};

#endif
