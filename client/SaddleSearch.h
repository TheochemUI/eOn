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
#include "LowestEigenmodeInterface.h" 

#include <string>

#include "Eigen.h"

using namespace std;

class Matter;
class Parameters;
class LowestEigenmodeInterface;

/** Find a saddle point using a minimum mode following method.
    Requires a function to determine the epicenter of atomic displacements from a minimum.
    Requires a method to find the minimum mode */
class SaddlePoint
{

public:

    // Return codes passed from server to client to indicate calculation status
    enum{
        STATUS_GOOD,
        STATUS_INIT,
        STATUS_BAD_NO_CONVEX,
        STATUS_BAD_HIGH_ENERGY,
        STATUS_BAD_MAX_CONCAVE_ITERATIONS,
        STATUS_BAD_MAX_ITERATIONS,
        STATUS_BAD_NOT_CONNECTED,
        STATUS_BAD_PREFACTOR,
        STATUS_BAD_HIGH_BARRIER,
        STATUS_BAD_MINIMA
    };

    // Constants used to displace atoms before a saddle search
    enum{
        DISP_LOAD,
        DISP_NOT_FCC_OR_HCP,
        DISP_MIN_COORDINATED,
        DISP_LAST_ATOM,
        DISP_RANDOM
    };

    // Methods for finding the minimum mode
    enum{
        MINMODE_DIMER,
        MINMODE_LANCZOS,
        #warning deprecated
        MINMODE_OLDLANCZOS, // remove when tests finished
        MINMODE_EXACT
    };
 
    SaddlePoint(); // The object shall be initialized later with SaddlePoint::initialize
 
    /** Constructor
    @param[in]  initial      initial state minimum
    @param[in]  saddle       conformation where to start the saddle point search; also used to return the final saddle point.
    @param[in]  *parameters  runtime parameters */
    SaddlePoint(Matter *initial, Matter *saddle, Parameters *parameters);

    ~SaddlePoint(); // destructor

    void initialize(Matter *initial, Matter *saddle, Parameters *parameters);

    void displaceAndSetMode(Matter *matter); // make a displacement of atoms centered on the EpiCenter atom and set the initial mode accordingly  

    /** Determine a nearby saddle point
    @param[out]  *min1     one of the minima connected to the saddle point
    @param[out]  *min2     the other minima connected to the saddle point
    The value returned is true if the calculation converges. */
    long locate();//(Matter *min1, Matter *min2);
    LowestEigenmodeInterface const * getLowestEigenmode() const;
    long getnFreeCoord() const;
	
	Matrix<double, Eigen::Dynamic, 3> getSaddlePositions();
    Matrix<double, Eigen::Dynamic, 3> getEigenMode();

    Matrix<double, Eigen::Dynamic, 3> mode;
    void loadMode(string filename);
    void loadMode(FILE * modeFile);
    void saveMode(FILE * modeFile);

    long forceCallsSaddlePointConcave;
    long forceCallsSaddlePointConvex;

private:
    Matter *initial;
    Matter *saddle;
    Parameters *parameters;
    LowestEigenmodeInterface *lowestEigenmode; // method used to determine the lowest eigenmode
 
    double eigenValue; // estimate for the lowest eigenvalue
    Matrix<double, Eigen::Dynamic, 3> eigenMode; // lowest eigenmode
    Matrix<double, Eigen::Dynamic, 3> initialDisplacement;
    long nFreeCoord; // number of free coordinates
    long status;

    void clean(); // clean up dynamical allocated memory

    Matrix<double, Eigen::Dynamic, 3> projectedForce(Matrix<double, Eigen::Dynamic, 3> force); // projected minmode force

    /** Determine the two minima connected to the saddle point by displacing forward and backward along the lowest eigenmode from the saddle and minimizing
    @param[out]  *min1   one minima connected to the saddle point
    @param[out]  *min2   the other minima connected to the saddle point */

    void displaceInConcaveRegion();

    void searchForSaddlePoint(double initialEnergy);
    void addForceCallsSaddlePoint(long fcalls, double eigenvalue);

};

#endif
