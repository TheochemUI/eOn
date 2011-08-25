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
class SaddleSearch
{

public:

    // Return codes passed from server to client to indicate calculation status
    // XXX: Shouldn't this be in ProcessSearchJob.cpp?
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
        STATUS_CHECKPOINT, //12
    };

    // Constants used to displace atoms before a saddle search
    static const char DISP_LOAD[];
    static const char DISP_NOT_FCC_OR_HCP[];
    static const char DISP_MIN_COORDINATED[];
    static const char DISP_LAST_ATOM[];
    static const char DISP_RANDOM[];

    // Methods for finding the minimum mode
    static const char MINMODE_DIMER[];
    static const char MINMODE_LANCZOS[];
    static const char MINMODE_EXACT[];
 
    SaddleSearch(); // The object shall be initialized later with SaddleSearch::initialize
 
    /** Constructor
    @param[in]  initial      initial state minimum
    @param[in]  saddle       conformation where to start the saddle point search; also used to return the final saddle point
    @param[in]  *parameters  runtime parameters */
    SaddleSearch(Matter *initial, Matter *saddle, Parameters *parameters);

    ~SaddleSearch(); // destructor

    void initialize(Matter *initial, Matter *saddle, Parameters *parameters);

    void displaceAndSetMode(Matter *matter); // make a displacement of atoms centered on the EpiCenter atom and set the initial mode accordingly

    /** Determine a nearby saddle point
    @param[out]  *min1     one of the minima connected to the saddle point
    @param[out]  *min2     the other minima connected to the saddle point
    The value returned is true if the calculation converges. */
    long locate();//(Matter *min1, Matter *min2);
    LowestEigenmodeInterface const * getLowestEigenmode() const;

    AtomMatrix getSaddlePositions();
    AtomMatrix getEigenMode();
    double getEigenValue();

    AtomMatrix mode;
    void loadMode(string filename);
    void loadMode(FILE * modeFile);
    void saveMode(FILE * modeFile);

    long forceCallsSaddleSearchConcave;
    long forceCallsSaddleSearchConvex;
    long iterations;

private:
    Matter *initial;
    Matter *saddle;
    Parameters *parameters;
    LowestEigenmodeInterface *lowestEigenmode; // method used to determine the lowest eigenmode

    double eigenValue; // estimate for the lowest eigenvalue
    AtomMatrix eigenMode; // lowest eigenmode
    AtomMatrix initialDisplacement;
    long status;

    void clean(); // clean up dynamically allocated memory

    AtomMatrix projectedForce(AtomMatrix force); // projected minmode force

    /** Determine the two minima connected to the saddle point by displacing forward and backward along the lowest eigenmode from the saddle and minimizing
    @param[out]  *min1   one minima connected to the saddle point
    @param[out]  *min2   the other minima connected to the saddle point */

    void displaceInConcaveRegion();

    void searchForSaddlePoint(double initialEnergy);
    void addForceCallsSaddleSearch(long fcalls, double eigenvalue);

};

#endif
