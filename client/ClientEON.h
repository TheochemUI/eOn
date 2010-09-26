/*
 *===============================================
 *  EON ClientEON.h
 *===============================================
 */

#ifndef CLIENT_EON_H
#define CLIENT_EON_H

#include "Matter.h"
#include "Prefactors.h"
#include "SaddlePoint.h"

/** The main program.*/
namespace client_eon {

    void doSaddleSearch(void);
    void loadDataAndRelax(char const parameters_passed[], char const reactant_passed[]);
    bool loadDisplacementAndMode(char const displacement_passed[], char const mode_passed[]);
    void initializeNewSearch(); // initialize all variables used in the saddle point search
    bool connectedToReactantState(); // true if a minimum connect to the reactant passed from the master
    void determineBarriers(); // determine and store the barrier heights for both the forward and reverse processes
    bool barriersWithinWindow(); // true if the barriers are within the window defined in the init.dat file
    bool prefactorsWithinWindow(); // true if the prefactors are within the window defined
    void saveData(); // save result data

    void printEndState(long status); // print the status of the saddle point calculation
    void printRequestedInfo(char *argv); // print help as requested by user
    void forcesOfConfig();  // print forces and energy to potentialInfo.txt
 
    int divineBundleSize(char const mode_passed[]); // Attempts to find the bundle size

    int rc; // BOINC related, return code from various functions
    int status; // return status of the saddle point calculation

    Matter *initial; // initial configuration.
    Matter *saddle; // configuration used during the saddle point search.
    Matter *displacement; // configuration used during the saddle point search.
    Matter *min1; // first minimum from the saddle
    Matter *min2; // second minimum from the saddle
    Matter *testConfig;  // test configuration from reactant.test

    SaddlePoint saddlePoint; // used to determine the saddle point
    double barriersValues[2]; // first element is from min1 to saddle, second is from min2 to saddle
    bool barriersOK; // set if the barriers fall into the window defined in Constants.h

    Prefactors prefactors; // used to determine the other determine the prefactors.
    double prefactorsValues[2]; // first element is from min1 to saddle, second is from min2 to saddle
    bool prefactorsOK; // is set if the prefactors fall into the window defined in Constants.h
    Parameters parameters; // store all the runtime parameters received from the server
    int result_pattern(char *filename);

}
#endif
