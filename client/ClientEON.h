/*
 *===============================================
 *  ClientEON.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/24/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef CLIENT_EON_H
#define CLIENT_EON_H

// standard c libraries
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cassert>

#ifdef BOINC
	#include <boinc/boinc_api.h>
	#include <boinc/diagnostics.h>     // boinc_init_diagnostics()
	#include <boinc/filesys.h>         // boinc_fopen(), etc...
#else
	#include "false_boinc.h"
#endif

// includes for boinc
#ifdef WIN32
	#include <boinc_win.h>
	#include <win_util.h>
#endif

#include "Constants.h"
#include "Parameters.h"
#include "Matter.h"
#include "ConjugateGradients.h"
#include "QMBox.h"
#include "Prefactors.h"
#include "SaddlePoint.h"
#include "LowestEigenmodeInterface.h"  

/** The main program.*/
namespace client_eon {

      /// Load the unperturbed configuration and the run time parameters received by BOINC from the master.
      void doSaddleSearch(void);
      void loadDataAndRelax(char const parameters_passed[], char const reactant_passed[]);
      bool loadDisplacementAndMode(char const displacement_passed[], char const mode_passed[]);
      void initializeNewSearch();///< Initialize all variables used in the saddle point search.
      bool connectedToReactantState();///< Return true if one of the two states separated by the saddle point is equivalent to the conf received from the master 
      void determineBarriers();///< Determine and store the barrier heights for both the forward and reverse processes.
      bool barriersWithinWindow();///< Return true if the barriers are within the window defined in the init.dat file.
      bool prefactorsWithinWindow();///< Return true if the prefactors are within the window defined.
      void saveData(); ///< Save result data.

      void printEndState(long state);///< Print out the end state of the saddle point determination.
      void printRequestedInfo(char *argv); ///< Print help as requested by user.
      void forcesOfConfig();  ///prints forces and energy to potentialInfo.txt
     
      int divineBundleSize(char const mode_passed[]); // Attempts to find the bundle size

      int rc;///< BOINC related, return code from various functions.
      int state;///< Variable is equal to GOOD if a converged saddle point is obtained.

      Matter *initial;///< The initial unperturbed configuration. 
      Matter *saddle;///< Configuration used during the saddle point search.
      Matter *displacement;///< Configuration used during the saddle point search.
      Matter *min1;///< Used to determine one of the stable states connected to the saddle point.
      Matter *min2;///< Used to determine the other stable state connected to the saddle point.
      Matter *testConfig;  /// test configuration from reactant.test

      SaddlePoint saddlePoint;///< Used to determine the saddle point.
      double barriersValues[2];///< First element is from min1 to saddle, second is from min2 to saddle.
      bool barriersOK;///< Is set if the barriers fall into the window defined in Constants.h.

      Prefactors prefactors;///< Used to determine the other determine the prefactors.
      double prefactorsValues[2];///< First element is from min1 to saddle, second is from min2 to saddle.
      bool prefactorsOK;///< Is set if the prefactors fall into the window defined in Constants.h.
      Parameters parameters;///< Store all the runtime parameters received from the server.

};
#endif
