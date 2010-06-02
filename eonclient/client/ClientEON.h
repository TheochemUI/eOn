/*
 *===============================================
 *  ClientEON.h
 *  eon2
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
#include <iostream>
#include <stdlib.h>
#include <time.h>

// includes for boinc
#ifdef WIN32
#include <boinc_win.h>
#include <win_util.h>
#endif

#include "common/Constants.h"
#include "common/Parameters.h"
#include "common/Matter.h"
#include "common/ConjugateGradients.h"
#include "Prefactors.h"
#include "SaddlePoint.h"
/** The main program.*/
namespace client_eon {
      struct Vector {
            ~Vector() {delete [] vector;}
            long size;
            double * vector;
      } eigenvector;
      /// Load the unpertubed configuration and the run time parameters recieved by BOINC from the master.
      void doSaddleSearch(void);
      void loadDataAndRelax(char const parameters_passed[], char const reactant_passed[]);
      bool loadDisplacementAndMode(char const displacement_passed[], char const mode_passed[]);
      void initializeNewSearch();///< Initialize all variables used in the saddle point search.
      bool connectedToReactantState();///< Return true if one of the two states seperated by the saddle point is equivalent to the conf. recieved from the master 
      void determineBarriers();///< Determine and store the barrier heights for both the forward and reverse process.
      bool barriersWithinWindow();///< Return true if the barriers lay within the window defined in the init.dat file.
      bool prefactorsWithinWindow();///< Return true if the prefactors lay within the window defined.
      void saveData(); ///< Save result data.

      void printEndState(long state);///< Print out the end state of the saddle point determination.
      void printRequestedInfo(char *argv); ///< Print help as requested by user.


      int rc;///< BOINC related, return code from various functions.
      int state;///< Variable is equal to GOOD if a converged saddle point is obtained.

      Matter *initial;///< The initial unpertubed configuration. 
      Matter *saddle;///< Configuration used during the saddle point search.
      Matter *min1;///< Used to determine one of the stable states connected to the saddle point.
      Matter *min2;///< Used to determine the other stable state connected to the saddle point.

      SaddlePoint saddlePoint;///< Used to determine the other determine saddle point.
      double barriersValues[2];///< First element is from min1 to saddle, second is from min2 to saddle.
      bool barriersOK;///< Is set if the barriers fall into the window defined in Constants.h.

      Prefactors prefactors;///< Used to determine the other determine the prefactors.
      double prefactorsValues[2];///< First element is from min1 to saddle, second is from min2 to saddle.
      bool prefactorsOK;///< Is set if the prefactors fall into the window defined in Constants.h.
      Parameters parameters;///< Store all the runtime parameters recieved from the server
};
#endif
