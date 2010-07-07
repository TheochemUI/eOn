/*
 *===============================================
 *  Constants.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/3/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

/** Collection of constants.*/

#define STRING_SIZE 512

namespace constants {
    const std::string VERSION_INFO_EON2("4.4; no Silicon");    

    long getStringSize();///< The size of all char arrays used.  
    
    // Constants used to determine which potential to use.
    long getPotentialNewPotential();///< If a new interface is specified in NewPotential_interface.cpp this constant tells that the new (user provided) potential should be used. Code will crash at runtime if used and no new potential has been defined. (0)
    long getPotentialLJ();///< Constant used to indicate that Lennard-Jones should be used as potential. (1)
    long getPotentialMorse();///< Constant used to indicate that Morse should be used as potential. (2)
    long getPotentialEMT();///< Constant used to indicate that EMT should be used as potential. (3)
    long getPotentialEAM();///< Constant used to indicate that EAM should be used as potential. (11)
    long getPotentialEDIP();///< Constant used to indicate that EDIP should be used as potential. (4)
    long getPotentialVASP();///< Constant used to indicate that VASP should be used as potential. (5)
    long getPotentialTersoff();///< Constant used to indicate that Tersoff should be used as potential. (6)
    long getPotentialSW();///< Constant used to indicate that SW should be used as potential. (7)
    long getPotentialLenosky();///< Constant used to indicate that Lenosky should be used as potential. (8)
    long getPotentialLJBinary();///< Constant used to indicate that LJBinary should be used as potential. (9)
    long getPotentialAluminum();///< Constant used to indicate that Aluminum should be used as potential. (9)
            
    // Values passed between server and clients
    long getStateGood();///< Constant used to indicate the saddle point and prefactor determination was succesful.
    long getStateInit();///< Constant used to indicate the calculation should continue.
    long getStateSaddleSearchNoConvexRegion();///< Constant used to indicate that it was not possible to reach the convex region when doing the saddle point search.
    long getStateSaddleSearchTerminatedBarrier();///< Constant used to indicate that no saddle point was loacted as the search was terminated due to too much energy has been pumped into the system. 
    long getStateSaddleSearchTerminatedConcaveIterations();///< Constant used to indicate that too many concave iterations was made in a series.
    long getStateSaddleSearchTerminatedTotalIterations();///< Constant used to indicate that too many iterations was made.
    long getStateNotConnected();///< Constant used to indicate that saddle point was not connected to the initial state.
    long getStateBadPrefactor();///< Constant used to indicate that either the forward or reverse process have a prefactor that is either bellow getPrefactorMin or above getPrefactorMax.
    long getStateBadBarrier();///< Constant used to indicate that either the forward or reverse process have a barrier that is above the run time parameter in Parameter.
    
    // Constants used in the client
    double getMaxDifferencePos();///< The maximal distance between atoms in two configurations to consider the configurations as equivalent.
    
    // Constants used in prefactor determination
    double getPrefactorMax();///< The maximal prefactor value.
    double getPrefactorMin();///< The minimal prefactor value.
    
    // Constants used in conjugated gradient
    double getCurvatureStep();///< The infinitesimal step size used in ConjugateGradients to determine the optimal step.
    double getMaxMoveFullRelax();///< The maximal norm of the displacement vector at one step during a full minimization.
    double getTimeStep();///< The time step size used in Quickmin.
        
    // Constants used in displacement prior saddle point search
    long getNoPerturbation();///< Constant used to indicate that no perturbation should be performed prior the saddle point search.
    long getPerturbateNotFccOrHcp();///< Constant used to indicate which type of atoms that should be considered in the displacement prior a saddle point search.
    long getPerturbateMinimalCoordinated();///< Constant used to indicate which type of atoms that should be considered in the displacement prior a saddle point search.
    long getPerturbateLastAtom();///< Constant used to indicate which type of atoms that should be considered in the displacement prior a saddle point search.
    double getNeighborCutoff();///< Radius for which atoms that should be used in the analysis determining the local atomic structure.
    
    // Constants used in saddle point determination to set which algorithm that is going to be used for the lowest eigenmode determination.
    long getLowestEigenmodeDimer();///< Constant used to indicate that the Dimer method should be used for the lowest eigenmode determination.
    long getLowestEigenmodeLanczos();///< Constant used to indicate that the Lanczos method should be used for the lowest eigenmode determination.

    // Constants used in saddle point search
//    long getMaximumJumpAttemptsSP();///< The maximal number of jumps allowed to perform before the search for the convex region is terminated.
    long getMaximumIterationsConcaveSeriesSP();///< The maximal number of iterations in a series allowed to perform during the saddle point.    
    //IO file names
    const std::string PARMS_PASSED_FILE_NAME("parameters_passed.dat");///< Name for file containing the run time parameters.
    const std::string REAC_PASSED_FILE_NAME("reactant_passed.con");///< Name of file containing the initial configuration. This structure is relaxed within the code.
    const std::string DISPLACEMENT_PASSED_FILE_NAME("displacement_passed.con");///< Approximate saddle point. The program will refine the saddle point.
    const std::string MODE_PASSED_FILE_NAME("mode_passed.dat");///< Eigenvector. When refining saddle point.    
    const std::string RESULTS_FILE_NAME("results.dat");///< Name for file describing the obtained result.
    const std::string REAC_FILE_NAME("reactant.con");///< Name of file containig the initial configuration. This structure is relaxed within the code.
    const std::string SEND_SP_CONF_FILE_NAME("saddle.con");///< Name of file containing the saddle point configuration if any was found.
    const std::string SEND_PROD_FILE_NAME("product.con");///< Name of file containing the other minima being connected to the saddle point. Note that REAC_FILE_NAME contains the other minima.
    const std::string MODE_FILE_NAME("mode.dat");///< Save last eigenvector.
}
#endif
