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
/* [GH: moved to Potential.h]
    // Constants to identify available potentials
    long potentialNewPotential = 0; // specifies that a user define potential should be defined in NewPotential_interface.cpp
    long potentialLJ       = 1; // Lennard-Jones potential
    long potentialMorse    = 2; // Morse potential
    long potentialEMT      = 3; // Effective Medium Theory potential
    long potentialEDIP     = 4; // EDIP potential
    long potentialVASP     = 5; // VASP DFT
    long potentialTersoff  = 6; // Tersoff Si potential
    long potentialSW       = 7; // Stillinger Webber Si potential
    long potentialLenosky  = 8; // Lenosky Si potential
//  long potentialLJBinary = 9; // Binary Lennard-Jones potential
    long potentialAluminum = 10; // Voter Embedded Atom Al potential and implemented by Hannes Jonsson's group
    long potentialEAM      = 11; // Embedded Atom Al potential
    long potentialQSC      = 12; // Quantum Sutton-Chen potential
    long potentialZpIce    = 13; // Zhu-Phillpott water/Pt potential
    long potentialTip4p    = 14; // TIP4P water potential
*/

/*  [GH: moved to SaddlePoint.h]
    // Return codes passed from server to client to indicate calculation status
    long statusGood = 0; // saddle and prefactor determination was succesful
    long statusInit = 1; // the calculation should continue.
    long statusBadNoConvex = 2; // the convex region was not reached in the saddle search
    long statusBadHighEnergy = 3; // the energy limit was reached and the search terminated (was TerminatedBarrier)
    long statusBadMaxConcaveIterations = 4; // the concave iteration limit was reached (was TerminatedConcaveIterations)
    long statusBadMaxIterations = 5; // the iteration limit was reached (was TerminatedTotalIterations) 
    long statusBadNotConnected = 6; // the saddle point was not connected to the initial state.
    long statusBadPrefactor = 7; // the forward or reverse process have a prefactor that is either below getPrefactorMin or above getPrefactorMax
    long statusBadHighBarrier = 8; // the forward or reverse process have a barrier that is over the energy limit
    long statusBadMinima = 9; // a minimization from the saddle did not converge and neither minima matched the reactant (was MinimaNotConverged)
*/

/*  [GH: moved to Parameters.h for use in Matter.cpp]
    // Constants used in the client
    double maxDifferencePos = 0.1; // The distance criterion for comparing geometries
*/ 

/*  [GH: moved to Parameters.h for use in SaddlePoint.cpp]
    // Constants used in prefactor determination
    double prefactorMax = 10e20; // max prefactor allowed
    double prefactorMin = 10e8; // min prefactor allowed
*/
    
/*  [GH: moved to Parameters.h for use in ConjugateGradients.cpp and Quickmin.cpp]
    // Constants used by the optimizers
    double cgCurvatureStep = 0.001; // finite difference step size used in conjugate gradients
    double cgMaxMoveFullRelax = 0.2; // maximum displacement vector for a step during minimization
    double qmTimeStep = 0.1; // time step size used in Quickmin
*/

/*  [GH: moved to SaddlePoint.h]
    // Constants used to displace atoms before a saddle search
    long dispNone = 0; // make no displacement before a saddle search
    long dispNotFccOrHcp = 1; // displace any atoms which are not HCP or FCC (e.g. in a grain boundary)
    long dispMinCoordinated = 2; // displace the minimum coordinated atoms in the initial configuration
    long dispLastAtom = 3; // displace the last atom in the configuration file
*/

/*  [GH: moved to Parameters.h for use in EpiCenters.cpp]
    double neighborCutoff = 0.33; // radius used in the local atomic structure analysis
*/

/*  [GH: moved to SaddlePoint.h]
    // Constants to define the eigenmode determination algorithm
    long minmodeDimer = 1; // the Dimer method
    long minmodeLanczos = 2; // the Lanczos method
*/

/*  [GH: moved to Parameters.h for use in SaddlePoint.cpp]
    // Constants used in saddle point search
    long maxIterations = 256; // maximum number of iterations in a saddle search
    long maxIterationsConcave = 256; // maximum number of iterations in the concave region
*/

    // File names
    const std::string PARMS_PASSED_FILE_NAME("parameters_passed.dat"); // run time parameters
    const std::string REAC_PASSED_FILE_NAME("reactant_passed.con"); // the reactant configuration which is relaxed within the code
    const std::string DISPLACEMENT_PASSED_FILE_NAME("displacement_passed.con"); // the initial configuration for the saddle search
    const std::string MODE_PASSED_FILE_NAME("mode_passed.dat"); // initial lowest Hessian eigenvector
    const std::string RESULTS_FILE_NAME("results.dat"); // results data
    const std::string REAC_FILE_NAME("reactant.con"); // reactant configuration found
    const std::string SEND_SP_CONF_FILE_NAME("saddle.con"); // saddle point configuration found
    const std::string SEND_PROD_FILE_NAME("product.con"); // the product configuration found
    const std::string MODE_FILE_NAME("mode.dat"); // the last eignvector found
    const std::string READ("r"); // set the file to read
    const std::string APPEND("a"); // set the file to append

}
#endif
