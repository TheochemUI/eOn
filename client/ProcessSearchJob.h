
#ifndef PROCESSSEARCHJOB_H
#define PROCESSSEARCHJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"
#include "MinModeSaddleSearch.h"

/**
 * @file
 * @ingroup Jobs
 *  
 * \brief Finds possible escape mecahnisms from a state. 
 *
 * The process search job implements one of the following types of saddle searches, as 
 * well as an \ref Optimizer "optimizer", as defined by the config.init file.
 *
 * <ul>
 * <li> ref MinModeSaddleSearch
 * <li> ref BasinHopingSaddleSearch
 * <li> ref DynamicsSaddleSearch
 * <li> ref BiasedGradientSquaredDescent
 * </ul>
 *
 */

/**
 * Decleration of the Process Search job
 */

class ProcessSearchJob : public Job {
    public:
	//! Process Search job Constructor
	/*!
 	 * \param *params defined by the config.init file
 	 */	
        ProcessSearchJob(Parameters *params);
	//! Process Search job Deconstructor
        ~ProcessSearchJob(void);
	//! Kicks off the Process Search	
	std::vector<std::string> run(void);

    private:
	//! Runs the correct saddle search; also checks if the run was successful
        int  doProcessSearch(void);
 	//! UNDEFINED 
 	/*! 
 	 * Function is not defined in the cpp file and is marked for deletion.
 	 */
        VectorXi movedAtoms(void);
	//! Logs the run status and makes sure the run was successful
        void printEndState(int status);
	//! Writes the results from the run to file
        void saveData(int status);

        //! Container for the results of the run
	std::vector<std::string> returnFiles;

	//! Parameters defined by the config.init file
        Parameters *parameters;
	//! Pulled from parameters
	/*!
 	 * sa ref BasinHoppingSaddleSearch DynamicsSaddleSearch BiasedGradientSquaredDescent MinModeSaddleSearch
 	 */ 	
        SaddleSearchMethod *saddleSearch; 
        //! Initial configuration
	Matter *initial;
	//! Configuration used during the saddle point search      
        Matter *saddle;
        //! Configuration used during the saddle point search
	Matter *displacement; 
        //! First minimum from the saddle
	Matter *min1;  
	//! Second minimum from the saddle       
        Matter *min2;   
        
	//! Array containg the value of the potential barriers from reactant to product and vice versa	
        double barriersValues[2];
	//! Array containg the value of the prefactors of the matter from reactant to product and vice versa 
        double prefactorsValues[2];

	//! Force calls to find the saddle
        int fCallsSaddle;
	//! Force calls to minimize
        long fCallsMin;
	//! Force calls to find the prefactors
        int fCallsPrefactors;
};

#endif
