#ifndef JOB_H
#define JOB_H
#include "Parameters.h"
#include <vector>
#include <string>

/** @defgroup Jobs
*
* \brief ClientEON main procedures
*
* This page provides links to all of the available jobs that can be run by the 
* ClientEON, as well as documentation on the job class, and the overview section relating
* the job structure to the rest of the program. 
* 
*/

/**
 * @file
 * @ingroup Jobs
 *
 * \brief The job class is used to serve as an abstract class for all jobs,
 *  as well as to call a job at runtime based off of the passed in parameters. 
 *
 * The Static members are used to tell at runtime which job to run as set by the parameters,
 * and therefore as set by the config.init file. About half of the jobs are standalone, while
 * others are run from routines with the same name. A certain subset of jobs do not run
 * optimizers (SEE OVERVIEW) and are documented in their own files accordingly.
 *
 */

/**
 * Decleration of job class
 */

class Job { 
    public:
	
	//! job deconstructor.
	virtual ~Job() {};
	//! Virtual run; used solely for dynamic dispatch 
	virtual std::vector<std::string> run()=0;
	//! Runs a \ref ProcessSearchJob "Process Search" job 
        static const char PROCESS_SEARCH[];
	//! Runs a \ref SaddleSearchJob "Saddle Search" job
        static const char SADDLE_SEARCH[];
	//! ref MinimizationJob
        static const char MINIMIZATION[];
        //! ref PointJob
	static const char POINT[];
        //! ref ParallelReplicaJob
	static const char PARALLEL_REPLICA[];
        //! ref SafeHyperdynamicsJob
	static const char SAFE_HYPER[];
	//! ref TADJob
        static const char TAD[];
	//! ref ReplicaExchangeJob
        static const char REPLICA_EXCHANGE[];
	//! ref BasinHopingJob
        static const char BASIN_HOPPING[];
	//! ref HessianJob
        static const char HESSIAN[];
	//! ref FiniteDifferenceJob
        static const char FINITE_DIFFERENCE[];
	//! ref NudgedElasticBandJob
        static const char NUDGED_ELASTIC_BAND[];
	//! ref DynamicsJob
        static const char DYNAMICS[];
	//! ref PrefactorJob
        static const char PREFACTOR[];
	//! ref GlobalOptimizationJob
	static const char GLOBAL_OPTIMIZATION[];
	//! ref StructureComparisonJob
	static const char STRUCTURE_COMPARISON[];
	//! ref MonteCarloJob
        static const char MONTE_CARLO[];
	//! ref TestJob
        static const char TEST[];
	//! Grabs the correct job as specified by the parameters
	/*!
 	  \param *parameters defined by the config.init file
	*/
        static Job *getJob(Parameters *parameters);
};
#endif
