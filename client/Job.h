#ifndef JOB_H
#define JOB_H

#include "Parameters.h"
#include "HelperFunctions.h"

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
        static const std::string PROCESS_SEARCH;
	//! Runs a \ref SaddleSearchJob "Saddle Search" job
        static const std::string SADDLE_SEARCH;
	//! ref MinimizationJob
        static const std::string MINIMIZATION;
        //! ref PointJob
	static const std::string POINT;
        //! ref ParallelReplicaJob
	static const std::string PARALLEL_REPLICA;
        //! ref SafeHyperdynamicsJob
	static const std::string SAFE_HYPER;
	//! ref TADJob
        static const std::string TAD;
	//! ref ReplicaExchangeJob
        static const std::string REPLICA_EXCHANGE;
	//! ref BasinHopingJob
        static const std::string BASIN_HOPPING;
	//! ref HessianJob
        static const std::string HESSIAN;
	//! ref FiniteDifferenceJob
        static const std::string FINITE_DIFFERENCE;
	//! ref NudgedElasticBandJob
        static const std::string NUDGED_ELASTIC_BAND;
	//! ref DynamicsJob
        static const std::string DYNAMICS;
	//! ref PrefactorJob
        static const std::string PREFACTOR;
	//! ref GlobalOptimizationJob
	static const std::string GLOBAL_OPTIMIZATION;
	//! ref StructureComparisonJob
	static const std::string STRUCTURE_COMPARISON;
	//! ref MonteCarloJob
        static const std::string MONTE_CARLO;
	//! ref TestJob
        static const std::string TEST;
	//! Grabs the correct job as specified by the parameters
	/*!
 	  \param *parameters defined by the config.init file
	*/
        static Job *getJob(Parameters *parameters);
};
#endif
