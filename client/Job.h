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
 * Declaration of job class
 */

class Job { 
    public:
	
	//! job deconstructor.
	virtual ~Job() {};
	//! Virtual run; used solely for dynamic dispatch 
	virtual std::vector<std::string> run()=0;
	//! Grabs the correct job as specified by the parameters
	/*!
 	  \param *parameters defined by the config.init file
	*/
        static Job *getJob(Parameters *parameters);
};
#endif
