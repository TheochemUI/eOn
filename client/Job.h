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
    private:
        // make const
        JobType jtype;
    public:
        Job(Parameters *parameters) : jtype{parameters->job} {}
        virtual ~Job() = default;
        //! Virtual run; used solely for dynamic dispatch
        virtual std::vector<std::string> run() = 0;
        JobType getType() { return this->jtype; };
};

namespace helper_functions {
    Job *makeJob(Parameters *params);
} // namespace helper_functions
#endif
