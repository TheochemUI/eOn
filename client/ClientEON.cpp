/*
 *===============================================
 *  EON Client
 *===============================================
 */


#include "Constants.h"
#include "Parameters.h"
#include "Job.h"
#include "ProcessSearch.h"
#include "Minimization.h"

#include <dirent.h>

#ifdef BOINC
    #include <boinc/boinc_api.h>
    #include <boinc/diagnostics.h>     // boinc_init_diagnostics()
    #include <boinc/filesys.h>         // boinc_fopen(), etc...
    #ifdef WIN32
        #include <boinc_win.h>
        #include <win_util.h>
    #endif
#else
    #include "false_boinc.h"
#endif

#ifdef COMPRESSION
#include "Compression.h"
const char BOINC_INPUT_ARCHIVE[] = "input.tgz";
const char BOINC_RESULT_ARCHIVE[] = "result.tgz";

// This will match any file that doesn't have the string "passed",
// has either a "con" or "dat" in it, and doesn't start with a period.
int result_pattern(char *filename)
{
        if (strstr(filename, "passed") != NULL) {
            return 0; 
        }else if (filename[0] == '.') {
            return 0;
        }else if (strstr(filename, "con") == NULL && 
                  strstr(filename, "dat") == NULL) {
            return 0;
        }
        return 1;
}
#endif

int getBundleSize(void) {
    DIR *dir;
    struct dirent *dp;
    int num_bundle=0;

    dir = opendir(".");

    while ((dp=readdir(dir))) {
        if (dp->d_name[0] == '.') {
            continue;
        }
        if (strstr(dp->d_name, "passed")==NULL) {
            continue;
        }
        char *ch = strrchr(dp->d_name, '_')+1;
        char *cch = strrchr(dp->d_name, '.');
        *cch = '\0';
        if (isdigit(*ch)) {
            int i=atoi(ch)+1;
            if (i>num_bundle) {
                num_bundle = i;
            }
        }
    }

    return num_bundle;
}

int main(int argc, char **argv) 
{
    int rc;
	// BOINC is started
	rc = boinc_init();
	if(rc){
		boinc_finish(rc);
	}

    #ifdef COMPRESSION
    //We want to uncompress our input file
    char resolved[STRING_SIZE];
    rc = boinc_resolve_filename(BOINC_INPUT_ARCHIVE, resolved, sizeof(resolved));
    if (extract_archive(resolved) != 0) {
        printf("error extracting input archive\n");
        boinc_finish(1);
    }
    #endif

    int bundleSize = getBundleSize();
    #ifndef NDEBUG
    printf("Bundle size of %i\n", bundleSize);
    #endif

    // store all the runtime parameters received from the server
    Parameters parameters; 
    string parameters_passed("parameters_passed.dat");
    int error = parameters.load(parameters_passed);
    if (error) {
        boinc_finish(1);
    }
 
    // Initialize random generator
    if(parameters.randomSeed < 0)
    {
        unsigned i = time(NULL);
        parameters.randomSeed = i;
        helper_functions::random(i);
    }else{
        helper_functions::random(parameters.randomSeed);
    }
    printf("Random seed is: %ld\n", parameters.randomSeed);

    // Determine what type of job we are running according 
    // to the parameters file. 
    Job *job;

    if (parameters.jobType == Parameters::PROCESS_SEARCH) {
        job = new ProcessSearch(&parameters);
    }else if (parameters.jobType == Parameters::MINIMIZATION) {
        job = new Minimization(&parameters);
    }

    //If no bundles run once; otherwise, run bundleSize number of times.
    if (bundleSize == 0) {
        job->run(-1);
        boinc_fraction_done(1.0);
    }else{
        for (int i=0; i<bundleSize; i++) {
            job->run(i);
            boinc_fraction_done((double)(i+1)/(bundleSize));
        }
    }

    delete job;

    #ifdef COMPRESSION
    //XXX: Error handling!
    rc = boinc_resolve_filename(BOINC_RESULT_ARCHIVE, resolved, sizeof(resolved));
    create_archive(resolved, ".", result_pattern); 
    #endif

    boinc_finish(0);
}
