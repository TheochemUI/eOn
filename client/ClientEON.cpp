//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Constants.h"
#include "Parameters.h"
#include "Job.h"
#include "ProcessSearchJob.h"
#include "MinimizationJob.h"
#include "HessianJob.h"
#include "ParallelReplicaJob.h"
#include "ReplicaExchangeJob.h"
#include "DimerDrJob.h"
#include "DimerRotationJob.h"
#include "DisplacementSamplingJob.h"
#include "TestJob.h"

#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <time.h>

//Includes for FPE trapping
#ifdef OSX
    #include <xmmintrin.h>
    #include <mach/mach_init.h>
    #include <mach/task.h>
#endif
#ifdef LINUX
    #include <fenv.h>
#endif
#ifdef WIN32
    #include <float.h>
#endif


#ifndef WIN32
    #include <sys/time.h>
    #include <sys/resource.h>
    #include <sys/utsname.h>
#endif

Parameters parameters;

#ifdef BOINC
    #include <boinc/boinc_api.h>
    #include <boinc/diagnostics.h>     // boinc_init_diagnostics()
    #include <boinc/filesys.h>         // boinc_fopen(), etc...
#ifdef WIN32
    #include <boinc/boinc_win.h>
    #include <boinc/win_util.h>
#endif
#else
    #include "false_boinc.h"
#endif

#ifdef BOINC
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
        }else if (strstr(filename+strlen(filename)-3, "con") == NULL && 
                  strstr(filename+strlen(filename)-3, "dat") == NULL) {
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

    // Floating Point Trapping. It is platform specific!
    // This causes the program to crash on divison by zero,
    // invalid operations, and overflows.
    #ifdef LINUX
        feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    #endif
    #ifdef OSX
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK()
                               & ~_MM_MASK_INVALID 
                               & ~_MM_MASK_DIV_ZERO
                               & ~_MM_MASK_OVERFLOW);
    #endif 
    #ifdef WIN32
        unsigned int cw;
        cw  = _controlfp(0,0) & _MCW_EM;
        cw &= ~(_EM_INVALID|_EM_ZERODIVIDE|_EM_OVERFLOW);
        _controlfp(cw,_MCW_EM);
    #endif

    #ifdef WIN32
    time_t beginTime = time(NULL);
    #else
    struct timeval beginTime;
    gettimeofday(&beginTime, NULL);
    #endif

    #ifdef BOINC
    //We want to uncompress our input file
    char resolved[STRING_SIZE];
    rc = boinc_resolve_filename(BOINC_INPUT_ARCHIVE, resolved, sizeof(resolved));
    if (rc) {
        // 
        fprintf(stderr, "error: cannot resolve file %s\n", BOINC_INPUT_ARCHIVE);
        boinc_finish(rc);
    };
    if (extract_archive(resolved) != 0) {
        printf("error extracting input archive\n");
        boinc_finish(1);
    }
    #endif

    int bundleSize = getBundleSize();

    #ifndef NDEBUG
    printf("Bundle size of %i\n", bundleSize);
    #endif

    int error = parameters.load("parameters_passed.dat");
    if (error) {
        fprintf(stderr, "problem loading parameters file:%s\n", strerror(errno));
        boinc_finish(1);
    }

    FILE *stdoutdat;
    if (parameters.saveStdout == true) {
        //This is kinda a hack, but its for debugging only.
        //A bundle size of 1 must be used for this to work.
        stdoutdat = freopen("stdout_0.dat", "w", stdout);
    }

    // System Information
    #ifdef WIN32
    printf("Windows\n");
    #else
    struct utsname systemInfo;
    int status = uname(&systemInfo);
    if (status == 0) {
        printf("%s %s %s %s %s\n", 
               systemInfo.sysname, systemInfo.nodename, systemInfo.release,
               systemInfo.version, systemInfo.machine);
    }else{
        printf("unknown\n");
    }
    #endif


    // Determine what type of job we are running according 
    // to the parameters file. 
    Job *job=NULL;

    if (parameters.jobType == Job::PROCESS_SEARCH) {
        job = new ProcessSearchJob(&parameters);
    }else if (parameters.jobType == Job::MINIMIZATION) {
        job = new MinimizationJob(&parameters);
    }else if (parameters.jobType == Job::HESSIAN) {
        job = new HessianJob(&parameters);
    }else if (parameters.jobType == Job::PARALLEL_REPLICA) {
        job =  new ParallelReplicaJob(&parameters);
    }else if (parameters.jobType == Job::REPLICA_EXCHANGE) {
        job =  new ReplicaExchangeJob(&parameters);
    }else if (parameters.jobType == Job::DIMER_DR) {
        job =  new DimerDrJob(&parameters);
    }else if (parameters.jobType == Job::DIMER_ROTATION) {
        job =  new DimerRotationJob(&parameters);
    }else if (parameters.jobType == Job::DISPLACEMENT_SAMPLING) {
        job =  new DisplacementSamplingJob(&parameters);
    }else if (parameters.jobType == Job::TEST) {
        job =  new TestJob(&parameters);
    }

    // If no bundles run once; otherwise, run bundleSize number of times.
    if (bundleSize == 0) {
        job->run(-1);
        boinc_fraction_done(1.0);
    }else{
        for (int i=0; i<bundleSize; i++) {
            char buff[100];
            snprintf(buff, 100, "parameters_passed_%i.dat", i);
            string parametersFilename(buff);
            parameters.load(parametersFilename);
            job->run(i);
            boinc_fraction_done((double)(i+1)/(bundleSize));
        }
    }

    delete job;

    // Timing Information
    double utime=0, stime=0, rtime=0;
    #ifdef WIN32
    time_t endTime = time(NULL);
    time_t realTime = endTime-beginTime;
    rtime = (double)realTime;
    #else
    struct timeval endTime;
    gettimeofday(&endTime, NULL);
    rtime = (double)(endTime.tv_sec-beginTime.tv_sec) + 
            (double)(endTime.tv_usec-beginTime.tv_usec)/1000000.0;

    struct rusage r_usage;
    if (getrusage(RUSAGE_SELF, &r_usage)!=0) {
        fprintf(stderr, "problem getting usage info: %s\n", strerror(errno));
    }
    utime = (double)r_usage.ru_utime.tv_sec + (double)r_usage.ru_utime.tv_usec/1000000.0;
    stime = (double)r_usage.ru_stime.tv_sec + (double)r_usage.ru_stime.tv_usec/1000000.0;
    #endif

    printf("\ntiming information:\nreal %10.3f seconds\nuser %10.3f seconds\nsys  %10.3f seconds\n",
           rtime,utime,stime);

    #ifdef OSX
    task_t task = MACH_PORT_NULL;
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
       TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    {
        return -1;
    }
    unsigned int rss = t_info.resident_size;
    unsigned int vs  = t_info.virtual_size;
    printf("\nmemory usage:\nresident size (MB): %8.2f\nvirtual size (MB):  %8.2f\n", (double)rss/1024/1024, (double)vs/1024/1024);
    #endif
    if (parameters.saveStdout == true) {
        fclose(stdoutdat);
    }

    #ifdef BOINC
    //XXX: Error handling!
    rc = boinc_resolve_filename(BOINC_RESULT_ARCHIVE, resolved, sizeof(resolved));
    create_archive(resolved, ".", result_pattern); 
    #endif

    boinc_finish(0);
}
