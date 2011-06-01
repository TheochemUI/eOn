//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Bundling.h"
#include "CommandLine.h"
#include "Constants.h"
#include "Parameters.h"
#include "Job.h"

#include <errno.h>
#include <string.h>
#include <time.h>

#ifdef EONMPI
    #include <mpi.h>
    #include <unistd.h>
    #include <fcntl.h>
#endif

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
#endif

void enableFPE(void)
{
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
}

void printSystemInfo()
{
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
}

#ifndef EONMPIBGP
int main(int argc, char **argv) 
#else
int client_main(int argc, char **argv)
#endif
{
    Parameters parameters;
    #ifdef EONMPI
        if (MPI::Is_initialized() == false) {
            MPI::Init();
        }

        //XXX: Barrier for gpaw-python
        MPI::COMM_WORLD.Barrier();

        int irank = MPI::COMM_WORLD.Get_rank();
        int isize = MPI::COMM_WORLD.Get_size();
        //printf("client rank: %i size: %i\n", irank, isize);
        
        int *process_types = new int[isize];
        int process_type = 1;

        MPI::COMM_WORLD.Allgather(&process_type,     1, MPI::INT, 
                                  &process_types[0], 1, MPI::INT);

        int i, servers=0, clients=0, potentials=0;
        int server_rank=-1;
        int my_client_number=-1;
        for (i=0;i<isize;i++) {
            switch (process_types[i]) {
                case 0:
                    servers++;
                    server_rank = i;
                    break;
                case 1:
                    if (i==irank) {
                        my_client_number = clients;
                    }
                    clients++;
                    break;
                case 2:
                    potentials++;
                    break;
            }
        }

        int potential_group_size = potentials/clients;
        parameters.MPIPotentialRank = potential_group_size*my_client_number;
        int *potential_ranks = new int[potential_group_size];
        int j;
        for (i=0,j=0;i<isize;i++) {
            if (process_types[i] == 2) {
                potential_ranks[j] = i;
                j++;
            }
        }

        for (i=0;i<clients;i++) {
            MPI::Group orig_group, new_group;
            orig_group = MPI::COMM_WORLD.Get_group();
            int offset = i*potential_group_size;
            //printf("client: rank: %i offset: %i pot_group_size: %i\n", irank, 
            //       offset,potential_group_size);
            new_group = orig_group.Incl(potential_group_size, 
                                        &potential_ranks[offset]);
            (void)MPI::COMM_WORLD.Create(new_group);
        }
    #endif

    if (argc > 1) {
        commandLine(argc, argv);
        return 0;
    }

    #ifdef BOINC
    // BOINC is started
    int rc;
    rc = boinc_init();
    if(rc){
        boinc_finish(rc);
    }

    //We want to uncompress our input file
    char resolved[STRING_SIZE];
    rc = boinc_resolve_filename(BOINC_INPUT_ARCHIVE, resolved, sizeof(resolved));
    if (rc) {
        fprintf(stderr, "error: cannot resolve file %s\n", BOINC_INPUT_ARCHIVE);
        boinc_finish(rc);
    };
    if (extract_archive(resolved) != 0) {
        printf("error extracting input archive\n");
        boinc_finish(1);
    }
    #endif

    enableFPE();

    #ifdef WIN32
    time_t beginTime = time(NULL);
    #else
    struct timeval beginTime;
    gettimeofday(&beginTime, NULL);
    #endif

    #ifdef EONMPI
    //XXX: When do we stop? The server should probably tell everyone 
    //     when to stop.
    char logfilename[1024];
    snprintf(logfilename, 1024, "eonclient_%i.log", irank);
    int outFd = open(logfilename, O_CREAT|O_WRONLY|O_TRUNC, 0644);
    dup2(outFd, 1);
    dup2(outFd, 2);
    char *orig_path = new char[1024];
    getcwd(orig_path, 1024);
    while (true) {
        chdir(orig_path);
        char *path = new char[1024];
        int ready=1;
        printf("client: is ready, posting Send!\n");
        MPI::COMM_WORLD.Send(&ready,      1, MPI::INT,  server_rank, 0);
        MPI::COMM_WORLD.Recv(&path[0], 1024, MPI::CHAR, server_rank, 0);
        if (strncmp("STOPCAR", path, 1024) == 0) {
            MPI::Finalize();
            return 0;
        }
        printf("client: rank: %i chdir to %s\n", irank, path);
        
        if (chdir(path) == -1) {
            fprintf(stderr, "error: %s\n", strerror(errno));
        }
    #endif

    printSystemInfo();

    bool bundlingEnabled = true;
    int bundleSize = getBundleSize();
    if (bundleSize == 0) {
        bundleSize = 1;
    }else if (bundleSize == -1) {
        //Not using bundling
        bundleSize = 1;
        bundlingEnabled = false;
    }

    std::vector<std::string> bundledFilenames;
    for (int i=0;i<bundleSize;i++) {

        printf("Beginning Job %d of %d\n", i, bundleSize);
        std::vector<std::string> unbundledFilenames;
        if (bundlingEnabled) {
            unbundledFilenames = unbundle(i);
        }

        int error = parameters.load("config_passed.ini");
        if (error) {
            fprintf(stderr, "\nproblem loading parameters file\n");
            boinc_finish(1);
        }

        // Determine what type of job we are running according 
        // to the parameters file. 
        Job *job = Job::getJob(&parameters);
        if (job == NULL) {
            printf("error: Unknown job: %s\n", parameters.potential.c_str());
            return 1;
        }
        std::vector<std::string> filenames = job->run();

        if (bundlingEnabled) {
            bundle(i, filenames, &bundledFilenames);
            deleteUnbundledFiles(unbundledFilenames);
        }else{
            bundledFilenames = filenames;
        }

        boinc_fraction_done((double)(i+1)/(bundleSize));
        delete job;
    }

    #ifdef EONMPI
    //End of MPI while loop
    }
    #endif

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
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
       TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    {
        return -1;
    }
    unsigned int rss = t_info.resident_size;
    unsigned int vs  = t_info.virtual_size;
    printf("\nmemory usage:\nresident size (MB): %8.2f\nvirtual size (MB):  %8.2f\n",
           (double)rss/1024/1024, (double)vs/1024/1024);
    #endif

    #ifdef BOINC
    //XXX: Error handling!
    rc = boinc_resolve_filename(BOINC_RESULT_ARCHIVE, resolved, sizeof(resolved));
    char dirToCompress[] = ".";
    create_archive(resolved, dirToCompress, bundledFilenames); 
    #endif

    #ifdef EONMPI
        MPI::Finalize();
    #endif

    boinc_finish(0);
}
