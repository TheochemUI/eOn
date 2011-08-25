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
    #include <Python.h>
    #include <stdlib.h>
    #include <sstream>
#endif

#ifdef EONMPIBGP
    #include <libgen.h>
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
//    #ifdef WIN32
//        unsigned int cw;
//        cw  = _controlfp(0,0) & _MCW_EM;
//        cw &= ~(_EM_INVALID|_EM_ZERODIVIDE|_EM_OVERFLOW);
//        _controlfp(cw,_MCW_EM);
//    #endif
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

int main(int argc, char **argv) 
{
    Parameters parameters;

    #ifdef EONMPI
        bool client_standalone=false;
        if (getenv("EON_CLIENT_STANDALONE") != NULL) {
            client_standalone = true;
        }
        char *eon_server = NULL;
        int number_of_clients;
        if (!client_standalone) {
            if (getenv("EON_SERVER") == NULL) {
                fprintf(stderr, "error: must set the env var EON_SERVER_PATH\n");
                return 1;
            }else{
                eon_server = getenv("EON_SERVER");
            }
            if (getenv("EON_NUMBER_OF_CLIENTS") == NULL) {
                fprintf(stderr, "error: must set the env var EON_NUMBER_OF_CLIENTS\n");
                return 1;
            }
            number_of_clients = atoi(getenv("EON_NUMBER_OF_CLIENTS"));
        }else{
            number_of_clients = 1;
        }

        if (MPI::Is_initialized() == false) {
            MPI::Init();
        }

        int error;
        if (client_standalone) {
            error = parameters.load("config_passed.ini");
        }else{
            error = parameters.load("config.ini");
        }
        if (error) {
            fprintf(stderr, "\nproblem loading config.ini file\n");
        }

        //XXX: Barrier for gpaw-python
        MPI::COMM_WORLD.Barrier();

        int irank = MPI::COMM_WORLD.Get_rank();
        int isize = MPI::COMM_WORLD.Get_size();
        
        int *process_types = new int[isize];
        int process_type;

        process_type = 1;

        MPI::COMM_WORLD.Allgather(&process_type,     1, MPI::INT, 
                                  &process_types[0], 1, MPI::INT);

        int i, servers=0, clients=0, potentials=0;
        int server_rank=-1;
        int my_client_number=-1;
        std::vector<int> client_ranks;
        for (i=0;i<isize;i++) {
            switch (process_types[i]) {
                case 0:
                    servers++;
                    break;
                case 1:
                    if (i==irank) {
                        my_client_number = clients;
                    }
                    clients++;
                    client_ranks.push_back(i);
                    break;
                case 2:
                    potentials++;
                    break;
            }
        }

        if (clients < number_of_clients) {
            fprintf(stderr, "didn't launch as many mpi client ranks as"
                            "specified in EON_NUMBER_OF_CLIENTS\n");
            return 1;
        }
        clients = number_of_clients;

        if (parameters.potential == "mpi") {
            int *potential_ranks = new int[potentials];
            int j;
            for (i=0,j=0;i<isize;i++) {
                if (process_types[i] == 2) {
                    potential_ranks[j] = i;
                    j++;
                }
            }
            int potential_group_size = potentials/clients;

            for (i=0;i<clients;i++) {
                MPI::Group orig_group, new_group;
                orig_group = MPI::COMM_WORLD.Get_group();
                int offset = i*potential_group_size;
                new_group = orig_group.Incl(potential_group_size, 
                                            &potential_ranks[offset]);
                MPI::COMM_WORLD.Create(new_group);
            }

            if (my_client_number < number_of_clients) {
                parameters.MPIPotentialRank = potential_ranks[my_client_number*potential_group_size];
                //fprintf(stderr,"MPIPotentialRank: %i\n", parameters.MPIPotentialRank);
            }
        }

        if (!client_standalone) {
            server_rank = client_ranks[number_of_clients];
            if (my_client_number == number_of_clients) {
                std::ostringstream oss;
                oss << client_ranks.at(0);
                for (i=1;i<number_of_clients;i++) {
                    oss << ":" << client_ranks.at(i);
                }
                setenv("EON_CLIENT_RANKS", oss.str().c_str(), 1);
                char **py_argv = (char **)malloc(sizeof(char **)*2);
                py_argv[0] = argv[0];
                py_argv[1] = getenv("EON_SERVER");
                //fprintf(stderr, "rank: %i becoming %s\n", irank, py_argv[1]);
                Py_Initialize();
                Py_Main(2, py_argv);
                Py_Finalize();
                MPI::Finalize();
                return 0;
            }else if (my_client_number > number_of_clients) {
                MPI::Finalize();
                return 0;
            }
        }
    #endif

    #ifndef EONMPI
        if (argc > 1) {
            commandLine(argc, argv);
            return 0;
        }
    #endif

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
    //int outFd = open("/dev/null", O_WRONLY);
    int outFd = open(logfilename, O_WRONLY|O_CREAT|O_TRUNC, 0644);
    if (!client_standalone) {
        dup2(outFd, 1);
    }
    dup2(outFd, 2);
    char *orig_path = new char[1024];
    getcwd(orig_path, 1024);
    while (true) {
        chdir(orig_path);
        char *path = new char[1024];
        int ready=1;
        if (!client_standalone) {
            printf("client: is ready, posting Send to server rank: %i!\n", server_rank);
            //Tag "0" is tell communicator we are ready
            MPI::COMM_WORLD.Isend(&ready,      1, MPI::INT,  server_rank, 0);
            //Tag "1" is to tell the main akmc loop that a client is ready
            MPI::COMM_WORLD.Isend(&ready,      1, MPI::INT,  server_rank, 1);
            MPI::COMM_WORLD.Recv(&path[0], 1024, MPI::CHAR, server_rank, 0);
            if (strncmp("STOPCAR", path, 1024) == 0) {
                MPI::Finalize();
                return 0;
            }
            printf("client: rank: %i chdir to %s\n", irank, path);
        
            if (chdir(path) == -1) {
                fprintf(stderr, "error: %s\n", strerror(errno));
            }
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
            printf("error: Unknown job: %s\n", parameters.job.c_str());
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
    if (client_standalone) {
        break;
    }
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
    if (client_standalone) {
        MPI::COMM_WORLD.Abort(0);
    }else{
        MPI::Finalize();
    }
    #endif

    boinc_finish(0);
}
