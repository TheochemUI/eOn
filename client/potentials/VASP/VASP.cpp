#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#ifdef WIN32
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
//#define popen _popen
#else
#include <sys/wait.h>
#include <fcntl.h>
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
    #include "../../false_boinc.h"
#endif

#include "VASP.h"
bool  VASP::firstRun = true;
long  VASP::vaspRunCount = 0;
pid_t VASP::vaspPID = 0;

VASP::VASP(void)
{
	vaspRunCount++;
    return;
}

void VASP::cleanMemory(void)
{
	vaspRunCount--;
	if(vaspRunCount < 1) {
		FILE *stopcar = fopen("STOPCAR", "w");
		fprintf(stopcar, "LABORT = .TRUE.");
		fclose(stopcar);
	}
    return;
}

VASP::~VASP()
{
	cleanMemory();
}

void VASP::spawnVASP()
{
    if ((vaspPID=fork()) == -1) {
        fprintf(stderr, "error forking for vasp: %s\n", strerror(errno));
        boinc_finish(1);
    }

    if (vaspPID) {
        /* We are the parent */
        setvbuf(stdout, (char*)NULL, _IONBF, 0); //non-buffered output
    }else{
        /* We are the child */
        int outFd = open("vaspout", O_APPEND|O_CREAT|O_WRONLY);
        dup2(outFd, 1);
        dup2(outFd, 2);

        char vaspPath[1024];
        if(boinc_resolve_filename("vasp", vaspPath, 1024)) {
            fprintf(stderr, "problem resolving vasp filename\n");
            boinc_finish(1);
        }

        if (execlp(vaspPath, "vasp", NULL) == -1) {
            fprintf(stderr, "error spawning vasp: %s\n", strerror(errno));
            boinc_finish(1);
        }
    }
}

bool VASP::vaspRunning()
{
    pid_t pid;
    int status;

    if (vaspPID == 0) {
        return false;
    }

    pid = waitpid(vaspPID, &status, WNOHANG);

    if (pid) {
        fprintf(stderr, "vasp died unexpectedly!\n");
        boinc_finish(1);
    }

    return true;
}


void VASP::force(long N, const double *R, const int *atomicNrs, double *F, 
                 double *U, const double *box)
{
    writeNEWCAR(N, R, atomicNrs, box);

    if(!vaspRunning())
    {
        spawnVASP();
    }

	printf("vasp force call");
	fflush(stdout);
    while(access("FU", F_OK) == -1)
    {
        sleep(1);
		printf(".");
		fflush(stdout);
		vaspRunning();
            
    }
	printf("\n");
    readFU(N, F, U);
    remove("FU");
    vaspRunCount++;
    return;
}


void VASP::writeNEWCAR(long N, const double *R, const int *atomicNrs, 
                       const double *box)
{
    // Positions are scaled 
    long i = 0;
    long i_old = 0;
    FILE *NEWCAR;
    
    if(firstRun) {
        NEWCAR = fopen("POSCAR","w");
        firstRun = false;
    }else{
        NEWCAR = fopen("NEWCAR","w");
    }

    // header line (treated as a comment)
    i_old = 0;
    fprintf(NEWCAR, "%d ", atomicNrs[0]);
    for(i = 0; i < N; i++)
    {
        if(atomicNrs[i] != atomicNrs[i_old])
        {
            fprintf(NEWCAR, "%d ", atomicNrs[i]);
            i_old = i;
        }
    }
    fprintf(NEWCAR, ": Atomic numbers\n");
    
    // boundary box
    fprintf(NEWCAR, "1.0\n");
    fprintf(NEWCAR, " %.8f\t%.8f\t%.8f\n", box[0], 0.0, 0.0);
    fprintf(NEWCAR, " %.8f\t%.8f\t%.8f\n", 0.0, box[1], 0.0);
    fprintf(NEWCAR, " %.8f\t%.8f\t%.8f\n", 0.0, 0.0, box[2]);

    // the number of atoms of each of the the different atomic types
    i_old = 0;
    for(i = 0; i < N; i++)
    {
        if(atomicNrs[i] != atomicNrs[i_old])
        {
            fprintf(NEWCAR, "%li ", i - i_old);
            i_old = i;
        }
    }
    fprintf(NEWCAR, "%li\n", N - i_old);

    // coordinates for all atoms
    fprintf(NEWCAR, "Cartesian\n");
    for(i = 0; i < N; i++)
    {
        fprintf(NEWCAR, "%.19f\t%.19f\t%.19f\t T T T\n", R[i * 3 + 0], R[i * 3 + 1],  R[i * 3 + 2]);
    }
    fclose(NEWCAR);
    return;
}


void VASP::readFU(long N, double *F, double *U)
{
    FILE *FU;
    FU = fopen("FU", "r");
    
    fscanf(FU, "%lf", U);
    
    for(int i = 0; i < N; i++)
    {
        fscanf(FU, "%lf %lf %lf", &F[i * 3 + 0], &F[i * 3 + 1], &F[i * 3 + 2]);
    }
    fclose(FU);
    return;
}
