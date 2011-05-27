#include <Python.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>

void run_pyscript(char *scriptname);
void usage(void);
int client_main(int, char**);

void wrapper_usage()
{
    fprintf(stderr, "wrapper <eon.py> <gpaw_sp.py> <# eonclients> <# GPAWs>\n");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int irank, isize;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    MPI_Comm_size(MPI_COMM_WORLD, &isize);

    if (argc<4 || argc>5) {
        fprintf(stderr,"incorrect number of arguments\n");
        wrapper_usage();
        return 1;
    }

    int   client_number = atoi(argv[3]);
    int   gpaw_number   = atoi(argv[4]);
    char *gpaw_path     = argv[2];
    char *server_path   = argv[1];

    if (irank < gpaw_number) {
        printf("rank: %i becoming potential\n", irank);
        run_pyscript(gpaw_path);
    }else if (irank < gpaw_number+client_number) {
        printf("rank: %i becoming client\n", irank);
        client_main(0, NULL);
    }else if (irank == gpaw_number+client_number) {
        char *eonpath = getenv("EONPATH");
        char *pythonpath = getenv("PYTHONPATH");
        char newpath[2048];
        snprintf(newpath, 2048, "PYTHONPATH=%s:%s", eonpath, pythonpath);
        printf("setting path: %s\n", newpath);
        putenv(newpath);
        printf("rank: %i becoming server\n", irank);
        run_pyscript(server_path);
    }

    MPI_Finalize();
    return 0;
}

void run_pyscript(char *scriptname) 
{
    Py_Initialize();
    char *argv[] = { basename(scriptname) };
    printf("running script: %s\n", argv[0]);
    PySys_SetArgv(1, argv);
    printf("python path: %s\n", Py_GetPath());
    FILE *fp = fopen(scriptname, "r");
    PyRun_SimpleFile(fp, scriptname);
    Py_Finalize();
}
