#include <Python.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>

void usage(void);
int client_main(int, char**);

void wrapper_usage()
{
    fprintf(stderr, "wrapper <eon.py> <# eonclients> <# GPAWs>\n");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int irank, isize;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    MPI_Comm_size(MPI_COMM_WORLD, &isize);
    printf("eonbgp: rank: %i\n", irank);

    if (argc != 4) {
        fprintf(stderr,"incorrect number of arguments: %i\n", argc);
        wrapper_usage();
        return 1;
    }

    int   client_number = atoi(argv[2]);
    int   gpaw_number   = atoi(argv[3]);
    char *server_path   = argv[1];

    if (irank < gpaw_number) {
        //shouldn't happen
    }else if (irank < gpaw_number+client_number) {
        printf("rank: %i becoming client\n", irank);
        client_main(0, NULL);
    }else if (irank == gpaw_number+client_number) {
        Py_Initialize();
        Py_Main(2, argv);
        Py_Finalize();
    }

    MPI_Finalize();
    return 0;
}
