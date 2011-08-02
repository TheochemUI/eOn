#include <Python.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>

#warning "Building Combined server/client binary for use with MPMD"

void usage(void);
int client_main(int, char**);

void wrapper_usage()
{
    fprintf(stderr, "wrapper <eon.py> <# eonclients>\n");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int irank, isize;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    MPI_Comm_size(MPI_COMM_WORLD, &isize);
    printf("eonbgp: rank: %i\n", irank);

    if (argc != 2) {
        fprintf(stderr,"incorrect number of arguments: %i\n", argc);
        wrapper_usage();
        return 1;
    }

    int   client_number = atoi(argv[2]);
    char *server_path   = argv[1];

    if (irank == 0) {
        printf("rank: %i becoming server\n", irank);
        Py_Initialize();
        Py_Main(2, argv);
        Py_Finalize();
        MPI_Finalize();
    }else{
        printf("rank: %i becoming client\n", irank);
        client_main(0, NULL);
        //client MPI_Finalizes its self
    }

    return 0;
}
