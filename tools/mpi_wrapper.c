#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <strings.h>

int main (int argc, char *argv[])
{
    int rank, size; 
    MPI_Init (&argc, &argv);	
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    int len = strlen(argv[1]) + strlen(argv[rank + 2]) + 6;
    char buffer[len];
    snprintf(buffer, sizeof(buffer), "cd %s; %s", argv[rank + 2], argv[1]);
    system(buffer);
    
    MPI_Finalize();
    return 0;
}

