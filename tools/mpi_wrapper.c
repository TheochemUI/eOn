#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/time.h>

double mysecond();

int main (int argc, char *argv[])
{
    int rank, size;
    double t1, t2, t3;
    double used, idle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    int len = strlen(argv[1]) + strlen(argv[rank + 2]) + 16;
    char buffer[len];
    snprintf(buffer, sizeof(buffer), "cd %s; %s>/dev/null", argv[rank + 2], argv[1]);

    t1 = mysecond();
    system(buffer);
    t2 = mysecond();

    MPI_Barrier(MPI_COMM_WORLD);
    t3 = mysecond();

    used = t2 - t1;
    idle = t3 - t2;

    double globused = 0.0;
    double globidle = 0.0;

    MPI_Reduce(&used, &globused, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&idle, &globidle, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        printf("Total time used: %0.1lf\n", globused);
        printf("Total time idle: %0.1lf\n", globidle);
        printf("Efficiency: %0.2lf%%\n", (used)/(idle+used)*100.0);
    }
    MPI_Finalize();
    return 0;
}

double mysecond()
{
   struct timeval tp;
   struct timezone tzp;
   int i;
   i = gettimeofday(&tp,&tzp);
   return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
