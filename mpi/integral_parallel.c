#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[]) {
    int n, myid, numprocs, i;
    double pi, h, sum, x, from = 0, to = 10000, start, stop;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    n = strtol(argv[1], NULL, 10);
    start = MPI_Wtime();
    h = (to - from) / (double)n;
    sum = 0.0;
    for (i = myid + 1; i <= n; i += numprocs) {
        x = h * ((double)i - 0.5);
        sum += sqrt(x);
    }
    sum = h * sum;
    MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (myid == 0) printf("%d,%lf,%d,%lf\n", numprocs, pi, n, stop-start);
    MPI_Finalize();
    return 0;
}
