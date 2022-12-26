#include <math.h>
#include <mpi-ext.h>
#include <mpi.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int *ranks_gc;
int nf = 0;

const float left_br = 0;      /* lower limit of integration */
const float right_br = 10000; /* upper limit of integration */

typedef double (*math_func)(double);

double integrate(math_func f, double a, int num, double h) {
    int myid, numprocs;
    double x, sum = 0.0;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    for (int i = myid + 1; i <= n; i += numprocs) {
        x = h * ((double)i - 0.5);
        sum += f(x);
    }
    return h * sum;
}

static void verbose_errhandler(MPI_Comm *pcomm, int *perr, ...) {
    free(ranks_gc);

    MPI_Comm comm = *pcomm;
    int err = *perr;
    char errstr[MPI_MAX_ERROR_STRING];
    int i, rank, size, len, eclass;
    MPI_Group group_c, group_f;
    int *ranks_gf;

    MPI_Error_class(err, &eclass);
    if (MPIX_ERR_PROC_FAILED != eclass) {
        MPI_Abort(comm, err);
    }

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPIX_Comm_failure_ack(comm);
    MPIX_Comm_failure_get_acked(comm, &group_f);
    MPI_Group_size(group_f, &nf);
    MPI_Error_string(err, errstr, &len);

    ranks_gf = (int *)malloc(nf * sizeof(int));
    ranks_gc = (int *)malloc(nf * sizeof(int));
    MPI_Comm_group(comm, &group_c);
    for (i = 0; i < nf; i++)
        ranks_gf[i] = i;
    MPI_Group_translate_ranks(group_f, nf, ranks_gf, group_c, ranks_gc);

    free(ranks_gf);
}

int main(int argc, char *argv[]) {
    int n, size, i, j, ierr, num;
    double h, result, a, b, pi;
    double my_a, my_range;
    double startwtime, my_time = 0.0, mintime = 0.0, time = 0.0;
    int rank, source, dest, tag, count;
    MPI_Status status;
    double my_result;
    a = 0.;
    b = 1;        
    n = (argc > 1) ? pow(2, atoi(argv[1]))
                   : 512; /* number of increment within each process */

    dest = 0;  /* define the process that computes the final result */
    tag = 123; /* set the tag to identify this particular job */

    /* Starts MPI processes ... */

    // int rank, size;
    MPI_Errhandler errh;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm_create_errhandler(verbose_errhandler, &errh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);
    MPI_Barrier(MPI_COMM_WORLD);

    double *vec_r = (double *)malloc(size * sizeof(double));

    h = (right_br - left_br) / n; /* length of increment */
    num = n / size; /* number of intervals calculated by each process*/
    my_range = (right_br - left_br) / size;
    my_a = left_br + rank * my_range;

    startwtime = MPI_Wtime();
    my_result = integrate(sqrt, my_a, num, h);
    my_time = MPI_Wtime() - startwtime;

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == (size - 1) || rank == (size / 2)) {
        printf("Rank : %d, my_a: %f, my_res: %f\n", rank, my_a, my_result);
        printf("Rank %d / %d: bye bye!\n\n", rank, size);
        raise(SIGKILL);
    }

    if (rank == 0) {
        vec_r[0] = my_result;
        time = my_time;

        for (int i = 1; i < size; i++) {
            source = i; /* MPI process number range is [0,size-1] */
            MPI_Recv(&my_result, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
                     &status);

            vec_r[i] = my_result;

            MPI_Recv(&my_time, 1, MPI_DOUBLE, source, tag - 1, MPI_COMM_WORLD,
                     &status);

            time = fmax(time, my_time);
        }
    } else
        MPI_Send(&my_result, 1, MPI_DOUBLE, dest, tag,
                 MPI_COMM_WORLD); /* send my_result to intended dest.*/
    MPI_Send(&my_time, 1, MPI_DOUBLE, dest, tag - 1,
             MPI_COMM_WORLD); /* send my_time to intended dest.*/

    if (rank == 0) {
        if (nf != 0) {
            for (int i = 0; i < nf; i++) {
                my_a = left_br + ranks_gc[i] * my_range;

                startwtime = MPI_Wtime();
                my_result = integrate(sqrt, my_a, num, h);
                my_time = MPI_Wtime() - startwtime;

                printf("Recalculate : %d, my_a: %f, my_res: %f\n", ranks_gc[i],
                       my_a, my_result);
                vec_r[ranks_gc[i]] = my_result;
                time = fmax(time, my_time);
            }
            nf = 0;
        }

        result = 0;
        for (int i = 0; i < size; i++) {
            result += vec_r[i];
        }

        printf("\nRESULT: %f\n", result);
        printf("Time: %f\n", time);
    }

    free(vec_r);

    MPI_Finalize();
    return 0;
}
