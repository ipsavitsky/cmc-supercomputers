#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 6

void report(const char *msg, int coords[2]) {
    if (coords[0] == 0 && coords[1] == 0) {
        printf("(0,0) %s!\n", msg);
    }
}

void send_coords_and_value(int coords[2], int value, int other_coords[2],
                           int other_rank, MPI_Comm comm) {
    MPI_Send(coords, 2, MPI_INT, other_rank, 0, comm);
    MPI_Send(&value, 1, MPI_INT, other_rank, 0, comm);
}

void receive_coords_and_value(int coords[2], int *value, int other_coords[2],
                              int other_rank, MPI_Comm comm) {
    MPI_Recv(other_coords, 2, MPI_INT, other_rank, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(value, 1, MPI_INT, other_rank, 0, comm, MPI_STATUS_IGNORE);
}

void collide_rows(int row_1, int row_2, int coords[2], int *a,
                  int best_coords[2], MPI_Comm comm) {
    int result = 0;
    int recieved_coords[2];
    int other_rank = 0;
    int other_coords[2];
    other_coords[1] = coords[1];
    if (coords[0] == row_1 || coords[0] == row_2) {
        other_coords[0] = coords[0] == row_1 ? coords[0] + 1 : coords[0] - 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        send_coords_and_value(best_coords, *a, other_coords, other_rank, comm);
    } else if (coords[0] == row_1 + 1 || coords[0] == row_2 - 1) {
        other_coords[0] =
            coords[0] == row_1 + 1 ? coords[0] - 1 : coords[0] + 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        receive_coords_and_value(recieved_coords, &result, recieved_coords,
                                 other_rank, comm);
        if (result > *a) {
            *a = result;
            best_coords[0] = recieved_coords[0];
            best_coords[1] = recieved_coords[1];
        }
    }
}

void compress_row(int row_n, int pl_l, int pl_r, int coords[2], int *a, int best_coords[2],
                  MPI_Comm comm) {
    int result = 0;
    int recieved_coords[2];
    int other_rank = 0;
    int other_coords[2];

    other_coords[0] = coords[0];

    if (coords[0] == row_n && (coords[1] == pl_l || coords[1] == pl_r)) {
        other_coords[1] = coords[1] == 0 ? coords[1] + 1 : coords[1] - 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        send_coords_and_value(best_coords, *a, other_coords, other_rank, comm);
    }
    if (coords[0] == row_n &&
        (coords[1] == pl_l + 1 || coords[1] == pl_r - 1)) {
        other_coords[1] = coords[1] == 1 ? coords[1] - 1 : coords[1] + 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        receive_coords_and_value(recieved_coords, &result, recieved_coords,
                                 other_rank, comm);
        if (result > *a) {
            *a = result;
            best_coords[0] = recieved_coords[0];
            best_coords[1] = recieved_coords[1];
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, tasks;
    MPI_Comm comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);

    int size[2] = {SIZE, SIZE};
    int periodic[2] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 0, &comm);
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    srand(rank + 6);
    int a = rand() % 1000;
    printf("Coordinates for process %d: (%d, %d)\n", rank, coords[0],
           coords[1]);
    printf("a[%d][%d] = %d\n", coords[0], coords[1], a);

    int result = 0;
    int other_coords[2];
    int recieved_coords[2];
    int best_coords[2];
    best_coords[0] = coords[0];
    best_coords[1] = coords[1];
    int other_rank = 0;

    other_coords[1] = coords[1];
    collide_rows(0, 5, coords, &a, best_coords, comm);
    MPI_Barrier(comm);
    report("step 1", coords);

    collide_rows(1, 4, coords, &a, best_coords, comm);
    MPI_Barrier(comm);
    report("step 2", coords);

    if (coords[0] == 3) {
        other_coords[0] = coords[0] - 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        send_coords_and_value(best_coords, a, other_coords, other_rank, comm);
    }
    if (coords[0] == 2) {
        other_coords[0] = coords[0] + 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        receive_coords_and_value(recieved_coords, &result, recieved_coords,
                                 other_rank, comm);
        if (result > a) {
            a = result;
            best_coords[0] = recieved_coords[0];
            best_coords[1] = recieved_coords[1];
        }
    }
    MPI_Barrier(comm);
    report("step 3", coords);

    if (coords[0] == 2 && (coords[1] == 0 || coords[1] == 5)) {
        other_coords[0] = coords[0];
        other_coords[1] = coords[1] == 0 ? coords[1] + 1 : coords[1] - 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        send_coords_and_value(best_coords, a, other_coords, other_rank, comm);
    }
    if (coords[0] == 2 && (coords[1] == 1 || coords[1] == 4)) {
        other_coords[0] = coords[0];
        other_coords[1] = coords[1] == 1 ? coords[1] - 1 : coords[1] + 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        receive_coords_and_value(recieved_coords, &result, recieved_coords,
                                 other_rank, comm);
        if (result > a) {
            a = result;
            best_coords[0] = recieved_coords[0];
            best_coords[1] = recieved_coords[1];
        }
    }
    MPI_Barrier(comm);
    report("step 4", coords);

    if (coords[0] == 2 && (coords[1] == 1 || coords[1] == 4)) {
        other_coords[0] = coords[0];
        other_coords[1] = coords[1] == 1 ? coords[1] + 1 : coords[1] - 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        send_coords_and_value(best_coords, a, other_coords, other_rank, comm);
    }
    if (coords[0] == 2 && (coords[1] == 2 || coords[1] == 3)) {
        other_coords[0] = coords[0];
        other_coords[1] = coords[1] == 2 ? coords[1] - 1 : coords[1] + 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        receive_coords_and_value(recieved_coords, &result, recieved_coords,
                                 other_rank, comm);
        if (result > a) {
            a = result;
            best_coords[0] = recieved_coords[0];
            best_coords[1] = recieved_coords[1];
        }
    }
    MPI_Barrier(comm);
    report("step 5", coords);

    if (coords[0] == 2 && coords[1] == 3) {
        other_coords[0] = coords[0];
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        send_coords_and_value(best_coords, a, other_coords, other_rank, comm);
    }
    if (coords[0] == 2 && coords[1] == 2) {
        other_coords[0] = coords[0];
        other_coords[1] = 3;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        receive_coords_and_value(recieved_coords, &result, recieved_coords,
                                 other_rank, comm);
        if (result > a) {
            a = result;
            best_coords[0] = recieved_coords[0];
            best_coords[1] = recieved_coords[1];
        }
    }
    MPI_Barrier(comm);
    report("step 6", coords);
    if (coords[0] == 2 && coords[1] == 2) {
        printf("Max result: %d on (%d, %d)\n", a, best_coords[0],
               best_coords[1]);
    }
    MPI_Finalize();
    return 0;
}
