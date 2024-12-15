#include <stdio.h>
#include <mpi.h>

#define N 5

int main(int argc, char *argv[]) {
    double t1, t2;
    int id, numProcs;
    int *A = new int[N];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (id == 0) {
         for (int i=0; i<N; i++)
             A[i] = i;
    }

    if (id == 0) printf("Before:\n");
    for (int i=0; i<numProcs; i++) {
        if (id == i) {
           printf("Proc %d: ", id);
           for (int i=0; i<N; i++) printf("%d ", A[i]);
           printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Bcast(A, N, MPI_INT, 0, MPI_COMM_WORLD);

    if (id == 0) printf("After:\n");
    for (int i=0; i<numProcs; i++) {
        if (id == i) {
           printf("Proc %d: ", id);
           if(id<N){
                printf(A[id]);
           }
           else{
                printf(A[N-1]);
           }
           printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
}
