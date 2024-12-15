#include <stdio.h>
#include <mpi.h>

#define npp 2

int main(int argc, char *argv[]) {
    double t1, t2;
    int id, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int *S = new int[npp];
    int *M;
    int *R;
    for (int i=0; i<npp; i++) S[i] = id;

    if (id == 0) {
        M = new int[numProcs*npp];
        R = new int[numProcs];
        for (int i=0; i<npp*numProcs; i++) M[i] = i;
    }

    if (id == 0) {
        printf("Before:\n");
        for (int i=0; i<npp*numProcs; i++) printf("%d ", M[i]);
        printf("\n");
        //MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Scatter(M, npp, MPI_INT, S, npp, MPI_INT, 0, MPI_COMM_WORLD);

    int *x=new int[1];
    x[0]=0;
    for(int i=0;i<npp;i++){
        x[0]+=S[i];
    }

    MPI_Gather(x, 1, MPI_INT, R, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int sum=0;
    if (id == 0) {
        printf("After:\n");
        for (int i=0; i<numProcs; i++){
                printf("%d ", R[i]);
                sum+=R[i];
        }
        printf("\n");
        //MPI_Barrier(MPI_COMM_WORLD);
         printf("sum: %d", sum);
    }



    MPI_Finalize();
}
