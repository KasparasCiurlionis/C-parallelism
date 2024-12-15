#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
        int id, numProcs, buff;

	MPI_Status stat;
        MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int token=0;
    if(id==0){
        MPI_Send(&token, 1, MPI_INT, id+1, 0, MPI_COMM_WORLD);
    }
    printf("test\n");
	MPI_Recv(&token, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
	printf("Message received from process %d\n", stat.MPI_SOURCE);
    token+=id;
    if(id!=0){
        if(id<numProcs-1){
            MPI_Send(&token, 1, MPI_INT, id+1, 0, MPI_COMM_WORLD);
        }
        else{
            MPI_Send(&token, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    if(id==0){
        printf("token: %d", token);
    }

    MPI_Finalize();
}

