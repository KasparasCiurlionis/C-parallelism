#include <mpi.h>

int main(int argc, char *argv[]) {
        int id, numProcs, buff;

	MPI_Status stat;
        MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int token=0;
    int i=1;
    if(id==0){
        while(token < 100){
            if(i==numProcs){
                i=1;
            }
            MPI_Send(&token, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Recv(&token, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            printf("token: %d\n", token);
            i++;
        }
        if(token >= 100){
            for(int i=1;i<numProcs;i++){
                MPI_Send(&token, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
    }
    else{
        while(token < 100){
        MPI_Recv(&token, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        if(token>=100){
            break;
        }
        token+=id;
        MPI_Send(&token, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
}
