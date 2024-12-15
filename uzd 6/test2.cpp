#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
        int id, numProcs, buff;

	MPI_Status stat;
        MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int masyv[numProcs];
    for(int i=0;i<numProcs; i++){
        masyv[i]=i;
    }
	if (id == 0) {
		for (int i=0; i<numProcs; i++) {
			MPI_Send(masyv[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Recv(&buff, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    MPI_Send(&buff+1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);


	if (id == 0) {
		for (int i=0; i<numProcs; i++) {
			MPI_Recv(&buff, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            printf("Message received from process %d\n", stat.MPI_SOURCE);
            masyv[stat.MPI_SOURCE]=buff;
		}
	}
        MPI_Finalize();
}

