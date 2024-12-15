// Each process prints its ID and and the total number of processes.

#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
	int id, numProcs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	printf("I am process %d out of %d processes\n", id, numProcs);
	MPI_Finalize();
}
