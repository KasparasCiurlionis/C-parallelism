#include <stdio.h>
#include <mpi.h>

void DoWork(double n) {
    double t = MPI_Wtime();
    while (MPI_Wtime()-t < n);

}


int main(int argc, char *argv[]) {
    int id, numProcs, count, buff;
    double ts, tf;
    MPI_Status stat;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    MPI_Barrier(MPI_COMM_WORLD);

    ts = MPI_Wtime();

    int flag;

    if (id == 0) {
        buff = 2;
        int num=0;
        int p=0;
        for(int i=1;i<numProcs;i++){
            MPI_Issend(&buff, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            //MPI_Wait(&request, &stat);
            MPI_Test(&request, &flag, &stat);
            buff++;
        }
        num=buff;
        while(p<5){
            MPI_Irecv(&buff, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &stat);
            if(buff!=0){
                printf("%d\n", buff);
                p++;
                if(p==5){
                    break;
                }
            }
            buff=num;
            MPI_Issend(&buff, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD, &request);
            num++;
        }
        buff=-1;
        for(int i=1;i<numProcs;i++){
            MPI_Issend(&buff, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            //MPI_Wait(&request, &stat);
            MPI_Test(&request, &flag, &stat);
        }
    }

     if (id != 0) {
        buff=1;
        while(buff!=-1){
            MPI_Irecv(&buff, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &stat);
            if(buff!=2){
                for(int i=2;i<buff;i++){
                    if(buff%i==0){
                        buff=0;
                        break;
                    }
                }
            }
            else if(buff==-1){
                break;
            }
            MPI_Issend(&buff, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        }
    }





    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) {
        tf = MPI_Wtime();
        printf("Total time: %.2f\n", tf-ts);
    }

    MPI_Finalize();
}
