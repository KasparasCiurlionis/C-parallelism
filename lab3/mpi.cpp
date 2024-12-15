#include <stdio.h>
#include <mpi.h>

int get_x(int x);

int main(int argc, char *argv[]) {
    for(int i=0;i<15;i++){
    printf("%d %d\n", i, get_x(i));
    }
}

int get_x(int x){
int sum=0;
int i=0;
while(x>sum){
    i++;
    sum+=i+1;
}
return i;
}
