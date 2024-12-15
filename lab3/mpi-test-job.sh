#!/bin/bash
#SBATCH -p main
#SBATCH -n4
module load openmpi
mpicc -o main main.c
mpirun main

$ sbatch mpi-test-job
$ squeue -j JOBID
$ scancel JOBID
$ squeue
