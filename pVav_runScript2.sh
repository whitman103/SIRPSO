#!/bin/sh
#SBATCH -n 20

mpirun -np 20 ./pVav_MPI_means.exe "20" "InputFolder//outputTwo.txt" "1"