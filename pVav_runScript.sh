#!/bin/sh
#SBATCH --ntasks=3
#SBATCH -n 60

mpirun -np 20 ./pVav_MPI_means.exe "20" "InputFolder//outputOne.txt" "0" &
mpirun -np 20 ./pVav_MPI_means.exe "20" "InputFolder//outputTwo.txt" "1" &
mpirun -np 20 ./pVav_MPI_means.exe "20" "InputFolder//outputThree.txt" "2" &

wait