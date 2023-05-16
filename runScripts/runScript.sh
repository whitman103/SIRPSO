#!/bin/sh
#SBATCH --ntasks=3
#SBATCH --nodes=3

mpirun -np 20 ./SIRPSO_MPI.exe "20" "InputFolder//outputOne.txt" "0" &
mpirun -np 20 ./SIRPSO_MPI.exe "20" "InputFolder//outputTwo.txt" "1" &
mpirun -np 20 ./SIRPSO_MPI.exe "20" "InputFolder//outputThree.txt" "2" &

wait