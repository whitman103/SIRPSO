#!/bin/sh
#SBATCH -n 20

mpirun -np 5 ./pVav_ODE_means.exe "5" "InputFolder//outputOne.txt" "0"