#!/bin/sh
#SBATCH -n 20

mpirun -np 5 ./pVav_ODE_GMM.exe "20" "InputFolder//outputOne.txt" "0"