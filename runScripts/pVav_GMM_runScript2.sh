#!/bin/sh
#SBATCH -n 20

mpirun -np 20 ./pVav_ODE_GMM.exe "20" "InputFolder//outputTwo.txt" "1"