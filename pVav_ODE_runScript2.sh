#!/bin/sh
#SBATCH -n 20

mpirun -np 20 ./experimentalData_ODE.exe "20" "InputFolder//outputTwo.txt" "1"