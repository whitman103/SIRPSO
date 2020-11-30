#!/bin/sh

load module GCC/9.3.0
load module GCC/9.3.0 OpenMPI/4.0.3 
load module GCC/9.3.0 OpenMPI/4.0.3 Boost/1.72.0
sbatch pVav_ODE_runScript1.sh
sbatch pVav_ODE_runScript2.sh
sbatch pVav_ODE_runScript3.sh

