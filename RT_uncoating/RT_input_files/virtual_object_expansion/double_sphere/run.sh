#!/bin/bash

#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -p <partition>
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=sp_3

# Load the MPI module(s) for your cluster, e.g.:
# module load <mpi-module>

nodes=1
ntasks=32
cores=$((nodes * ntasks))

# Path to the LAMMPS executable, built with the custom gauss/wall pair style (29Aug2024 used here)
lmp=/path/to/lammps/src/lmp_mpi

mpirun -np ${cores} ${lmp} -in inCapsid -var SEED ${RANDOM}
