#!/bin/bash

#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node=64
#SBATCH --job-name=cone_2

module unload intel-oneapi
module load intel-mpi

nodes=1
ntasks=64
cores=$((nodes * ntasks))

mpirun -np ${cores} /ocean/projects/mca94p017p/kuntalg/softwares/LAMMPS/lammps-29Aug2024/src/lmp_mpi -in inCapsid -var SEED ${RANDOM}
