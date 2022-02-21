#!/bin/bash -l
#SBATCH -J ajo_job
#SBATCH -o ajo%J.out
#SBATCH -e ajo%J.err
#SBATCH -t 12:00:00
#SBATCH -N 3
#SBATCH -p small
#SBATCH --ntasks-per-node=24
#number of cores reserved
(( ncores = SLURM_NNODES * 24 ))
echo "Running namd with $SLURM_NNODES nodes containing total of $ncores cores"
module load intel/14.0.4.211
module load cray-mpich/7.2.4

aprun -n $ncores   /omapolku/paradisgeneral/bin/paradis /omapolku/$1

