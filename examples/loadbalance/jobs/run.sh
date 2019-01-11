#!/bin/bash
#SBATCH -A SNIC2018-5-16
#SBATCH -J lbal
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 0-10:00:00
#SBATCH -n 100

# export module libraries
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/corgi/lib

# activate threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# go to working directory
cd /pfs/nobackup/home/n/natj/corgi/examples/loadbalance/

mpirun python sim.py
