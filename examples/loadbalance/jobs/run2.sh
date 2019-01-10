#!/bin/bash
#SBATCH -A SNIC2018-2-41
#SBATCH -J lbal4
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 0-10:00:00
#SBATCH -n 1000

# export module libraries
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/corgi/lib

# activate threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# go to working directory
cd /pfs/nobackup/home/n/natj/corgi/examples/loadbalance/

mpirun python sim.py
