#!/bin/bash
#SBATCH -A SNIC2018-5-16
#SBATCH -J merge
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 1-00:00:00
#SBATCH -n 1

# export module libraries
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/corgi/lib

# activate threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# go to working directory
cd /pfs/nobackup/home/n/natj/corgi/examples/loadbalance/

srun python merge_output.py
