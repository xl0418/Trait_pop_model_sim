#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --partition=gelifes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --cpus-per-task=32
#SBATCH --job-name=perf32
#SBATCH --output=perf32

module load Python/3.6.4-foss-2018a
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python3 abcpp/perf.py >perf32.txt
