#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --partition=gelifes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --cpus-per-task=32

module load Python/3.6.4-foss-2018a
srun python3 abcpp/evo.py
