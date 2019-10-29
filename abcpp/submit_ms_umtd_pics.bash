#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --partition=gelifes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=20
#SBATCH --job-name=MSumtd

python3 BaleenWhale_MS_umtd_pics.py $*
