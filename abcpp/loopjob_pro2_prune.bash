#!/bin/bash
#SBATCH --partition=gelifes

#gamma=(0 0.001 0.01 0.1 0.5 1)
#a=(0 0.001 0.01 0.1 0.5 1)

for((s = 0; s <= 5; s++))
do
for((m = 0; m <= 5; m++))
do
echo "#!/bin/bash" > zMLjob$s$m
echo "#SBATCH --time=2-23:59:00" >> zMLjob$s$m
echo "#SBATCH --partition=gelifes" >> zMLjob$s$m
echo "#SBATCH --ntasks=1" >> zMLjob$s$m
echo "#SBATCH --nodes=1" >> zMLjob$s$m
echo "#SBATCH --cpus-per-task=32" >> zMLjob$s$m
echo "export OMP_NUM_THREADS="'$SLURM_CPUS_PER_TASK'  >> zMLjob$s$m
echo "module load Python/3.6.4-foss-2018a" >> zMLjob$s$m
echo "srun python3 evo_loop_prune.py "$s" "$m"" >> zMLjob$s$m
echo "rm zMLjob$s$m" >> zMLjob$s$m

sbatch --partition=gelifes --mem-per-cpu=1GB --job-name=smc$s$m --output=smc$s$m.log --mail-type=FAIL,TIME_LIMIT --mail-user=xl0418@gmail.com zMLjob$s$m

done
done
