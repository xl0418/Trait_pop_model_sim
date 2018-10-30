# The simulation code for the trait-population coevolution model

## Install

Follow the instruction below on cluster:

```bash
$ module purge
$ module load Python/3.6.4-foss-2018a
$ module load tbb
$ cd abcpp/
$ abcpp chmod +x ./build_pg_gelifes.sh
$ abcpp ./build_pg_gelifes.sh
```
Now, the compiling is done and you can submit the job to the cluster.

The script loopjob_pro2_nv.bash can submit multiple jobs to the cluster. 

```bash
$ sbatch loopjob_pro2_nv.bash
```

Here we want to test the model under 36 parameter combinations for one tree.
All the tree data is given in the folder *tree_data*. Change the tree number in the script evo_loop_nv.py to test different trees. 

```python
dir_path = '~'
files = dir_path + 'tree_data/example17/'
```


## Questions

- If you encountered this error message:

```bash
./build_pg_gelifes.sh: /bin/bash^M: bad interpreter: No such file or directory
```

you need to convert the script to linux coding by doing 

```bash
$ dos2unix build_pg_gelifes.sh
```
