# The simulation code for the trait-population coevolution model

## Running on Windows

Open the terminal and check if python is installed correctly.

```bash
Trait_pop_model_sim>python --version
Python 3.7.0
```

Python 3.7.0 is expected as the output.

Then, go into the folder abcpp and run 

````bash
Trait_pop_model_sim\abcpp>python BaleenWhale_MS_pics.py --treedata treedata/ --result output_file_name --num_threads 6
````

The first argument indicates the input data of the phylogenetic tree. The second argument gives the output name. The format of the output is set to be npy files. The third argument determines how many threads to use. The default setting will use all available threads. 

## Running on the cluster (Linux)

### Upload the scripts

Upload the whole repositary to the cluster.

```bash
scp -r Trait_pop_model_sim p-user@peregrine.hpc.rug.nl:/home/p-user
```



### Install on the cluster

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


## Questions

- If you encountered this error message:

```bash
./build_pg_gelifes.sh: /bin/bash^M: bad interpreter: No such file or directory
```

you need to convert the script to linux coding by doing 

```bash
$ dos2unix build_pg_gelifes.sh
```
