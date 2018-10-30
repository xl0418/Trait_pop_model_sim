import os
import sys
from time import perf_counter
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp


def normalize(v):
    norm = np.linalg.norm(v, ord=1, axis=0)
    return v/norm


# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_abserr(x, y):
    abs_err = np.sum(np.abs(normalize(x) - y), axis=1)
    max_err = np.nanmax(abs_err)
    return abs_err / max_err


# prep data
dir_path = os.path.dirname(os.path.realpath(__file__))
files = dir_path + '/../tree_data/example1/'
td = DVTreeData(path=files, scalar=1000)


# let's try to find a true simulation:
obs_param = DVParam(gamma=0.01, a=0.9, K=10000, nu=0.0001, r=1, theta=0, Vmax=1, inittrait=0, initpop=500, split_stddev=0.2)
print('try to find a completed true simulation with gamma =', obs_param[0], 'and a =', obs_param[1], '...')
obs = dvcpp.DVSim(td, np.tile(obs_param, (1000,1)))
valid_obs = np.where(obs['sim_time'] == td.evo_time)[0]
if valid_obs.size < 5:
    print('hopeless, does not compute.')
    sys.exit(-1)
print("found", valid_obs.size, "valid simulations")
i, j = argsort2D(obs['Z'][valid_obs])
obsN = normalize(np.mean(obs['N'][valid_obs][i,j], axis=0))
obsZ = normalize(np.mean(obs['Z'][valid_obs][i,j], axis=0))
obsV = normalize(np.mean(obs['V'][valid_obs][i,j], axis=0))

start = perf_counter()
population = 1000
generations = 10
params = np.tile(obs_param, (population, 1))                    # duplicate
params[:,0] = 0.01 #np.random.uniform(0.0, 1.0, params.shape[0])      # randomize 'gamma'
params[:,1] = 0.9 #np.random.uniform(0.0, 1.0, params.shape[0])      # randomize 'a'
for g in range(generations):
    pop = dvcpp.DVSim(td, params)
    
    # access fitness
    fitness = np.zeros(population)
    valid = np.where(pop['sim_time'] == td.evo_time)[0]
    if valid.size > 0:
        Z = pop['Z'][valid]
        i, j = argsort2D(Z)
        Z = Z[i, j]
        V = pop['V'][valid][i, j]
        fitness_z = 1.0 - normalized_abserr(Z, obsZ)
        fitness_v = 1.0 - normalized_abserr(V, obsV)
        fitness[valid] = 0.5 * (fitness_z + fitness_v)
    rank = np.argsort(fitness)
    q5 = rank[-population // 20:]     # best 5%

    # print something...
    print('Generation = %d  valid = %d  gamma = %f  a = %f  fitness = %f' % \
        (g, valid.size, np.mean(params[q5,0]), np.mean(params[q5,1]), np.mean(fitness[q5])))

    # reproduce
    offspring = dvcpp.discrete_distribution(fitness, population)
    params = params[offspring];

    # mutate
    params[:,0] = np.clip(params[:,0] + dvcpp.cauchy_distribution(b=0.01, num=population), 0.0, 1.0)
    params[:,1] = np.clip(params[:,1] + dvcpp.cauchy_distribution(b=0.01, num=population), 0.0, 1.0)
print("finished in", perf_counter() - start, "sec.")