import os
import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_master/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
from scipy.stats import norm

gamma = 0.001
a = 0.1
#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


# prep data
# dir_path = os.path.dirname(os.path.realpath(__file__))
# files = dir_path + '/../tree_data/example1/'
dir_path = 'c:/Liang/Googlebox/Research/Project2'
files = dir_path + '/treesim_newexp/example3/'

td = DVTreeData(path=files, scalar=1000)

prior = [0.5, 0.5, 0.5, 0.5,1e-11,1e-11]
gamma_prior_mean = prior[0]
gamma_prior_var = prior[1]
a_prior_mean = prior[2]
a_prior_var = prior[3]
nu_prior_mean = prior[4]
nu_prior_var = prior[5]
K=10e8
nu=1/(100*K)
# let's try to find a true simulation:
obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
             initpop_sigma = 10.0)
print('try to find a completed true simulation with gamma =', obs_param[0], 'and a =', obs_param[1], 'and nu =', obs_param[3],'...')
for r in range(1000):
    print(r)
    obs = dvcpp.DVSim(td, obs_param)
    pos_1pop = np.where(obs['N'] == 1)[0]
    print(obs['V'][pos_1pop])
    if obs['sim_time'] == td.sim_evo_time:
        break
if obs['sim_time'] < td.sim_evo_time:
    print('hopeless, does not compute.')
    sys.exit(-1)
s = np.argsort(obs['Z'])
obsN = obs['N'][s]
obsZ = obs['Z'][s]
obsV = obs['V'][s]
obsN = np.nan_to_num(obsN)
obsZ = np.nan_to_num(obsZ)
obsV = np.nan_to_num(obsV)
