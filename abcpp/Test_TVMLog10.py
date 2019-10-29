import sys, os
sys.path.append('C:/Liang/abcpp_ms6/abcpp')
from Dvtraitsim_TVM import DVSimTVMLog10
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
from matplotlib.pylab import *
import numpy as np
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
theta = 3  # optimum of natural selection
gamma = 1.006e-07  # intensity of natural selection
r = 1  # growth rate
a = 2.346e-06  # intensity of competition
K = 10e12  # carrying capacity
nu=1.1e-05
Vmax = 1
scalar = 20000


#full tree and pruned tree directory
dir_path = './'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)

# parameter settings
obs_param = DVParamLiang(gamma=gamma, a=a, K=K,h=1, nu=nu, r=1, theta=theta,V00=.1,V01=.1, Vmax=Vmax, inittrait=theta, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)


population = 1000
obs_param_TVMlog10 = np.tile(obs_param, (population, 1))
simmodelTVM = dvcpp.DVSimTVMLog10(td, obs_param_TVMlog10)
valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0].size
print(valid_TVM);


for rep in range(100):
    simresult = DVSimTVMLog10(td,obs_param)
    if simresult['sim_time'] == td.sim_evo_time:
        break
    else:
        print('%d simulations are all junks! Try more!' % rep)

simresult['Z']

