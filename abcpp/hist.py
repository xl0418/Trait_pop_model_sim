import os
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path

dir_path = os.path.dirname(os.path.realpath(__file__))
files = dir_path + '/../tree_data/example1/'
td = DVTreeData(path=files, scalar=10000)

param = DVParam(gamma=0.001, a=1, K=10000, nu=0.000001, r=1, theta=0, Vmax=1, inittrait=0, initpop=500, initpop_sigma=10)
params = np.tile(param, (1000, 1))       # duplicate
obs = dvcpp.DVSim(td, params)
idx = np.where(obs['sim_time'] == td.sim_evo_time)[0]
print("found", idx.size, "valid simulations")
obsN = obs['N'][idx]
obsV = obs['V'][idx]
obsZ = obs['Z'][idx]

fig, ax = plt.subplots(3, 2)
bins = 50;
ax[0,0].hist(obsZ[:,0], bins, density=1)
ax[1,0].hist(obsZ[:,td.total_species // 2], bins, density=1)
ax[2,0].hist(obsZ[:,-1], bins, density=1)
ax[0,1].hist(obsV[:,0], bins, density=1)
ax[1,1].hist(obsV[:,td.total_species // 2], bins, density=1)
ax[2,1].hist(obsV[:,-1], bins, density=1)
fig.tight_layout()
plt.show()
