import os
from time import perf_counter
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang, DVParamDrury
import dvtraitsim_cpp as dvcpp
import dvtraitsim_py as dvpy

dir_path = os.path.dirname(os.path.realpath(__file__))
files = dir_path + '/../tree_data/empirical_test2/full_data/'
td = DVTreeData(path=files, scalar=1000)


for r in range(0,1):
	K=100000000
	nu=1.0/(100.0*K)
	V0=(1/td.total_species)
	param = DVParamLiang(gamma=0.1, a=0.1, K=K, h=1, nu=nu, r=1.0, theta=0.0, V00=V0, V01=V0, Vmax=1.0, inittrait=0.0, initpop=500, initpop_sigma=10.0, break_on_mu=False)
	params = np.tile(param, (100, 1))       # duplicate

	print('running', params.shape[0], 'simulations (parallel):')
	start = perf_counter()
	dispatch = dvcpp.DVSimTVMLog10(td, params)
	time = perf_counter() - start
	print("valid simulations:", np.where(dispatch['sim_time'] == td.sim_evo_time)[0].size)
	print("in", time, "s")

	print('running', 1, 'simulation (python):')
	start = perf_counter()
	for p in range(0,1):
		R = dvpy.DVSimTVMLog10(td, param)

	time = perf_counter() - start
	print("in", time, "s")
