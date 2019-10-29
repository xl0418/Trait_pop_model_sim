import numpy as np
from dvtraitsim_shared import DVTreeData
import dvtraitsim_cpp as dvcpp


def competition_functions_Liang(a, zi, nj):
	""" competition functions, Liang's model.

	returns beta = Sum_j( exp(-a(zi-zj)^2) * Nj)
			sigma = Sum_j( 2a * (zi-zj) * exp(-a(zi-zj)^2) * Nj)
			sigmaSqr = Sum_j( 4a^2 * (zi-zj)^2 * exp(-a(zi-zj)^2) * Nj)
	"""
	T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
	t1 = np.exp(-a * T ** 2) * nj
	t2 = (2 * a) * T
	beta = np.sum(t1, axis=1)
	sigma = np.sum(t2 * t1, axis=1)
	sigmasqr = np.sum(t2 ** 2 * t1, axis=1)
	return beta, sigma, sigmasqr


def DVSimTVP(td, param):
	# parameters from DVParamLiang
	gamma = param[0]
	a = param[1]
	K = param[2]
	h2 = param[3] * param[3]
	nu = param[4]
	r = param[5]
	theta = param[6]
	V00 = param[7]
	V01 = param[8]
	Vmax = param[9]
	inittrait = param[10]
	initpop = param[11]
	initpop_sigma = param[12]
	break_on_mu = bool(param[13])

	sim_evo_time = td.sim_evo_time
	events = td.sim_events

	# Initialize trait evolution and population evolution matrices
	trait_RI_dr = np.zeros((sim_evo_time + 1, td.total_species))  # trait
	population_RI_dr = np.zeros((sim_evo_time + 1, td.total_species)).astype(np.int32)  # population
	existing_species = td.traittable
	idx = np.where(existing_species[0] == 1)[0]    # existing species
	assert(idx.size == 2), "number of existing species shall be 2 at simulation start"

	# Initialize trait variances
	V = np.zeros((sim_evo_time + 1, td.total_species))
	V[0, idx] = [V00, V01];

	#  initialize condition for species trait and population
	trait_RI_dr[0, idx] = inittrait  # trait for species
	population_RI_dr[0, idx] = np.random.normal(initpop, initpop_sigma, 2).astype(np.int32)
	node = 0;
	next_event = events[node];

	# trait-population co-evolution model, Liang
	for i in range(sim_evo_time):
		# pull current state
		Ni = population_RI_dr[i, idx]
		Vi = V[i, idx]
		zi = trait_RI_dr[i, idx]
		Ki = K
		dtz = theta - zi
		beta, sigma, sigmasqr = competition_functions_Liang(a, zi, Ni)

		# update
		var_trait = Vi / (2.0 * Ni)
		trait_RI_dr[i + 1, idx] = zi + h2 * Vi * (2.0 * gamma * dtz + 1 / Ki * sigma) + np.random.normal(0.0, var_trait)
		mu = Ni * r * np.exp(-gamma * dtz**2 + (1 - beta / Ki)) # un-truncated mean
		if np.any(mu <= 1.0):       # mu < 1.0 + 1.11e-16
			if (break_on_mu):
				print(i, "invalid mean population size")
				break
		ztp_lambda = dvcpp.ztp_lambda_from_untruncated_mean(mu)
		population_RI_dr[i + 1, idx] = dvcpp.ztpoisson(ztp_lambda)
		V[i + 1, idx] = (1-h2/2.0)*Vi + 2.0*h2 * Ni * nu * Vmax / (1.0 + 4.0 * Ni * nu) \
						+ h2/2.0 * Vi**2 * (
							-2.0 * gamma + 4.0 * gamma**2 * dtz**2 +
								1.0 / Ki * (2.0 * a * beta - sigmasqr) + 4.0 * gamma / Ki *
								dtz * sigma + sigma**2 / Ki**2
							)

		# events
		while (i + 1) == next_event[0]:
			daughter = next_event[2]
			if (daughter == -1):
				# extinction
				extinct_species = next_event[1]
				V[i + 1, extinct_species] = None
				trait_RI_dr[i + 1, extinct_species] = None
				population_RI_dr[i + 1, extinct_species] = 0
			else:
				# speciation
				parent = next_event[1]
				parentN = population_RI_dr[i + 1, parent]
				if parentN <= 1:
					print(i, "attempt to split singleton")
					# results in split <- 0, will be trapped by sanity check below
				split = dvcpp.split_binomial50(parentN)
				population_RI_dr[i + 1, daughter] = parentN - split
				population_RI_dr[i + 1, parent] = split
				V[i + 1, parent] *= 0.5
				V[i + 1, daughter] = V[i + 1, parent]
				trait_RI_dr[i + 1, daughter] = trait_RI_dr[i + 1, parent]
			# advance to next event/node
			node = node + 1
			next_event = events[node];
			idx = np.where(existing_species[node] == 1)[0]

		# sanity check
		if np.any(population_RI_dr[i + 1, idx] < 1):
			print(i, 'Inconsistent extinction')
			break
		if np.any(V[i + 1, idx] < 0.0) or np.any(V[i + 1, idx] > 100000.0):
			print(i, 'runaway variance')
			break

	row_ext = np.where(population_RI_dr == 0)[0]
	col_ext = np.where(population_RI_dr == 0)[1]
	V[row_ext, col_ext] = None
	trait_RI_dr[row_ext, col_ext] = None
	return { 'sim_time': i + 1, 'N': population_RI_dr[-1], 'Z': trait_RI_dr[-1], 'V': V[-1] }

