import argparse
import sys
import os
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
from scipy.stats import norm
import csv
import dendropy
from dendropy.model import continuous
from itertools import repeat
from itertools import starmap
from pic_compute import pic_compute
from tp_update_theta import tp_update


# parse arguments
parser = argparse.ArgumentParser(description='BaleenWhale arguments')
parser.add_argument("--treedata", required=True, type=str, help="treedata directories")
parser.add_argument("--result", required=True, type=str, help="result npy file")
parser.add_argument("--num_threads", default=-1, required=False, type=int, help="number of threads")
args = parser.parse_args()
files = args.treedata
output = args.result
num_threads = args.num_threads

#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


K_TVP=1e6
K_TV = 1e6
K_TVM = 1e12
nu=1e-4


td = DVTreeData(path=files, scalar=20000)

with open(files + 'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(files + 'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:, 1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:, 0] == label)[0][0])

logTL = lengthdata_array[length_index, 1].astype(np.float)
length = logTL

# reorder the traits according to the simulation order
sim_species_label = ["Balaena_mysticetus", "Balaenoptera_acutorostrata", "Caperea_marginata",
                     "Balaenoptera_borealis",
                     "Balaenoptera_physalus", "Eschrichtius_robustus", "Balaenoptera_musculus",
                     "Balaenoptera_omurai",
                     "Eubalaena_australis", "Megaptera_novaeangliae", "Balaenoptera_bonaerensis",
                     "Balaenoptera_brydei",
                     "Balaenoptera_edeni", "Eubalaena_glacialis", "Eubalaena_japonica"]

obsZ_ordered_sim = length[
    [np.where(sim_species_label[i] == extantlabels_array)[0][0] for i in range(15)]]
obsZ = obsZ_ordered_sim
meantrait = np.mean(obsZ)
# PIC calculation
taxa1 = dendropy.TaxonNamespace()
dataset_combined = dendropy.DataSet.get(path=files + "bw_char.nex", schema="nexus")
tree_emp = dataset_combined.tree_lists[0][0]
chars_emp = dataset_combined.char_matrices[0]
pic_emp = continuous.PhylogeneticIndependentConstrasts(tree=tree_emp,
                                                       char_matrix=chars_emp)
ctree_emp = pic_emp.contrasts_tree(character_index=0,
                                   annotate_pic_statistics=True, state_values_as_node_labels=False,
                                   corrected_edge_lengths=False)
emp_pic = []
label = []
for nd in ctree_emp.postorder_internal_node_iter():
    emp_pic.append(nd.pic_contrast_standardized)
    label.append(int(nd.label))

emp_pic_orded_node = abs(np.array(emp_pic)[np.argsort(label)])

tree_sim = dendropy.Tree.get(
    path=files + "bw.nex", schema="nexus",
    taxon_namespace=taxa1)
print('trying to estimate the parameters', '...')

# let's try to find a true simulation:
sampleparam_TVP = DVParamLiang(gamma=1, a=1, K=K_TVP, h=1, nu=nu, r=1, theta=0, V00=.0001,
                               V01=.0001, Vmax=100, inittrait=meantrait, initpop=1e5,
                               initpop_sigma=10.0, break_on_mu=False, num_threads=num_threads)
sampleparam_TV = DVParamLiang(gamma=1, a=1, K=K_TV, h=1, nu=nu, r=1, theta=0, V00=.0001, V01=.0001,
                              Vmax=100, inittrait=meantrait, initpop=1e5,
                              initpop_sigma=10.0, break_on_mu=False, num_threads=num_threads)
sampleparam_TVM = DVParamLiang(gamma=1, a=1, K=K_TVM, h=1, nu=nu, r=1, theta=0, V00=.0001,
                               V01=.0001, Vmax=100.0, inittrait=meantrait, initpop=1e5,
                               initpop_sigma=10.0, break_on_mu=False, num_threads=num_threads)

# pop = dvcpp.DVSim(td, obs_param)

population = 40000
generations = 30
total_population = population * 3

lefttrait = np.min(obsZ)
righttrait = np.max(obsZ)

prior = [3, 5.1, 80, 270, 0.0, nu * 100, 0, 1e-3, lefttrait, righttrait]
prior_TVM = [0.5, 5, 800, 1500, 0.0, nu * 100, 0, 0.5e-3, lefttrait, righttrait]

params_TVP = np.tile(sampleparam_TVP, (population, 1))  # duplicate
params_TVP[:, 0] = np.random.uniform(prior[0], prior[1], params_TVP.shape[0])  # randomize 'gamma'
params_TVP[:, 1] = np.random.uniform(prior[2], prior[3], params_TVP.shape[0])  # randomize 'a'
params_TVP[:, 4] = np.random.uniform(prior[4], prior[5], params_TVP.shape[0])  # randomize 'nu'
params_TVP[:, 6] = np.random.uniform(prior[8], prior[9], params_TVP.shape[0])  # randomize 'theta'
params_TVP[:, 9] = np.random.uniform(prior[6], prior[7], params_TVP.shape[0])  # randomize 'Vm'

params_TV = np.tile(sampleparam_TV, (population, 1))  # duplicate
params_TV[:, 0] = np.random.uniform(prior[0], prior[1], params_TV.shape[0])  # randomize 'gamma'
params_TV[:, 1] = np.random.uniform(prior[2], prior[3], params_TV.shape[0])  # randomize 'a'
params_TV[:, 4] = np.random.uniform(prior[4], prior[5], params_TV.shape[0])  # randomize 'nu'
params_TV[:, 6] = np.random.uniform(prior[8], prior[9], params_TV.shape[0])  # randomize 'theta'
params_TV[:, 9] = np.random.uniform(prior[6], prior[7], params_TV.shape[0])  # randomize 'Vm'

params_TVM = np.tile(sampleparam_TVM, (population, 1))  # duplicate
params_TVM[:, 0] = np.random.uniform(prior_TVM[0], prior_TVM[1],
                                     params_TVM.shape[0])  # randomize 'gamma'
params_TVM[:, 1] = np.random.uniform(prior_TVM[2], prior_TVM[3],
                                     params_TVM.shape[0])  # randomize 'a'
params_TVM[:, 4] = np.random.uniform(prior_TVM[4], prior_TVM[5],
                                     params_TVM.shape[0])  # randomize 'nu'
params_TVM[:, 6] = np.random.uniform(prior_TVM[8], prior_TVM[9],
                                     params_TVM.shape[0])  # randomize 'theta'
params_TVM[:, 9] = np.random.uniform(prior_TVM[6], prior_TVM[7],
                                     params_TVM.shape[0])  # randomize 'Vm'

# model choice
model_index = np.array([0, 1, 2])
model_params = np.repeat(model_index, repeats=population)
model_data = np.zeros(shape=(generations + 1, total_population))
model_data[0, :] = model_params
propose_model = model_params

# store parameters used
# TVP
gamma_data_TVP = np.zeros(shape=(generations + 1, population))
a_data_TVP = np.zeros(shape=(generations + 1, population))
nu_data_TVP = np.zeros(shape=(generations + 1, population))
vm_data_TVP = np.zeros(shape=(generations + 1, population))
theta_data_TVP = np.zeros(shape=(generations + 1, population))

# TV
gamma_data_TV = np.zeros(shape=(generations + 1, population))
a_data_TV = np.zeros(shape=(generations + 1, population))
nu_data_TV = np.zeros(shape=(generations + 1, population))
vm_data_TV = np.zeros(shape=(generations + 1, population))
theta_data_TV = np.zeros(shape=(generations + 1, population))

# TVM
gamma_data_TVM = np.zeros(shape=(generations + 1, population))
a_data_TVM = np.zeros(shape=(generations + 1, population))
nu_data_TVM = np.zeros(shape=(generations + 1, population))
vm_data_TVM = np.zeros(shape=(generations + 1, population))
theta_data_TVM = np.zeros(shape=(generations + 1, population))

gamma_data_TVP[0, :] = params_TVP[:, 0]
a_data_TVP[0, :] = params_TVP[:, 1]
nu_data_TVP[0, :] = params_TVP[:, 4]
vm_data_TVP[0, :] = params_TVP[:, 9]
theta_data_TVP[0, :] = params_TVP[:, 6]

gamma_data_TV[0, :] = params_TV[:, 0]
a_data_TV[0, :] = params_TV[:, 1]
nu_data_TV[0, :] = params_TV[:, 4]
vm_data_TV[0, :] = params_TV[:, 9]
theta_data_TV[0, :] = params_TV[:, 6]

gamma_data_TVM[0, :] = params_TVM[:, 0]
a_data_TVM[0, :] = params_TVM[:, 1]
nu_data_TVM[0, :] = params_TVM[:, 4]
vm_data_TVM[0, :] = params_TVM[:, 9]
theta_data_TVM[0, :] = params_TVM[:, 6]

fitness = np.zeros(shape=(generations, total_population))
# Initialize the weights.
weight_model = np.zeros(total_population)
weight_model.fill(1 / total_population)

# weights for paras of TVP
weight_gamma_TVP = np.zeros(population)
weight_gamma_TVP.fill(1 / population)
weight_a_TVP = np.zeros(population)
weight_a_TVP.fill(1 / population)
weight_nu_TVP = np.zeros(population)
weight_nu_TVP.fill(1 / population)
weight_vm_TVP = np.zeros(population)
weight_vm_TVP.fill(1 / population)
weight_theta_TVP = np.zeros(population)
weight_theta_TVP.fill(1 / population)

# weights for paras of TV
weight_gamma_TV = np.zeros(population)
weight_gamma_TV.fill(1 / population)
weight_a_TV = np.zeros(population)
weight_a_TV.fill(1 / population)
weight_nu_TV = np.zeros(population)
weight_nu_TV.fill(1 / population)
weight_vm_TV = np.zeros(population)
weight_vm_TV.fill(1 / population)
weight_theta_TV = np.zeros(population)
weight_theta_TV.fill(1 / population)

# weights for paras of TVM
weight_gamma_TVM = np.zeros(population)
weight_gamma_TVM.fill(1 / population)
weight_a_TVM = np.zeros(population)
weight_a_TVM.fill(1 / population)
weight_nu_TVM = np.zeros(population)
weight_nu_TVM.fill(1 / population)
weight_vm_TVM = np.zeros(population)
weight_vm_TVM.fill(1 / population)
weight_theta_TVM = np.zeros(population)
weight_theta_TVM.fill(1 / population)

for g in range(generations):
    # model 0
    TVP_sample_length = len(np.where(model_data[g, :] == 0)[0])
    TV_sample_length = len(np.where(model_data[g, :] == 1)[0])
    TVM_sample_length = len(np.where(model_data[g, :] == 2)[0])
    Z = np.zeros((1, td.total_species))
    pics = np.zeros((1, td.total_species - 1))

    if TVP_sample_length > 0:
        print('TVP simulations start...')

        simmodelTVP = dvcpp.DVSimTVP(td, params_TVP)
        valid_TVP = np.where(simmodelTVP['sim_time'] == td.sim_evo_time)[0]
        num_valid_sims_TVP = len(valid_TVP)
        if num_valid_sims_TVP < 20:
            print("WARNING:Valid simulations are too scarce!")
        if num_valid_sims_TVP > 0:
            Z_modelTVP = simmodelTVP['Z'][valid_TVP]
            pic_ordered_list = list(starmap(pic_compute,
                                                 zip(repeat(tree_sim), Z_modelTVP, repeat(taxa1),
                                                     range(num_valid_sims_TVP))))
            # in case the parallel computation returns disordered output
            order_list = []
            contrast_list = []
            for i in range(num_valid_sims_TVP):
                order_list.append(pic_ordered_list[i][1])
                contrast_list.append(pic_ordered_list[i][0])
            ordered_contrast_list = [contrast_list[item] for item in np.argsort(order_list)]
            contrast_array = abs(np.vstack(ordered_contrast_list))
            # i, j = argsort2D(Z_modelTVP)
            # Z_modelTVP = Z_modelTVP[i, j]
            # V = pop['V'][valid][i, j]
            Z_modelTVP = np.nan_to_num(Z_modelTVP)
            Z = np.vstack([Z, Z_modelTVP])
            pics = np.vstack([pics, contrast_array])
            # V = np.nan_to_num(V)
            # GOF: Goodness of fit
        if len(valid_TVP) == 0:
            print('No complete results from TVP model ')

    if TV_sample_length > 0:
        print('TV simulations start...')
        # model 1
        # for param_drury in params_DR:
        simmodelTV = dvcpp.DVSimTV(td, params_TV)
        valid_TV = np.where(simmodelTV['sim_time'] == td.sim_evo_time)[0]
        num_valid_sims_TV = len(valid_TV)
        Z_modelTV = simmodelTV['Z'][valid_TV]
        pic_ordered_list = list(starmap(pic_compute,
                                             zip(repeat(tree_sim), Z_modelTV, repeat(taxa1),
                                                 range(num_valid_sims_TV))))
        # in case the parallel computation returns disordered output
        order_list = []
        contrast_list = []
        for i in range(num_valid_sims_TV):
            order_list.append(pic_ordered_list[i][1])
            contrast_list.append(pic_ordered_list[i][0])
        ordered_contrast_list = [contrast_list[item] for item in np.argsort(order_list)]
        contrast_array = abs(np.vstack(ordered_contrast_list))
        # i, j = argsort2D(Z_modelTV)
        # Z_modelTV = Z_modelTV[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTV = np.nan_to_num(Z_modelTV)
        Z = np.vstack([Z, Z_modelTV])
        pics = np.vstack([pics, contrast_array])

    if TVM_sample_length > 0:
        print('TVM simulations start...')
        simmodelTVM = dvcpp.DVSimTVMLog10(td, params_TVM)
        valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0]
        num_valid_sims_TVM = len(valid_TVM)
        Z_modelTVM = simmodelTVM['Z'][valid_TVM]
        pic_ordered_list = list(starmap(pic_compute,
                                             zip(repeat(tree_sim), Z_modelTVM, repeat(taxa1),
                                                 range(num_valid_sims_TVM))))
        # in case the parallel computation returns disordered output
        order_list = []
        contrast_list = []
        for i in range(num_valid_sims_TVM):
            order_list.append(pic_ordered_list[i][1])
            contrast_list.append(pic_ordered_list[i][0])
        ordered_contrast_list = [contrast_list[item] for item in np.argsort(order_list)]
        contrast_array = abs(np.vstack(ordered_contrast_list))
        # i, j = argsort2D(Z_modelTVM)
        # Z_modelTVM = Z_modelTVM[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTVM = np.nan_to_num(Z_modelTVM)
        Z = np.vstack([Z, Z_modelTVM])
        pics = np.vstack([pics, contrast_array])

    Z = Z[1:, ]
    pics = pics[1:, ]

    valid = np.concatenate([valid_TVP, np.array(valid_TV) + TVP_sample_length,
                            np.array(valid_TVM) + TVP_sample_length + TV_sample_length
                            ]).astype(int)

    eudis = normalized_norm(Z, obsZ)
    eudis_pic = normalized_norm(pics, emp_pic_orded_node)
    # eudis = eudistance(Z, obsZ)

    fitness[g, valid] += 1 - eudis
    fitness[g, valid] += 1 - eudis_pic

    # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g, :])[-int(total_population // 4)]  # best 25%
    fit_index = np.where(fitness[g, :] > fitness[g, q5])[0]

    modelTVPperc = len(np.where(propose_model[fit_index] == 0)[0]) / len(fit_index)
    modelTVperc = len(np.where(propose_model[fit_index] == 1)[0]) / len(fit_index)
    modelTVMperc = len(np.where(propose_model[fit_index] == 2)[0]) / len(fit_index)

    print('Iteration = %d 25th Model TVP: %.1f%% ;  Model TV: %.1f%% ; Model TVM: %.1f%%...'
          % (g, modelTVPperc * 100, modelTVperc * 100, modelTVMperc * 100))
    print('Average fitness: %f' % np.mean(fitness[g, fit_index]))
    # reevaluate the weight of the best fitted  models
    weight_model_bestfitted = weight_model[fit_index] * fitness[g, fit_index] / sum(
        weight_model[fit_index] * fitness[g, fit_index])
    # sample new models from the fitness of previous best fitted models
    sample_model_index = sorted(
        np.random.choice(fit_index, total_population, p=weight_model_bestfitted))

    propose_model = model_params

    q5_TVP = np.argsort(fitness[g, :population])[-int(population // 200)]  # best 5%
    q5_TV = np.argsort(fitness[g, population:2 * population])[
                -int(population // 200)] + population  # best 5%
    q5_TVM = np.argsort(fitness[g, 2 * population:3 * population])[
                 -int(population // 200)] + 2 * population  # best 5%

    fit_index_TVP = np.where(fitness[g, :population] > fitness[g, q5_TVP])[0]
    fit_index_TV = np.where(fitness[g, population:2 * population] > fitness[g, q5_TV])[
                       0] + population
    fit_index_TVM = np.where(fitness[g, 2 * population:] > fitness[g, q5_TVM])[0] + 2 * population

    previous_bestfitted_index_TVP = fit_index_TVP
    previous_bestfitted_index_TV = fit_index_TV - population
    previous_bestfitted_index_TVM = fit_index_TVM - 2 * population

    chosengamma_TVP, chosena_TVP, chosennu_TVP, chosenvm_TVP, chosentheta_TVP = np.mean(
        params_TVP[previous_bestfitted_index_TVP, 0]), \
                                                                                np.mean(params_TVP[
                                                                                            previous_bestfitted_index_TVP, 1]), \
                                                                                np.mean(params_TVP[
                                                                                            previous_bestfitted_index_TVP, 4]), \
                                                                                np.mean(params_TVP[
                                                                                            previous_bestfitted_index_TVP, 9]), \
                                                                                np.mean(params_TVP[
                                                                                            previous_bestfitted_index_TVP, 6])

    chosengamma_TV, chosena_TV, chosennu_TV, chosenvm_TV, chosentheta_TV = np.mean(
        params_TV[previous_bestfitted_index_TV, 0]), \
                                                                           np.mean(params_TV[
                                                                                       previous_bestfitted_index_TV, 1]), \
                                                                           np.mean(params_TV[
                                                                                       previous_bestfitted_index_TV, 4]), \
                                                                           np.mean(params_TV[
                                                                                       previous_bestfitted_index_TV, 9]), \
                                                                           np.mean(params_TV[
                                                                                       previous_bestfitted_index_TV, 6])

    chosengamma_TVM, chosena_TVM, chosennu_TVM, chosenvm_TVM, chosentheta_TVM = np.mean(
        params_TVM[previous_bestfitted_index_TVM, 0]), \
                                                                                np.mean(params_TVM[
                                                                                            previous_bestfitted_index_TVM, 1]), \
                                                                                np.mean(params_TVM[
                                                                                            previous_bestfitted_index_TVM, 4]), \
                                                                                np.mean(params_TVM[
                                                                                            previous_bestfitted_index_TVM, 9]), \
                                                                                np.mean(params_TVM[
                                                                                            previous_bestfitted_index_TVM, 6])

    print('Mean estimates: TVP gamma: %.3e ; a: %.3e ; nu: %.3e ; Vm : %f; theta : %f' % (
    chosengamma_TVP, chosena_TVP, chosennu_TVP, chosenvm_TVP, chosentheta_TVP))
    print('Mean estimates: TV gamma: %.3e ; a: %.3e ; nu: %.3e ; Vm : %f; theta : %f' % (
    chosengamma_TV, chosena_TV, chosennu_TV, chosenvm_TV, chosentheta_TV))
    print('Mean estimates: TVM gamma: %.3e ; a: %.3e ; nu: %.3e ; Vm : %f; theta : %f' % (
    chosengamma_TVM, chosena_TVM, chosennu_TVM, chosenvm_TVM, chosentheta_TVM))

    model_data[g + 1, :] = propose_model
    gamma_data_TVP[g + 1, :] = params_TVP[:, 0]
    a_data_TVP[g + 1, :] = params_TVP[:, 1]
    nu_data_TVP[g + 1, :] = params_TVP[:, 4]
    vm_data_TVP[g + 1, :] = params_TVP[:, 9]
    theta_data_TVP[g + 1, :] = params_TVP[:, 6]

    gamma_data_TV[g + 1, :] = params_TV[:, 0]
    a_data_TV[g + 1, :] = params_TV[:, 1]
    nu_data_TV[g + 1, :] = params_TV[:, 4]
    vm_data_TV[g + 1, :] = params_TV[:, 9]
    theta_data_TV[g + 1, :] = params_TV[:, 6]

    gamma_data_TVM[g + 1, :] = params_TVM[:, 0]
    a_data_TVM[g + 1, :] = params_TVM[:, 1]
    nu_data_TVM[g + 1, :] = params_TVM[:, 4]
    vm_data_TVM[g + 1, :] = params_TVM[:, 9]
    theta_data_TVM[g + 1, :] = params_TVM[:, 6]

    if len(np.where(propose_model == 0)[0]) > 0:
        params_TVP_update = params_TVP[:, [0, 1, 4, 9, 6]]
        modelinex = 0
        # update TP paras and weights
        weight_gamma_TVP, weight_a_TVP, weight_nu_TVP, weight_vm_TVP, weight_theta_TVP, \
        propose_gamma_TVP, propose_a_TVP, propose_nu_TVP, propose_vm_TVP, propose_theta_TVP = \
            tp_update(previous_bestfitted_index_TVP, propose_model, params_TVP_update,
                      weight_gamma_TVP,
                      weight_a_TVP, weight_nu_TVP, weight_vm_TVP, weight_theta_TVP, modelinex)
        modelTVP = np.where(propose_model == modelinex)
        params_TVP = np.tile(sampleparam_TVP, (len(modelTVP[0]), 1))
        params_TVP[:, 0] = propose_gamma_TVP
        params_TVP[:, 1] = propose_a_TVP
        params_TVP[:, 4] = propose_nu_TVP
        params_TVP[:, 9] = propose_vm_TVP
        params_TVP[:, 6] = propose_theta_TVP

    if len(np.where(propose_model == 1)[0]) > 0:
        params_TV_update = params_TV[:, [0, 1, 4, 9, 6]]
        modelinex = 1
        if len(valid_TV) > 0:
            # update TP paras and weights
            weight_gamma_TV, weight_a_TV, weight_nu_TV, weight_vm_TV, weight_theta_TV, \
            propose_gamma_TV, propose_a_TV, propose_nu_TV, propose_vm_TV, propose_theta_TV = \
                tp_update(previous_bestfitted_index_TV, propose_model, params_TV_update,
                          weight_gamma_TV,
                          weight_a_TV, weight_nu_TV, weight_vm_TV, weight_theta_TV, modelinex)
            modelTV = np.where(propose_model == modelinex)
            params_TV = np.tile(sampleparam_TV, (len(modelTV[0]), 1))
            params_TV[:, 0] = propose_gamma_TV
            params_TV[:, 1] = propose_a_TV
            params_TV[:, 4] = propose_nu_TV
            params_TV[:, 9] = propose_vm_TV
            params_TV[:, 6] = propose_theta_TV

    if len(np.where(propose_model == 2)[0]) > 0:
        params_TVM_update = params_TVM[:, [0, 1, 4, 9, 6]]
        modelinex = 2
        # update TP paras and weights
        weight_gamma_TVM, weight_a_TVM, weight_nu_TVM, weight_vm_TVM, weight_theta_TVM, \
        propose_gamma_TVM, propose_a_TVM, propose_nu_TVM, propose_vm_TVM, propose_theta_TVM = \
            tp_update(previous_bestfitted_index_TVM, propose_model, params_TVM_update,
                      weight_gamma_TVM,
                      weight_a_TVM, weight_nu_TVM, weight_vm_TVM, weight_theta_TVM, modelinex)
        modelTVM = np.where(propose_model == modelinex)
        params_TVM = np.tile(sampleparam_TVM, (len(modelTVM[0]), 1))
        params_TVM[:, 0] = propose_gamma_TVM
        params_TVM[:, 1] = propose_a_TVM
        params_TVM[:, 4] = propose_nu_TVM
        params_TVM[:, 9] = propose_vm_TVM
        params_TVM[:, 6] = propose_theta_TVM

    #

    para_data = {'model_data': model_data, 'fitness': fitness, 'Z': Z,
                 'gamma_data_TVP': gamma_data_TVP, 'a_data_TVP': a_data_TVP,
                 'nu_data_TVP': nu_data_TVP,
                 'vm_data_TVP': vm_data_TVP, 'theta_data_TVP': theta_data_TVP,
                 'gamma_data_TV': gamma_data_TV, 'a_data_TV': a_data_TV, 'nu_data_TV': nu_data_TV,
                 'vm_data_TV': vm_data_TV, 'theta_data_TV': theta_data_TV,
                 'gamma_data_TVM': gamma_data_TVM, 'a_data_TVM': a_data_TVM,
                 'nu_data_TVM': nu_data_TVM,
                 'vm_data_TVM': vm_data_TVM, 'theta_data_TVM': theta_data_TVM
                 }

    np.save(output, para_data)
