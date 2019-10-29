import numpy as np
from scipy.stats import norm
def tp_update(previous_bestfitted_index,propose_model,params,weight_gamma,weight_a,
              weight_nu,weight_vm,weight_theta,modelindex):
    # Room in TP model to sample new paras
    previous_gamma = params[previous_bestfitted_index, 0]
    previous_a = params[previous_bestfitted_index, 1]
    previous_nu = params[previous_bestfitted_index, 2]
    previous_vm = params[previous_bestfitted_index, 3]
    previous_theta = params[previous_bestfitted_index, 4]

    weight_gamma = weight_gamma[previous_bestfitted_index] / sum(weight_gamma[previous_bestfitted_index])
    weight_a = weight_a[previous_bestfitted_index] / sum(weight_a[previous_bestfitted_index])
    weight_nu = weight_nu[previous_bestfitted_index] / sum(weight_nu[previous_bestfitted_index])
    weight_vm = weight_vm[previous_bestfitted_index] / sum(weight_vm[previous_bestfitted_index])
    weight_theta = weight_theta[previous_bestfitted_index] / sum(weight_theta[previous_bestfitted_index])

    gamma_pre_mean = np.sum(previous_gamma * weight_gamma)
    gamma_pre_var = np.sum((previous_gamma - gamma_pre_mean) ** 2 * weight_gamma)
    a_pre_mean = np.sum(previous_a * weight_a)
    a_pre_var = np.sum((previous_a - a_pre_mean) ** 2 * weight_a)
    nu_pre_mean = np.sum(previous_nu * weight_nu)
    nu_pre_var = np.sum((previous_nu - nu_pre_mean) ** 2 * weight_nu)
    vm_pre_mean = np.sum(previous_vm * weight_vm)
    vm_pre_var = np.sum((previous_vm - vm_pre_mean) ** 2 * weight_vm)
    theta_pre_mean = np.sum(previous_theta * weight_theta)
    theta_pre_var = np.sum((previous_theta - theta_pre_mean) ** 2 * weight_theta)

    # sample parameters by the weights computed in last loop.
    population = len(np.where(propose_model == modelindex)[0])
    sample_gamma_index = np.random.choice(previous_bestfitted_index, population, p=weight_gamma)
    sample_a_index = np.random.choice(previous_bestfitted_index, population, p=weight_a)
    sample_nu_index = np.random.choice(previous_bestfitted_index, population, p=weight_nu)
    sample_vm_index = np.random.choice(previous_bestfitted_index, population, p=weight_vm)
    sample_theta_index = np.random.choice(previous_bestfitted_index, population, p=weight_theta)

    # mean of the sample for gamma
    propose_gamma0 = params[sample_gamma_index, 0]
    # draw new gamma with mean and variance
    propose_gamma = abs(np.random.normal(propose_gamma0, np.sqrt(2 * gamma_pre_var)))
    # mean of the sample for a
    propose_a0 = params[sample_a_index, 1]
    # draw new a with mean and variance
    propose_a = abs(np.random.normal(propose_a0, np.sqrt(2 * a_pre_var)))
    # mean of the sample for nu
    propose_nu0 = params[sample_nu_index, 2]
    # draw new nu with mean and variance
    propose_nu = abs(np.random.normal(propose_nu0, np.sqrt(2 * nu_pre_var)))
    # mean of the sample for Vm
    propose_vm0 = params[sample_vm_index, 3]
    # draw new Vm with mean and variance
    propose_vm = abs(np.random.normal(propose_vm0, np.sqrt(2 * vm_pre_var)))
    # mean of the sample for theta
    propose_theta0 = params[sample_theta_index, 4]
    # draw new theta with mean and variance
    propose_theta = abs(np.random.normal(propose_theta0, np.sqrt(2 * theta_pre_var)))

    extend_weight_gamma = weight_gamma[previous_bestfitted_index.searchsorted(sample_gamma_index)]
    extend_weight_a = weight_a[previous_bestfitted_index.searchsorted(sample_a_index)]
    extend_weight_nu = weight_nu[previous_bestfitted_index.searchsorted(sample_nu_index)]
    extend_weight_vm = weight_vm[previous_bestfitted_index.searchsorted(sample_vm_index)]
    extend_weight_theta = weight_theta[previous_bestfitted_index.searchsorted(sample_theta_index)]

    # compute new weights for gamma and a
    weight_gamma_denominator = np.sum(extend_weight_gamma * norm.pdf(propose_gamma, propose_gamma0,
                                                                           np.sqrt(2 * gamma_pre_var)))
    weight_gamma_numerator = norm.pdf(propose_gamma, gamma_pre_mean, np.sqrt(2 * gamma_pre_var))  # Estimation script is using gamma_prior_mean and gamma_prior_var
    weight_gamma = weight_gamma_numerator / weight_gamma_denominator

    weight_a_denominator = np.sum(extend_weight_a * norm.pdf(propose_a, propose_a0,
                                                                   np.sqrt(2 * a_pre_var)))
    weight_a_numerator = norm.pdf(propose_a, a_pre_mean, np.sqrt(2 * a_pre_var))
    weight_a = weight_a_numerator / weight_a_denominator

    weight_nu_denominator = np.sum(extend_weight_nu * norm.pdf(propose_nu, propose_nu0,
                                                                     np.sqrt(2 * nu_pre_var)))
    weight_nu_numerator = norm.pdf(propose_nu, nu_pre_mean, np.sqrt(2 *nu_pre_var))
    weight_nu = weight_nu_numerator / weight_nu_denominator

    weight_vm_denominator = np.sum(extend_weight_vm * norm.pdf(propose_vm, propose_vm0,
                                                                     np.sqrt(2 * vm_pre_var)))
    weight_vm_numerator = norm.pdf(propose_vm, vm_pre_mean, np.sqrt(2 *vm_pre_var))
    weight_vm = weight_vm_numerator / weight_vm_denominator

    weight_theta_denominator = np.sum(extend_weight_theta * norm.pdf(propose_theta, propose_theta0,
                                                                     np.sqrt(2 * theta_pre_var)))
    weight_theta_numerator = norm.pdf(propose_theta, theta_pre_mean, np.sqrt(2 *theta_pre_var))
    weight_theta = weight_theta_numerator / weight_theta_denominator

    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma = weight_gamma / sum(weight_gamma)
    weight_a = weight_a / sum(weight_a)
    weight_nu = weight_nu / sum(weight_nu)
    weight_vm = weight_vm / sum(weight_vm)
    weight_theta = weight_theta / sum(weight_theta)

    return weight_gamma,weight_a,weight_nu,weight_vm,weight_theta,\
           propose_gamma,propose_a,propose_nu,propose_vm,propose_theta