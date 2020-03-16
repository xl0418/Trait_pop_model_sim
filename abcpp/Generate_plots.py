def generate_plots(result,plots_dir):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import sys
    import matplotlib
    matplotlib.use('Agg')
    est_data = np.load(result,allow_pickle=True).item()
    population = int(len(est_data['model_data'][0]) / 3)
    fitness = est_data['fitness'][-2]

    q5_TVP = np.argsort(fitness[:population])[-int(population // 200 + 1)]  # best 5%
    q5_TV = np.argsort(fitness[population:2 * population])[
                -int(population // 200 + 1)] + population  # best 5%
    q5_TVM = np.argsort(fitness[2 * population:])[-int(population // 200 + 1)] + 2 * population  #
    # best 5%

    fit_index_TVP = np.where(fitness[:population] > fitness[q5_TVP])[0]
    fit_index_TV = np.where(fitness[population:2 * population] > fitness[q5_TV])[0] + population
    fit_index_TVM = np.where(fitness[2 * population:] > fitness[q5_TVM])[0] + 2 * population

    previous_bestfitted_index_TVP = fit_index_TVP
    previous_bestfitted_index_TV = fit_index_TV - population
    previous_bestfitted_index_TVM = fit_index_TVM - 2 * population

    top_pop = [len(previous_bestfitted_index_TVP),len(previous_bestfitted_index_TV),
               len(previous_bestfitted_index_TVM)]

    gamma_vec_TVP = est_data['gamma_data_TVP'][-2, previous_bestfitted_index_TVP]
    a_vec_TVP = est_data['a_data_TVP'][-2, previous_bestfitted_index_TVP]
    nu_vec_TVP = est_data['nu_data_TVP'][-2, previous_bestfitted_index_TVP]
    vm_vec_TVP = est_data['vm_data_TVP'][-2, previous_bestfitted_index_TVP]
    theta_vec_TVP = est_data['theta_data_TVP'][-2, previous_bestfitted_index_TVP]

    gamma_vec_TV = est_data['gamma_data_TV'][-2, previous_bestfitted_index_TV]
    a_vec_TV = est_data['a_data_TV'][-2, previous_bestfitted_index_TV]
    nu_vec_TV = est_data['nu_data_TV'][-2, previous_bestfitted_index_TV]
    vm_vec_TV = est_data['vm_data_TV'][-2, previous_bestfitted_index_TV]
    theta_vec_TV = est_data['theta_data_TV'][-2, previous_bestfitted_index_TV]

    gamma_vec_TVM = est_data['gamma_data_TVM'][-2, previous_bestfitted_index_TVM]
    a_vec_TVM = est_data['a_data_TVM'][-2, previous_bestfitted_index_TVM]
    nu_vec_TVM = est_data['nu_data_TVM'][-2, previous_bestfitted_index_TVM]
    vm_vec_TVM = est_data['vm_data_TVM'][-2, previous_bestfitted_index_TVM]
    theta_vec_TVM = est_data['theta_data_TVM'][-2, previous_bestfitted_index_TVM]

    # Log the progress
    output_log = sys.stdout
    f = open('Plotting_log.txt', 'w')
    sys.stdout = f

    print('TVP 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
        np.mean(gamma_vec_TVP) ,
        np.mean(a_vec_TVP) , np.mean(nu_vec_TVP) ,
        np.mean(vm_vec_TVP), np.mean(theta_vec_TVP)))
    print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
        np.var(gamma_vec_TVP) ,
        np.var(a_vec_TVP) , np.var(nu_vec_TVP) ,
        np.var(vm_vec_TVP), np.var(theta_vec_TVP)))

    print('TV 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
        np.mean(gamma_vec_TV) ,
        np.mean(a_vec_TV) , np.mean(nu_vec_TV),
        np.mean(vm_vec_TV), np.mean(theta_vec_TV)))
    print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
        np.var(gamma_vec_TV) ,
        np.var(a_vec_TV) , np.var(nu_vec_TV) ,
        np.var(vm_vec_TV), np.var(theta_vec_TV)))

    print('TVM 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
        np.mean(gamma_vec_TVM) ,
        np.mean(a_vec_TVM) , np.mean(nu_vec_TVM) ,
        np.mean(vm_vec_TVM), np.mean(theta_vec_TVM)))
    print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
        np.var(gamma_vec_TVM),
        np.var(a_vec_TVM) , np.var(nu_vec_TVM) ,
        np.var(vm_vec_TVM), np.var(theta_vec_TVM)))
    print("")
    print("Generating plots is done!!")

    sys.stdout = output_log
    f.close()
    est_para = ['gamma', 'alpha', 'nu', 'vm', 'theta']  # ,'vm'
    model_para = ['AWC', 'UWC', 'MWC']
    est_array = np.concatenate([gamma_vec_TVP, a_vec_TVP, nu_vec_TVP, vm_vec_TVP, theta_vec_TVP,
                                gamma_vec_TV, a_vec_TV, nu_vec_TV, vm_vec_TV, theta_vec_TV,
                                gamma_vec_TVM, a_vec_TVM, nu_vec_TVM, vm_vec_TVM, theta_vec_TVM])
    est_label = np.repeat(np.tile(est_para,3),np.tile(top_pop,5))
    model_label = np.repeat(model_para, np.array(top_pop) * len(est_para))

    ss_list = {'est': est_array, 'est_label': est_label,
               'model_label': model_label}
    ss_df = pd.DataFrame(ss_list)
    # plt.ioff()
    vioplot = sns.catplot(x="model_label", y="est", col="est_label",palette=["#CB1B45","#FAD689","#0D5661"],
                          data=ss_df, kind="box", height=5, aspect=.6, sharey=False)

    vioplot.set_axis_labels("", "Estimate value")
    # vioplot._legend.set_title('$h^2$')
    axes = vioplot.axes.flatten()
    axes[0].set_title("$\gamma $")
    axes[1].set_title("$\\alpha $")
    axes[2].set_title("$\\nu $")
    axes[3].set_title("$V_m $")
    axes[4].set_title("$\\theta $")

    for no_plots in range(0,4):
        axes[no_plots].ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)

    vioplot.savefig(plots_dir+'/ms_plots.png')

