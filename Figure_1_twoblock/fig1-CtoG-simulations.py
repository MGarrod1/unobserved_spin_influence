"""

Derive full graph and block level
mean-field Ising controls for two
block SBM for a range of field
budgets and inverse temperatures.

We then evaluate the controls for each
of the different values.

M Garrod, Jan 2020.

"""


import random
import numpy as np
from ising_block_level_influence import N_Block_sbm_class as NBlock

#Seed the random number generators:
seed = 1
random.seed(seed)
np.random.seed(seed)

N_Block=250
coupling_matrix = np.asarray([[10.0,2.5],[2.5,2.5]])
block_sizes=[N_Block,N_Block]
block_background_field=np.asarray([0.0,0.0])

block_system = NBlock.block_mf_ising_system(coupling_matrix,block_sizes,block_background_field)


sbm_graph=block_system.make_sbm() #Sample a particular SBM from the ensemble

ising_analysis = NBlock.ising_analysis(sbm_graph, coupling_matrix, block_sizes, block_background_field)

MC_Sims=15
spin_alignment=1.0 #Start at fully aligned metastable state.

#H_vals =  [  N_Block*(10**k) for k in np.arange(0.0, 1.8, 0.1)]
H_vals =  [  N_Block*(10**k) for k in np.linspace(0.0,2.0,20)]

save_file_prefix="two_block_markup"

#The figure only compares block level and full IIM controls.
ising_analysis.controls_to_get = {'no control':True,
                                'uniform control':True,
                                'NC control':False,
                                'SV control':False,
                                'Block control':True,
                                'Full control':True,
                                'AOB control':False,
                                'Degree control':False}


for beta_factor in [0.5,1.2,1.5] :

    #Block control optimization:
    ising_analysis.gamma = 1.0
    ising_analysis.tol = 1E-5
    ising_analysis.max_mf_fp_iterations = 10000
    ising_analysis.mf_fp_init_state = spin_alignment*np.ones(len(ising_analysis.block_sizes))
    ising_analysis.mf_fp_noisy = False

    ising_analysis.max_mf_iim_iterations = 3000
    ising_analysis.mf_iim_tolerance = 1E-8
    ising_analysis.mf_iim_step_size = 1.0
    ising_analysis.mf_iim_init_control = 'uniform'
    ising_analysis.mf_iim_noisy = True

    # Full control optimization:
    ising_analysis.full_mf_system.gamma = 1.0
    ising_analysis.full_mf_system.tol = 1E-5
    ising_analysis.full_mf_system.max_mf_fp_iterations = 10000
    ising_analysis.full_mf_system.mf_fp_init_state = spin_alignment * np.ones(len(ising_analysis.full_graph))
    ising_analysis.full_mf_system.mf_fp_noisy = False

    ising_analysis.full_mf_system.max_mf_iim_iterations = 1000
    ising_analysis.full_mf_system.mf_iim_step_size = 50.0
    ising_analysis.full_mf_system.mf_iim_tolerance = 1E-6
    ising_analysis.full_mf_system.mf_iim_init_control = 'uniform'
    ising_analysis.full_mf_system.mf_iim_noisy = False

    #MC Parameters
    ising_analysis.T = 20000
    ising_analysis.T_Burn = 10000
    ising_analysis.MC_Runs = MC_Sims
    ising_analysis.eval_initial_state = spin_alignment * np.ones(len(ising_analysis.full_graph))

    ising_analysis.H_sweep_data_fname = "Data/{}_data_spins{}_bf_{}".format(save_file_prefix,
        round(spin_alignment, 0),beta_factor).replace('.', '-')
    ising_analysis.H_sweep_diagnostics_fname = "Data/{}_diagnostics_spins{}_bf_{}".format(save_file_prefix,
        round(spin_alignment, 0),beta_factor).replace('.', '-')

    ising_analysis.make_h_sweep_data(beta_factor, H_vals)

ising_analysis.save_iim_eval_parameters("Data/{}_params.csv".format(save_file_prefix))