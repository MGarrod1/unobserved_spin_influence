"""

Carries out simulations for figures 2c,2d &2e.

This figure explores the behaviour of the optimal
block level control as a function of the field budget (H)
and inverse temperature (beta).

Running the iterative scheme with different initial conditions
allows us to find magnetisations (and hence optimal controls)
associated with the different metastable states of the system.

M Garrod, Dec 2019.

"""

from tqdm import tqdm as tqdm
import numpy as np
import pandas as pd
import networkx as nx

from ising_block_level_influence import N_Block_sbm_class as NBlock
from spatial_spin_monte_carlo import spatial_spin_monte_carlo as Spins


def block_system_param_save(block_system,fname) :
    parameters={}
    parameters["$\gamma_{\mathrm{block}}$"]=block_system.gamma
    parameters["$\\tau_{\mathrm{block}}$"] = block_system.tol
    parameters['$\epsilon_{\mathrm{block}}$'] = block_system.mf_iim_step_size
    parameters['$a_{\mathrm{block}}$'] = block_system.mf_iim_tolerance
    param_df = pd.DataFrame({'Parameter' : list(parameters.keys()), 'Value' : list(parameters.values())} )
    param_df.to_csv(fname,index=False)

#Parameters:
coupling_matrix=np.asarray([[10.0, 2.5,0.0], [2.5, 7.5 ,2.5],[0.0,2.5,10.0 ]])
block_sizes=[400,400,400]
B=5.0
block_background_field=np.asarray([B,0.0,-B])


block_system = NBlock.block_mf_ising_system(coupling_matrix,block_sizes,block_background_field)

#Computing the critical temp through a sampled SBM
#Allows us to mathc up values with the phase diagram
#plot
sbm_graph=block_system.make_sbm()
#beta_c = Spins.crit_beta(nx.from_numpy_matrix(coupling_matrix))
beta_c = Spins.crit_beta(sbm_graph)


#Fixed point iteration parameters:
block_system.gamma=1.0
block_system.tol=1E-5 #Was 10^{-10}???
block_system.max_mf_fp_iterations=10000
block_system.mf_fp_init_state=np.ones(len(block_system.block_sizes))
block_system.mf_fp_noisy=False

#IIM parameters
block_system.max_mf_iim_iterations=1000
block_system.mf_iim_step_size=1.0
block_system.mf_iim_tolerance=1E-10
block_system.mf_iim_init_control='uniform'
block_system.mf_iim_noisy=False

H_Vals = [10**k for k in np.linspace(1.5,4.5,25)]
beta_f_vals = [ 5.0,2.5,1.0]

block_cont_data= pd.DataFrame()

num_block_betas=25
beta_f_vals = np.linspace(0.5, 5.0, num_block_betas)

#beta_f_vals = [ 1.0,2.5,5.0]

Init_state_samps=25
for pp in tqdm(range(Init_state_samps),position=0) :

    block_system.mf_fp_init_state=np.random.uniform(-1.0,1.0,3)


    for beta_factor in tqdm(beta_f_vals,position=1) :

        beta=beta_c*beta_factor

        A=nx.to_numpy_matrix(block_system.graph)
        m, m_seq=block_system.mf_magnetization(block_background_field ,beta,return_sequence=True)
        sus_vec=block_system.mf_magnetization_gradient(m,beta)

        current_df = pd.DataFrame({'beta_factor': [beta_factor],
                                   'init_cond': [ str(list(np.asarray(block_system.mf_fp_init_state))) ],
                                   'magnetisation' : [np.mean(m)],
                                   'm1': [m[0]],
                                   'm2': [m[1]],
                                   'm3': [m[2]],
                                   'susB1':[sus_vec[0]],
                                   'susB2':[sus_vec[1]],
                                   'susB3':[sus_vec[2]],
                                   'm_fp_convergence':[m_seq]})

        block_cont_data = block_cont_data.append(current_df)

block_cont_data.to_csv("Data/three_block_sus_data.csv")

block_system_param_save(block_system,"Data/three_block_control_params.csv")
