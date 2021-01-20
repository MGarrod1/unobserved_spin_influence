"""

Computes the magnetization as a function of the
inverse temperature for a three block SBM with
external fields using:
a) Block level mean-field approximation
b) Full graph mean-field approximation
c) Monte Carlo simulations with Metropolis
Hastings dynamics

Variant which focuses on studying the impact of
different realisations of the full graph.

M Garrod, Dec 2019.

"""

import three_block_sbm_class as ThreeBlock
import numpy as np
import random

def block_aligned_state(block_size,a1,a2,a3) :
    vector = np.concatenate(( a1*np.ones(block_size), a2*np.ones(block_size), a3*np.ones(block_size)))
    return vector

def spin_block_aligned_state(block_size,a1,a2,a3) :
    vector = np.concatenate(( np.random.choice([-1,1],block_size,p=[a1,1.0-a1])  , np.random.choice([-1,1],block_size,p=[a2,1.0-a2]) ,
                              np.random.choice([-1,1],block_size,p=[a3,1.0-a3]) ) )
    return vector

class three_block_phase_simulator :

    def __init__(self,coupling_matrix,block_size,B) :

        self.block_size=block_size
        self.three_block = ThreeBlock.three_block_sbm_analysis(coupling_matrix=coupling_matrix, block_size=block_size)
        self.three_block.set_background_field(B)
        self.sbm_graph = self.three_block.make_sbm()
        self.three_block.get_critical_temp(self.sbm_graph)

        self.min_beta_factor = 0.5
        self.max_beta_factor = 5.0

        self.N=len(self.sbm_graph)

    def make_block_level_phase_diag_data(self,f_path) :

        num_block_inits = 25
        num_block_betas = 100
        initial_conditions = [ np.random.uniform(-1.0, 1.0, 3) for k in range(num_block_inits) ]
        beta_vals_block = np.linspace(self.max_beta_factor, self.min_beta_factor, num_block_betas)
        mag_data_block = self.three_block.sample_block_level_phase_diagram(beta_vals_block, initial_conditions)
        mag_data_block.to_csv(f_path)

    def make_full_level_phase_diag_data(self,f_path) :
        num_full_initials = 50
        num_full_betas = 100
        a_vals = [np.random.uniform(-1.0, 1.0, 3) for k in range(num_full_initials)]
        initials = [block_aligned_state(self.block_size, a[0], a[1], a[2]) for a in a_vals]
        beta_vals_full = np.linspace(self.max_beta_factor, self.min_beta_factor, num_full_betas)
        mag_data_MF = self.three_block.sample_full_mf_phase_diagram(self.sbm_graph,beta_vals_full,initials,f_path=f_path)
        mag_data_MF.to_csv(f_path)

    def make_mc_phase_diag_data(self,f_path) :

        num_MC_initials = 25
        T = 10000
        T_Burn = 120000
        MC_Runs = 2
        num_bins = 50

        # Discrete equivalent of the sampling strategy used above:
        a_vals = [np.random.uniform(0.0, 1.0, 3) for k in range(num_MC_initials)]
        initial_conditions = [spin_block_aligned_state(self.block_size, a[0], a[1], a[2]) for a in a_vals]

        beta_vals = np.linspace(self.max_beta_factor, self.min_beta_factor, 20)
        mag_data_mc, mean_mag_hists = self.three_block.sample_mc_phase_diagram(self.sbm_graph, beta_vals, initial_conditions, T,
                                                                          T_Burn, MC_Runs, num_bins=num_bins)
        mag_data_mc.to_csv(f_path)


if __name__ == "__main__" :

    # Seed the random number generators:
    seed = 1
    random.seed(seed)
    np.random.seed(seed)

    f_path_block = "Data/block_level_phase_data.csv"
    f_path_full_graph = "Data/full_MF_phase_data.csv"
    f_path_mc = "Data/MC_phase_data.csv"

    #Parameters:
    min_beta_factor = 0.5
    max_beta_factor = 5.0
    B=5.0
    block_size=400
    coupling_matrix=np.asarray([[10.0, 2.5,0.0], [2.5, 7.5 ,2.5],[0.0,2.5,10.0 ]])

    three_block_phase_simulator = three_block_phase_simulator(coupling_matrix,block_size,B)

    # three_block_phase_simulator.make_full_level_phase_diag_data(f_path_full_graph)
    # three_block_phase_simulator.make_block_level_phase_diag_data(f_path_block)
    three_block_phase_simulator.make_mc_phase_diag_data(f_path_mc)




