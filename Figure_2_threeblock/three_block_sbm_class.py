"""

Code to perform simulations of the Ising model on a three block
SBM with chain-like structure.

Opposing external fields are on the two end blocks.

M Garrod, Sep 2019.

"""

#Python functions:
import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm as tqdm

#My own functions:
from ising_block_level_influence import mean_field_IIM
from spatial_spin_monte_carlo import spatial_spin_monte_carlo as Spins


class three_block_sbm_analysis :


    def __init__(self,coupling_matrix=np.asarray([[10.0, 2.5,0.0], [2.5, 7.5 ,2.5],[0.0,2.5,10.0 ]]),block_size=400) :

        self.N1=block_size
        self.N2=block_size
        self.N3=block_size
        self.sizes = [self.N1, self.N2, self.N3]
        self.N = np.sum(self.sizes)
        self.N_blocks = len(self.sizes)
        self.coupling_matrix = coupling_matrix
        self.coupling_graph = nx.from_numpy_matrix(self.coupling_matrix)



    def set_background_field(self,Magnitude) :
        self.Field_Magnetide = Magnitude
        B1 = self.Field_Magnetide
        B2 = 0.0
        B3 = -self.Field_Magnetide
        self.background_field = np.concatenate((B1 * np.ones(self.N1), B2 * np.ones(self.N2), B3 * np.ones(self.N3)))
        self.background_field_block = np.asarray([B1, B2, B3])


    def prob_matrix_from_coupling(self) :
        prob_mat = np.zeros( ( len(self.sizes),len(self.sizes) ))
        for i in range(self.N_blocks) :
            for j in range(self.N_blocks) :

                if i == j :
                    prob_mat[i][j] = self.coupling_matrix[i][j]/self.sizes[i]
                else :
                    prob_mat[i][j] = (self.coupling_matrix[i][j]*self.N)/(2.0*self.sizes[i]*self.sizes[j])

        return prob_mat


    def make_sbm(self) :
        prob_mat = self.prob_matrix_from_coupling()
        sbm_graph = nx.stochastic_block_model(self.sizes, prob_mat)

        return sbm_graph

    def get_critical_temp(self,sbm_graph):
        self.beta_c = Spins.crit_beta(sbm_graph)
        print("beta c = {}".format(self.beta_c))

    def sample_block_level_magnetization_series(self,sbm_graph,T,T_Burn,beta_factor ,Initial_State=None) :


        """

        Returns time series of magnetisations obtained
        from Monte Carlo simulations for the three
        different blocks in the system.

        Parameters
        ------------

        sbm_graph : networkx graph

        T : int

        Number of timesteps to run for.

        T_Burn : int

        Burn-in time for the Monte-Carlo simulations

        beta_factor : float

        multiple of the critical temperature to
        sample at.

        Initial_State : numpy array

        state to initialize the Ising simulations in.



        """

        beta = beta_factor * self.beta_c
        spin_series = Spins.Run_MonteCarlo(sbm_graph, T, beta, T_Burn=T_Burn, positions=None, Initial_State=Initial_State,control_field=self.background_field)
        Mag_series_B1 = [ np.mean(k) for k in np.transpose( np.transpose(spin_series)[0:self.N1] )  ]
        Mag_series_B2 = [np.mean(k) for k in  np.transpose( np.transpose(spin_series)[self.N1:self.N1+self.N2] )  ]
        Mag_series_B3 = [np.mean(k) for k in np.transpose( np.transpose(spin_series)[self.N1 + self.N2 : self.N1 + self.N2 + self.N3] )  ]
        return Mag_series_B1 , Mag_series_B2 , Mag_series_B3



    def sample_block_level_magnetization_means(self,sbm_graph,T,T_Burn,beta_factor ,Initial_State=None):

        """

        Use Monte Carlo simulations to estimate the magnetization
        at the level of blocks.


        Parameters
        ------------

        sbm_graph : networkx graph

        T : int

        Number of timesteps to run for.

        T_Burn : int

        Burn-in time for the Monte-Carlo simulations

        beta_factor : float

        multiple of the critical temperature to
        sample at.

        Initial_State : numpy array

        state to initialize the Ising simulations in.

        Returns
        -------------




        """

        Mag_series_B1, Mag_series_B2, Mag_series_B3 = self.sample_block_level_magnetization_series(sbm_graph, T, T_Burn,
                                                                                                   beta_factor,
                                                                                                   Initial_State=Initial_State)
        block_magnetizations = [np.mean(Mag_series_B1), np.mean(Mag_series_B2), np.mean(Mag_series_B3)]

        return block_magnetizations


    def mean_field_block_level_averages(self,sbm_graph,beta_factor) :

        """

        Estimate the magnetization of network using the mean-field
        approximation.

        Use full-mean field and block level approximations.


        Returns
        ----------

        mf_magnetization_data : pandas dataframe.

        """

        beta = beta_factor*self.beta_c

        mf_magnetization_data = pd.DataFrame()



        #m_full = mean_field_IIM.mean_field_magnetization_weighted_graph(sbm_graph, self.background_field, beta)
        mf_full_system = mean_field_IIM.mean_field_ising_system(sbm_graph, self.background_field)
        m_full = mf_full_system.mf_magnetization(self.background_field,beta)

        #m_block = mean_field_IIM.mean_field_magnetization_weighted_graph(self.coupling_graph,self.background_field_block,beta)

        mf_block_system = mean_field_IIM.mean_field_ising_system(self.coupling_graph,self.background_field_block)
        m_block = mf_block_system.mf_magnetization(self.background_field_block,beta)

        MF_Full_B1 = np.mean(m_full[0:self.N1])
        MF_Full_B2 = np.mean(m_full[self.N1:self.N1 + self.N2])
        MF_Full_B3 = np.mean(m_full[self.N1 + self.N2:self.N1 + self.N2 + self.N3])


        mf_magnetization_data = mf_magnetization_data.append(pd.DataFrame({'MF_Full_B1' : [MF_Full_B1],'MF_Full_B2' : [MF_Full_B2],
                                                                           'MF_Full_B3' : [MF_Full_B3],'MF_Block_B1' : [m_block[0]],
                                                                           'MF_Block_B2' : [ m_block[1]] , 'MF_Block_B3' : [m_block[2]],
                                                                           'beta_factor' : [ beta_factor]}))

        return mf_magnetization_data


    def sample_block_level_phase_diagram(self,beta_vals_block,initial_conditions) :

        """

        Compute the magnetisation using the block level
        mean field approximation for a range of
        different beta values and initial conditions
        in order to contain a phase diagram for the system.

        Parameters
        -------------

        beta_vals_block : numpy array

        list or array of values of the multiplier
        of the critical temperature.

        This is multiplied by the critical temp
        of the full graph, rather than the critical
        temperature associated with the block coupling matrix.

        initial_conditions : list

        list of arrays specifiying initial conditions for the
        iterative procedure.

        Returns
        ------------

        mag_data_block : pandas dataframe

        Data about the magnetisations for each beta
        and each initial condition.

        """
        mag_data_block = pd.DataFrame()

        mf_block_system = mean_field_IIM.mean_field_ising_system(self.coupling_graph, self.background_field_block)
        mf_block_system.mf_fp_noisy=False

        for beta_f in tqdm(beta_vals_block,position=0,leave=True):
            beta = beta_f*self.beta_c
            # For each beta use all possible initial conditions:
            for index, initial_state in enumerate(initial_conditions):
                mf_block_system.mf_fp_init_state=initial_state
                m_block = mf_block_system.mf_magnetization(self.background_field_block, beta)

                """
                m_block = mean_field_IIM.mean_field_magnetization_weighted_graph(self.coupling_graph,
                                                                                 self.background_field_block,
                                                                                 beta, initial_state=initial_state,
                                                                                 gamma=1.0,noisy=False)
                
                """


                df_to_append = pd.DataFrame({'beta_factor': [beta_f], 'Mean_mag_block': [np.mean(m_block)],
                                             'MB1_block': [m_block[0]], 'MB2_block': [m_block[1]],
                                             'MB3_block': [m_block[2]], 'init_state_num': [index],
                                             'initial_state': [initial_state]})

                mag_data_block = mag_data_block.append(df_to_append)

        return mag_data_block



    def sample_full_mf_phase_diagram(self,sbm_graph,beta_vals,initial_conditions,f_path="Data/full_MF_phase_data.csv") :


        """

        Compute the magnetisation using the full graph
        mean field approximation for a range of
        different beta values and initial conditions
        in order to contain a phase diagram for the system.

        Parameters
        -------------

        beta_vals_block : numpy array

        list or array of values of the multiplier
        of the critical temperature.

        initial_conditions : list

        list of arrays specifiying initial conditions for the
        iterative procedure.

        Returns
        ------------

        mag_data_block : pandas dataframe

        Data about the magnetisations for each beta
        and each initial condition.

        """
        mag_data = pd.DataFrame()
        A = nx.to_numpy_matrix(sbm_graph)
        sbm_graph = nx.from_numpy_matrix(A)
        mf_full_system = mean_field_IIM.mean_field_ising_system(sbm_graph, self.background_field)

        #Fixed point iteration parameters:
        mf_full_system.mf_fp_noisy=False
        mf_full_system.gamma = 1.0
        mf_full_system.tol = 1E-5
        mf_full_system.max_mf_fp_iterations = 10000
        mf_full_system.mf_fp_init_state = 'aligned'

        for beta_f in tqdm(beta_vals,position=0,leave=True,desc='Beta vals'):
            beta = beta_f*self.beta_c
            # For each beta use all possible initial conditions:
            index = 0
            for initial_state in tqdm( initial_conditions ,desc = 'Random initial conds' ):

                mf_full_system.mf_fp_init_state=initial_state
                m_full = mf_full_system.mf_magnetization(self.background_field, beta)

                df_to_append = pd.DataFrame({'beta_factor': [beta_f], 'Mean_mag_MF': [np.mean(m_full)],
                                             'MB1_MF': [np.mean(m_full[0:self.N1])], 'MB2_MF': [np.mean(m_full[self.N1:self.N1+self.N2])],
                                             'MB3_MF': [np.mean(m_full[self.N1+self.N2:self.N1+self.N2+self.N3])], 'init_state_num': [index],
                                             'initial_state': [initial_state]})

                mag_data = mag_data.append(df_to_append)
                index+=1
                mag_data.to_csv(f_path)

        return mag_data


    def sample_mc_phase_diagram(self,sbm_graph,beta_vals,initial_conditions,T,T_Burn,MC_Runs,num_bins=50):

        """

        Compute the magnetisation using Monte Carlo
        simulations for a range of
        different beta values and initial conditions
        in order to contain a phase diagram for the system.

        Returns
        ---------

        mag_data_block : pandas dataframe

        Data about the magnetisations for each beta
        and each initial condition.

        mean_mag_hists : list

        list of arrays conntaining a histograms
        of the magnetisation for each value of
        beta considered.


        """
        mag_data_mc = pd.DataFrame()

        mean_mag_hists = []
        for beta_f in tqdm(beta_vals,leave=True,position=0):

            #List to store all of the sampled magnetisations
            #For a particular beta value to make the hsitogram
            fixed_beta_magnetisation_samps = []

            for run in tqdm(range(MC_Runs),leave=True,position=1):
                for index, Initial_State in enumerate(initial_conditions):
                    B1, B2, B3 = self.sample_block_level_magnetization_series(sbm_graph, T, T_Burn, beta_f,Initial_State=Initial_State)

                    #Magnetisation time series for the whole graph:
                    overall_mag = [np.mean(p) for p in np.transpose(np.asarray([B1, B2, B3]))]

                    df_to_append = pd.DataFrame({'beta_factor': [beta_f], 'Mean_mag_MC': [np.mean(overall_mag)],
                                                 'MB1_MC': [np.mean(B1)], 'MB2_MC': [np.mean(B2)],
                                                 'MB3_MC': [np.mean(B3)], 'init_state_num': [index],
                                                 'initial_state': [Initial_State],'magnetisation_series' : [overall_mag]})

                    mag_data_mc = mag_data_mc.append(df_to_append)
                    fixed_beta_magnetisation_samps = np.concatenate((fixed_beta_magnetisation_samps, overall_mag))

            hist, bins = np.histogram(fixed_beta_magnetisation_samps, bins=np.linspace(-1.0, 1.0, num_bins))
            mean_mag_hists.append(hist)

        return mag_data_mc , mean_mag_hists






