"""

Utilities for running Monte-Carlo simulations
on the Pokec social network.

M Garrod, Oct 2020.

"""

import numpy as np
import tqdm
import math
from spatial_spin_monte_carlo import spatial_spin_monte_carlo as Spins

def Run_MonteCarlo_Average(input_graph, T, beta_factor,beta_c, T_Burn=0,addition_control=None,sampling_method="Metropolis",initial_state=None,full_graph_field=None):

    """

    Samples a sequence of spin states on an Ising system
    for the graph supplied.
    
    More efficient than the version in the class as we don't
    store all of the spin states.

    Parameters
    -------------

    graph : networkx graph

    network structure

    T : int

    number of time steps to run the simulation for

    beta : float

    Inverse temperature

    T_Burn : int (opt)

    Burn in time. We run the dynamics for T+T_Burn timesteps
    and only record samples after T_Burn.

    positions : numpy ndarray (optional)

    Positions of nodes. Needs to be in the same order as the
    node list. If this is specified we use the spatial spins
    class.

    Returns
    ------------

    Spin_Series : ndarray (N_Block x T)

    Array continaing time series of spin values for each
    of the blocks.

    """


    spin_system = Spins.spins(input_graph)

    # Set params:
    spin_system.sampling_method = sampling_method
    spin_system.Beta = beta_factor*beta_c

    # Set the initial state:
    spin_system.node_states = np.copy(initial_state)


    #Set the field:
    spin_system.field = spin_system.applied_field_to_each_node

    if addition_control is not None :
        spin_system.scalar_field = full_graph_field + addition_control
    else :
        spin_system.scalar_field = full_graph_field

    current_block_level_mags = []
    for p in tqdm.tqdm_notebook( range(T_Burn + T) , miniters=int(float(T)/10) ) :
        #for p in range(T_Burn + T) :
        spin_system.do_spin_flip()
        # Copy the array so it is not a pointer:
        current_state = np.copy(spin_system.node_states)
        if p > T_Burn - 1:
            current_block_level_mags.append( np.mean(current_state) )

    return current_block_level_mags


    
def Run_MonteCarlo_Block(input_graph,groups_node_ids ,T, beta_factor,beta_c, T_Burn=0,addition_control=None,sampling_method="Metropolis",initial_state=None,full_graph_field=None):

    """

    Samples a sequence of spin states on an Ising system
    for the graph supplied. Averages the result across
    blocks with supplied group node ids



    Parameters
    -------------

    input_graph : networkx graph

    network structure

    groups_node_ids : list

    List of lists containing the ids for the different
    groups of nodes within the graph. The ids used must
    match those witin the graph.

    T : int

    number of time steps to run the simulation for

    beta : float

    Inverse temperature

    T_Burn : int (opt)

    Burn in time. We run the dynamics for T+T_Burn timesteps
    and only record samples after T_Burn.

    positions : numpy ndarray (optional)

    Positions of nodes. Needs to be in the same order as the
    node list. If this is specified we use the spatial spins
    class.

    Returns
    ------------

    Spin_Series : ndarray (N_Block x T)

    Array continaing time series of spin values for each
    of the blocks.

    """


    spin_system = Spins.spins(input_graph)

    # Set params:
    spin_system.sampling_method = sampling_method
    spin_system.Beta = beta_factor*beta_c

    # Set the initial state:
    spin_system.node_states = np.copy(initial_state)


    #Set the field:
    spin_system.field = spin_system.applied_field_to_each_node

    if addition_control is not None :
        spin_system.scalar_field = full_graph_field + addition_control
    else :
        spin_system.scalar_field = full_graph_field

    current_block_level_mags = []
    for p in tqdm.tqdm_notebook( range(T_Burn + T) , miniters=int(float(T)/10) ) :
        #for p in range(T_Burn + T) :
        spin_system.do_spin_flip()
        # Copy the array so it is not a pointer:
        current_state = np.copy(spin_system.node_states)
        if p > T_Burn - 1:
            #block averages using the indices of nodes:
            block_level_means = [ np.mean(np.take(current_state,indices)) for indices in groups_node_ids ] 
            current_block_level_mags.append( block_level_means )

    block_mag_series = np.asarray(current_block_level_mags)
            

    return current_block_level_mags


def Run_MonteCarlo_Snapshot(input_graph,groups_node_ids ,T, beta_factor,beta_c,frac_to_sample=1.0, T_Burn=0,addition_control=None,sampling_method="Metropolis",initial_state=None,full_graph_field=None):

    """

    Samples the spin state at a random point in time at the level of blocks by sampling from blocks.

    Parameters
    -------------

    input_graph : networkx graph

    network structure

    groups_node_ids : list

    List of lists containing the ids for the different
    groups of nodes within the graph. The ids used must
    match those witin the graph.

    T : int

    number of time steps to run the simulation for

    beta : float

    Inverse temperature

    T_Burn : int (opt)

    Burn in time. We run the dynamics for T+T_Burn timesteps
    and only record samples after T_Burn.

    positions : numpy ndarray (optional)

    Positions of nodes. Needs to be in the same order as the
    node list. If this is specified we use the spatial spins
    class.

    Returns
    ------------

    Spin_Series : ndarray (N_Block x T)

    Array continaing time series of spin values for each
    of the blocks.

    """


    spin_system = Spins.spins(input_graph)

    # Set params:
    spin_system.sampling_method = sampling_method
    spin_system.Beta = beta_factor*beta_c

    # Set the initial state:
    spin_system.node_states = np.copy(initial_state)


    #Set the field:
    spin_system.field = spin_system.applied_field_to_each_node

    if addition_control is not None :
        spin_system.scalar_field = full_graph_field + addition_control
    else :
        spin_system.scalar_field = full_graph_field

    current_block_level_mags = []
    for p in tqdm.tqdm_notebook( range(T_Burn + T) , miniters=int(float(T)/10) ) :
        #for p in range(T_Burn + T) :
        spin_system.do_spin_flip()
        # Copy the array so it is not a pointer:
        current_state = np.copy(spin_system.node_states)
        
        
        
    block_level_means = [ ]
    for indices in groups_node_ids :
        # Compute number of blocks to sample.
        # ceil so that we always sample at least 1 node per block.
        num_to_sample = math.ceil(frac_to_sample*len(indices))
        sample_indices=np.random.choice(indices,num_to_sample,replace=False)
        block_level_means.append(np.mean(np.take(current_state,sample_indices)))
            

    return block_level_means

