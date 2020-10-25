"""

Functions used for processing the
Pokec social network dataset.

Contains functions for reading in the
data and extracting block structure of
the graph based on specified age and
regions of users profiles.


M Garrod, Dec 2019.

"""

import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm
import geopy.distance


def read_reduced_profile_info(path_to_profiles) :
    """
    Read in the reduced profile information which has been pre-processed to only
    include the information about the node attributes which we are interested in.

    """
    print("Reading profile information...")
    profile_info = pd.read_csv(path_to_profiles)
    return profile_info


def read_edge_list(path_to_edge_list) :
    print("Reading edge data...")
    edge_data = pd.read_csv(path_to_edge_list,sep='\t',names=['id_1','id_2'],dtype=np.int32)
    return edge_data

def geo_dist(coords_1,coords_2) :
    return geopy.distance.geodesic(coords_1, coords_2).km



def make_region_counts_data(profile_info) :

    """

    Counts the number of users within each
    region within the dataset.

    Parameters
    -----------

    profile_info : pandas dataframe

    Dataframe containing a 'region' attribute
    as profile information.

    Returns
    ----------

    region_count_data : pandas dataframe

    dataframe containing the region names and
    counts of users within each region.

    """

    region_names = list(profile_info['region'])
    region_name_set = list(set(region_names))

    #Create dataframe containing the user counts within each region:
    region_count_data = pd.DataFrame(columns=["region", "num users"])
    for name in tqdm(region_name_set):
        df_to_append = pd.DataFrame({"region": [name], "num users": [region_names.count(name)]})
        region_count_data = region_count_data.append(df_to_append)

    return region_count_data


def make_edge_data_table(profile_info,edge_data) :
    """
    Given the edge list and a specified set of
    user profiles we construct a dataframe containing
    the list of edges and the attributes of the nodes
    at each end of the respective edges.

    Parameters
    ------------

    profile_info : pandas dataframe

    Specified set of user profiles with corresponding user
    ids.

    edge_edge : pandas dataframe

    edge list in pandas dataframe format.

    Parameters
    ------------

    edge_data_with_info : pandas dataframe

    Each corresponds to an edge in the graph and
    contains te attributes of each of the users
    involved.

    """
    print("Merging edges and profile information...")

    edge_data_with_info = edge_data.merge(profile_info,left_on="id_1",right_on="user_id")
    edge_data_with_info = edge_data_with_info.merge(profile_info,left_on='id_2',right_on='user_id')
    edge_data_with_info=edge_data_with_info.drop(columns=['Unnamed: 0_x','Unnamed: 0_y','user_id_x','user_id_y'])
    edge_data_with_info.columns = ['id_1','id_2','gender_1','region_1','age_1','u1_lat','u1_lon','gender_2','region_2','age_2','u2_lat','u2_lon']

    return edge_data_with_info


def make_graph_from_edge_data(profile_info,edge_data_with_info,add_disconnected_nodes=True) :

    """
    Given a pandas dataframe with each row corresponding
    to an edge we can make a networkx graph.

    Parameters
    -------------

    edge_data_with_info : pandas dataframe

    Each corresponds to an edge in the graph and
    contains te attributes of each of the users
    involved.

    add_disconnected_nodes : bool

    If False then we exclude nodes which
    are not involved in any connections.

    Returns
    -------------

    graph : networkx graph

    """

    graph = nx.DiGraph()

    graph_user_ids = list(profile_info['user_id'])

    if add_disconnected_nodes == True :
        for id_val in graph_user_ids :
            graph.add_node(id_val)

    for id_1,id_2 in zip( list(edge_data_with_info['id_1']) , list(edge_data_with_info['id_2'])) :
        graph.add_edge(id_1,id_2)

    for id_1,id_2 in zip( list(edge_data_with_info['id_2']) , list(edge_data_with_info['id_1'])) :
        graph.add_edge(id_1,id_2)

    graph = graph.to_undirected(reciprocal=False)

    return graph


def get_undirected_edges_between(edge_data_with_info, user_ids_1, user_ids_2):
    """

    Given a pandas dataframe containing the
    edges associated with a (perhaps) directed
    graph and two sets of users
    returns the set of edges of the
    equivalent undirected graph.

    Parameters
    ------------

    edge_data_with_info : pandas dataframe

    user_ids_1 : list

    list of user ids for the first group

    user_ids_2 : list

    list of user ids for the second group


    Returns
    ----------

    shared_pairs : list

    List of tuples assoicated with the
    edges in the network.

    edges_12 : list

    List of tuples for edges from block
    1 to block 2

    edges_21 : list

    List of tuples for edges from block
    2 to block 1

    """

    edges_12 = edge_data_with_info.loc[
        (edge_data_with_info['id_1'].isin(user_ids_1)) &
        (edge_data_with_info['id_2'].isin(user_ids_2))]
    id_pairs_12 = [(i, j) for i, j in zip(list(edges_12['id_1']), list(edges_12['id_2']))]
    edges_21 = edge_data_with_info.loc[
        (edge_data_with_info['id_1'].isin(user_ids_2)) &
        (edge_data_with_info['id_2'].isin(user_ids_1))]
    id_pairs_21 = [(j, i) for i, j in zip(list(edges_21['id_1']), list(edges_21['id_2']))]

    shared_pairs = id_pairs_12
    for pair in id_pairs_21:
        if pair not in shared_pairs:
            shared_pairs.append(pair)

    return id_pairs_12, id_pairs_21, shared_pairs


def block_label_tuple_to_string(block_tuple) :
    """

    Parameters
    ------------

    block_tuple : tuple

    Tuple containing the region name
    (str) and (pair of ages)

    Returns
    --------------

    block_str : str

    String containing reduced info.

    """

    name = block_tuple[0]
    name.split('-')[1]

    ages = block_tuple[1]


    return name.split('-')[1].title() + "\n{}-{}".format(int(ages[0]),int(ages[1]))


def make_blau_coord_combinations(region_names,age_ranges) :
    blau_coord_combinations = []
    for i in region_names:
        for j in age_ranges:
            blau_coord_combinations.append((i, j))
    return blau_coord_combinations


def make_coupling_graph(block_interaction_data,blau_coord_combinations,on_attribute='coupling'):

    """

    Make a networkx graph from block interaction
    data.

    Assumes pre computed blau coord combinations.

    Parameters
    -------------

    block_interaction_data : pandas dataframe

    Dataframe containing numbers of edges and interaction
    strengths between different blocks within a larger graph.

    blau:coord_combinations : list

    List of tuples containing the different region name age range
    paris

    on_attribute : str (opt)

    String specifying which variable to take as the edge weight.


    """
    #coupling_graph = nx.Graph()
    coupling_graph = nx.OrderedGraph()

    for pair_1 in tqdm(blau_coord_combinations):
        for pair_2 in blau_coord_combinations:
            try :
                coupling = list(block_interaction_data.loc[(block_interaction_data['cat_1'] == pair_1) & (block_interaction_data['cat_2'] == pair_2)][on_attribute])[0]
            except :
                pdb.set_trace()

            coupling_graph.add_edge(pair_1, pair_2, weight=coupling)

    return coupling_graph


def get_region_age_profile_ids(reduced_profiles,region_name, age_range):

    """

    Extract the ids of the profiles with a specified region
    name and age range.

    Paramrters
    -----------

    region_name : str

    Name of the region to be extracted

    age_range : tuple

    tuple containing a minimum and maximum
    age.

    Returns
    ----------

    user_ids : list

    list of user ids for the specific region and age range.

    """


    profile_set = reduced_profiles[(reduced_profiles['region'] == region_name) & (reduced_profiles['age'].between(age_range[0], age_range[1]))]

    user_ids = list(profile_set['user_id'])

    return user_ids

class pokec_graph:

    def __init__(self, region_names, age_ranges, subsample_size=None):

        """

        Initialize the class by reading
        in the data.

        Orders the profiles within the respective blocks.
        This makes processing later on easier.

        """

        path_to_reduced_profiles = 'Data/reduced_pokec_profiles.csv'
        path_to_edge_list = 'Data/raw_data/soc-pokec-relationships.txt'

        # Read in relevant data:
        profile_info = read_reduced_profile_info(path_to_reduced_profiles)
        self.edge_data = read_edge_list(path_to_edge_list)

        # Find the combinations for the reduced profiles:
        self.make_blau_coord_combinations(region_names, age_ranges)

        # Take only the nodes containing valid ages and regions:
        self.profiles_w_data = profile_info.loc[(profile_info['age'] > 0.0) & (profile_info['user_lat'] > 0.0)]

        self.reduced_profiles = pd.DataFrame()
        print("Extracting blocks...")
        for coord in self.blau_coord_combinations:
            this_df = self.profiles_w_data[(self.profiles_w_data['region'] == coord[0]) & (
                self.profiles_w_data['age'].between(coord[1][0], coord[1][1]))]
            self.reduced_profiles = self.reduced_profiles.append(this_df)

        if subsample_size is not None:
            if len(self.reduced_profiles) < subsample_size:
                print(
                    "Found only {} profiles (less than {}). Taking all in list.".format(len(self.reduced_profiles),
                                                                                        subsample_size))
            else:
                self.reduced_profiles = self.reduced_profiles.sample(n=subsample_size)

    def make_graph_on_regions(self, region_names, add_disconnected_nodes=True, lcc_only=False):

        """

        Make the subgraph of users within a specified
        set of regions in the data.

        Parameters
        ------------

        region_names : list

        List of strings containing region names


        Returns
        ------------

        graph : networkx graph

        graph of the edges

        """

        edge_data_with_info = make_edge_data_table(self.reduced_profiles, self.edge_data)
        graph  = make_graph_from_edge_data(self.reduced_profiles, edge_data_with_info,
                                              add_disconnected_nodes=add_disconnected_nodes)

        if lcc_only == True:
            print("Extracting LCC...")

            #lcc_graph = max(nx.connected_component_subgraphs(graph), key=len) #Now depreciated.
            Gcc = sorted(nx.connected_components(graph), key=len, reverse=True)
            lcc_graph = graph.subgraph(Gcc[0])


            lcc_set = list(lcc_graph.nodes())
            self.reduced_profiles = self.reduced_profiles.loc[self.reduced_profiles['user_id'].isin(lcc_set)]
            edge_data_with_info = make_edge_data_table(self.reduced_profiles, self.edge_data)
            graph = make_graph_from_edge_data(self.reduced_profiles, edge_data_with_info,
                                                  add_disconnected_nodes=add_disconnected_nodes)

        return graph, edge_data_with_info

    def get_region_age_profile_ids(self, reduced_profiles, region_name, age_range):

        """

        Extract the ids of the profiles with a specified region
        name and age range.

        Paramrters
        -----------

        region_name : str

        Name of the region to be extracted

        age_range : tuple

        tuple containing a minimum and maximum
        age.

        Returns
        ----------

        user_ids : list

        list of user ids for the specific region and age range.

        """

        profile_set = reduced_profiles[
            (reduced_profiles['region'] == region_name) & (reduced_profiles['age'].between(age_range[0], age_range[1]))]

        user_ids = list(profile_set['user_id'])

        return user_ids

    def make_blau_coord_combinations(self, region_names, age_ranges):
        self.blau_coord_combinations = []
        for i in region_names:
            for j in age_ranges:
                self.blau_coord_combinations.append((i, j))

        return self.blau_coord_combinations

    def block_sizes(self, reduced_profiles):
        block_sizes = [len(self.get_region_age_profile_ids(reduced_profiles, coord[0], coord[1])) for coord in
                       self.blau_coord_combinations]
        return block_sizes

    def make_block_interaction_data(self, graph, reduced_profiles, edge_data_with_info, region_names, age_ranges):

        """

        Get the block level interaction data for
        a specific set of profiles given a set of
        region names and age bins.


        Parameters
        -------------

        reduced_profiles : pandas dataframe

        profile information

        edge_data_with_info : pandas dataframe

        profile information

        region_names : list

        list of strings

        age_ranges : list

        list of tuples containing pairs of numbers

        Returns
        -------------

        block_interaction_data : pandas dataframe


        """

        block_interaction_data = pd.DataFrame()
        N = len(reduced_profiles)

        blau_coord_combinations = self.make_blau_coord_combinations(region_names, age_ranges)

        for pair_1 in tqdm(blau_coord_combinations):  # (original symmetric assuming code...)
            for pair_2 in blau_coord_combinations:

                user_ids_1 = self.get_region_age_profile_ids(reduced_profiles, pair_1[0], pair_1[1])
                user_ids_2 = self.get_region_age_profile_ids(reduced_profiles, pair_2[0], pair_2[1])
                N1 = len(user_ids_1)
                N2 = len(user_ids_2)

                average_age_block_1 = np.mean(
                    list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_1)]['age']))
                average_age_block_2 = np.mean(
                    list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_2)]['age']))
                age_diff = abs(average_age_block_2 - average_age_block_1)

                lat_b1 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_1)]['user_lat'])[0]
                lng_b1 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_1)]['user_lon'])[0]

                lat_b2 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_2)]['user_lat'])[0]
                lng_b2 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_2)]['user_lon'])[0]

                spatial_dist = geo_dist((lat_b1, lng_b1), (lat_b2, lng_b2))


                edges_12, edges_21, shared_pairs = get_undirected_edges_between(edge_data_with_info, user_ids_1,user_ids_2)
                if pair_1 == pair_2:
                    edges_between = int(len(shared_pairs) / 2)
                else:
                    edges_between = int(len(shared_pairs))

                if pair_1 == pair_2:
                    edge_density = (2.0 * edges_between) / (N1 * (N1 - 1.0))
                    edges_per_node = N1 * edge_density
                    coupling = N1 * edge_density

                else:
                    edge_density = float(edges_between) / (N1 * N2)
                    edges_per_node = ((N1 * N2 * edge_density) / N)

                    coupling = N1 * edge_density

                edges_between_w_graph = [p for p in nx.edge_boundary(graph, user_ids_1, nbunch2=user_ids_2)]
                edges_bet_graph = len(edges_between_w_graph)

                current_df = pd.DataFrame({'cat_1': [pair_1], 'cat_2': [pair_2], 'N_C1': [N1], 'N_C2': [N2],
                                           'edges_between': [edges_between], 'edges_between_graph': [edges_bet_graph],
                                           'frac_between': [edge_density], 'edges_per_node': [edges_per_node]
                                              , 'coupling': [coupling],
                                           'edges_12': [len(edges_12)], 'edges_21': [len(edges_21)], 'age difference':
                                               [age_diff], 'distance (km)': [spatial_dist]
                                              , 'b1 mean age': [average_age_block_1],
                                           'b2 mean age': [average_age_block_2],
                                           'b1_lat': [lat_b1], 'b1_lng': [lng_b1], 'b2_lat': [lat_b2],
                                           'b2_lng': [lng_b2]})

                block_interaction_data = block_interaction_data.append(current_df)

        return block_interaction_data

    def make_age_interaction_data(self, reduced_profiles, edge_data_with_info, age_ranges):

        """

        Parameters
        -------------

        reduced_profiles : pandas dataframe

        profile information

        edge_data_with_info : pandas dataframe

        profile information

        age_ranges : list

        list of tuples containing pairs of numbers

        Returns
        -------------

        age_interaction_data : pandas dataframe

        """

        age_interaction_data = pd.DataFrame()
        N = len(reduced_profiles)

        for pair_1 in tqdm(age_ranges):
            for pair_2 in age_ranges:

                user_ids_1 = list(
                    reduced_profiles.loc[reduced_profiles['age'].between(pair_1[0], pair_1[1])]['user_id'])
                user_ids_2 = list(
                    reduced_profiles.loc[reduced_profiles['age'].between(pair_2[0], pair_2[1])]['user_id'])

                N1 = len(user_ids_1)
                N2 = len(user_ids_2)

                age_diff = abs(np.mean(pair_1) - np.mean(pair_2))

                edges_between_data = edge_data_with_info.loc[
                    (edge_data_with_info['id_1'].isin(user_ids_1)) &
                    (edge_data_with_info['id_2'].isin(user_ids_2))]

                edges_between = len(edges_between_data)

                if pair_1 == pair_2:
                    edge_density = (2.0 * edges_between) / (N1 * (N1 - 1.0))
                    edges_per_node = N1 * edge_density
                    coupling = N1 * edge_density

                else:
                    edge_density = float(edges_between) / (N1 * N2)
                    edges_per_node = ((N1 * N2 * edge_density) / N)

                    coupling = N1 * edge_density

                current_df = pd.DataFrame({'cat_1': [pair_1], 'cat_2': [pair_2], 'N_C1': [N1], 'N_C2': [N2],
                                           'edges_between': [edges_between],
                                           'frac_between': [edge_density], 'edges_per_node': [edges_per_node]
                                              , 'coupling': [coupling], 'age difference': [age_diff]})

                age_interaction_data = age_interaction_data.append(current_df)

        return age_interaction_data

    def make_distance_interaction_data(self, reduced_profiles, edge_data_with_info, region_names):

        """

        Get the block level interaction data for
        a specific set of profiles given a set of
        region names and age bins.


        Parameters
        -------------

        reduced_profiles : pandas dataframe

        profile information

        edge_data_with_info : pandas dataframe

        profile information

        region_names : list

        list of strings

        age_ranges : list

        list of tuples containing pairs of numbers

        Returns
        -------------

        block_interaction_data : pandas dataframe


        """

        distance_interaction_data = pd.DataFrame()
        N = len(reduced_profiles)

        for pair_1 in tqdm(region_names):
            for pair_2 in region_names:

                user_ids_1 = list(reduced_profiles.loc[reduced_profiles['region'] == pair_1]['user_id'])
                user_ids_2 = list(reduced_profiles.loc[reduced_profiles['region'] == pair_2]['user_id'])
                N1 = len(user_ids_1)
                N2 = len(user_ids_2)

                lat_b1 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_1)]['user_lat'])[0]
                lng_b1 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_1)]['user_lon'])[0]

                lat_b2 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_2)]['user_lat'])[0]
                lng_b2 = list(reduced_profiles.loc[reduced_profiles['user_id'].isin(user_ids_2)]['user_lon'])[0]

                spatial_dist = geo_dist((lat_b1, lng_b1), (lat_b2, lng_b2))

                edges_between_data = edge_data_with_info.loc[
                    (edge_data_with_info['id_1'].isin(user_ids_1)) &
                    (edge_data_with_info['id_2'].isin(user_ids_2))]

                edges_between = len(edges_between_data)

                if pair_1 == pair_2:
                    edge_density = (2.0 * edges_between) / (N1 * (N1 - 1.0))
                    edges_per_node = N1 * edge_density
                    coupling = N1 * edge_density

                else:
                    edge_density = float(edges_between) / (N1 * N2)
                    edges_per_node = ((N1 * N2 * edge_density) / N)

                    coupling = N1 * edge_density

                current_df = pd.DataFrame({'cat_1': [pair_1], 'cat_2': [pair_2], 'N_C1': [N1], 'N_C2': [N2],
                                           'edges_between': [edges_between],
                                           'frac_between': [edge_density], 'edges_per_node': [edges_per_node]
                                              , 'coupling': [coupling], 'distance (km)': [spatial_dist],
                                           'b1_lat': [lat_b1], 'b1_lng': [lng_b1], 'b2_lat': [lat_b2],
                                           'b2_lng': [lng_b2]})

                distance_interaction_data = distance_interaction_data.append(current_df)

        return distance_interaction_data



    def make_coupling_graph(self, block_interaction_data, on_attribute='coupling'):
        """

        Assumes pre computed blau coord combinations.

        """
        # coupling_graph = nx.Graph()

        coupling_graph = nx.OrderedGraph()

        for pair_1 in tqdm(self.blau_coord_combinations):
            for pair_2 in self.blau_coord_combinations:


                coupling = list(block_interaction_data.loc[(block_interaction_data['cat_1'] == pair_1) & (
                                block_interaction_data['cat_2'] == pair_2)][on_attribute])[0]

                coupling_graph.add_edge(pair_1, pair_2, weight=coupling)

        return coupling_graph


    def make_block_info(self) :

        """

        Construct a dataframe containing information
        about the blocks (e.g average ages)

        """
        block_ids = [ pair for pair in self.blau_coord_combinations]
        block_sizes = [len(self.get_region_age_profile_ids(self.reduced_profiles, pair_1[0], pair_1[1]))
                       for pair_1 in self.blau_coord_combinations]
        profile_id_sets = [self.get_region_age_profile_ids(self.reduced_profiles, pair_1[0], pair_1[1])
                           for pair_1 in self.blau_coord_combinations]
        mean_block_ages = [np.mean(list(self.reduced_profiles.loc[self.reduced_profiles['user_id'].isin(id_set)]['age']))
            for id_set in profile_id_sets]

        block_data = pd.DataFrame({'Block' : block_ids, 'Size' : block_sizes, 'Mean Age' : mean_block_ages})
        return block_data



def generate_3_region_Bratislava(age_ranges,user_threshold=10000) :

    """

    Constructs the LCC graph and the block interaction data
    matrix for regions in Bratislava with specifies age
    ranges.

    Dataframes are saved to .csv files.

    Parameters
    -------------

    age_ranges : list

    List of tuples containing pairs of ascending age ranges.
    Assume that they should be non overlapping

    user_threshold : int

    Only take regions which contain more than this threshold of users.

    """

    path_to_reduced_profiles = 'Data/reduced_pokec_profiles.csv'

    #Profile read in:
    all_profiles = read_reduced_profile_info(path_to_reduced_profiles)
    all_profiles=all_profiles.loc[(all_profiles['age'] > 0.0) & (all_profiles['user_lat'] > 0.0)]

    #Get profiles in Bratislava and find the corresponding region names:
    bratislav_profiles = all_profiles.loc[ all_profiles['region'].str.contains("bratislava") ]
    bratislav_region_counts = make_region_counts_data(bratislav_profiles)
    region_names = list(bratislav_region_counts.loc[bratislav_region_counts['num users']>user_threshold]['region'])

    #Construct the Pokec class in order to make the graph:
    pokec_class = pokec_graph(region_names,age_ranges)
    graph, edge_data_with_info = pokec_class.make_graph_on_regions(region_names,lcc_only=True)
    print("Number of nodes in LCC graph = {}".format(len(graph)))

    block_interaction_data = pokec_class.make_block_interaction_data(graph,pokec_class.reduced_profiles,edge_data_with_info,region_names,age_ranges)

    nx.write_graphml(graph,"Data/Bratislava_graph.graphml")
    block_interaction_data.to_csv("Data/Bratislava_Blocks_agecats_{}.csv".format(len(age_ranges)))

    block_data=pokec_class.make_block_info()
    block_data.to_csv("Data/Bratislava_Block_Data.csv")

    pokec_class.reduced_profiles.to_csv("Data/Bratislava_Profiles.csv")

    #Make the coupling graph:
    coupling_graph = pokec_class.make_coupling_graph(block_interaction_data,on_attribute='coupling')
    nx.write_graphml(coupling_graph,'Data/Bratislava_coupling.graphml')
