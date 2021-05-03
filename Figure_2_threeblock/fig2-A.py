"""

Code for figure 2A.

Visualize a 3 block SBM with
fixed background fields.

M Garrod, Dec 2019.

"""

import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt
import three_block_sbm_class as ThreeBlock


def plot_three_block(sbm_graph,three_block,fname="three_block_plot",color_on='background_field',label=None) :

    N1 = N2 = N3 = 400
    pos_block_1 = [ [ np.random.uniform(0,1) , np.random.uniform(0,1.0)] for k in range(N1)]
    pos_block_2 = [ [ np.random.uniform(1.5,3.0) , np.random.uniform(0,1.0)] for k in range(N2)]
    pos_block_3 = [ [ np.random.uniform(4.0,5.0) , np.random.uniform(0,1.0)] for k in range(N2)]
    pos=np.concatenate((pos_block_1,pos_block_2,pos_block_3))

    pos_dict={}
    for k in range(len(sbm_graph)) :
        pos_dict[k]=pos[k]

    optimal_distance = 0.01/(len(sbm_graph)**0.5)
    pos_spring = nx.spring_layout(sbm_graph,iterations=1,pos=pos_dict,k=optimal_distance)
    #pos_spring=pos_dict

    if color_on=='background_field':
        node_cols=three_block.background_field
    else :
        node_cols=color_on

    fig,ax = plt.subplots(1,figsize=(8,6))
    nx.draw_networkx_nodes(sbm_graph,pos=pos_spring,node_size=25,node_color=node_cols,alpha=0.95,cmap='coolwarm')
    nx.draw_networkx_edges(sbm_graph,pos=pos_spring,alpha=0.2)
    plt.axis('off')

    vmin = min(node_cols)
    vmax = max(node_cols)

    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin = vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm,pad=-0.01)

    #Set colourbar labels:
    cbar.set_ticks([-5, -2.5,0.0,2.5 ,5.0])
    cbar.set_ticklabels(['-5','-2.5' ,'0.0', '2.5', '5'])
    cbar.ax.tick_params(labelsize=20)

    if label is not None:
        plt.text(-1.0, 0.18, label, fontsize=25)

    cbar.ax.set_ylabel('Ambient Field', rotation=90,fontsize=16,labelpad=15.0)
    plt.savefig(fname + ".jpg",dpi=100)


if __name__ == "__main__" :
    # Seed the random number generators:
    seed = 2
    random.seed(seed)
    np.random.seed(seed)
    B = 5.0
    three_block = ThreeBlock.three_block_sbm_analysis()
    three_block.set_background_field(B)
    sbm_graph = three_block.make_sbm()
    three_block.get_critical_temp(sbm_graph)
    plot_three_block(sbm_graph, three_block, fname="Plots/three_blocks_background",label='(a)')

