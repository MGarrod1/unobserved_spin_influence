"""

Code used to create figures 1a and 1b.
1) Distribution of control values in a network plot
2) Plot of control against value of the degree.
3) Histograms of the control.

M Garrod, Dec 2019.

"""

import matplotlib.pyplot as plt
from ising_block_level_influence import two_block_sbm_class as two_block

#Seed the random number generators:
import random
import numpy as np
seed = 1
random.seed(seed)
np.random.seed(seed)

two_block_class = two_block.two_block_sbm_analysis(nodes_per_block=250)
sbm_graph = two_block_class.make_sbm()
two_block_class.get_critical_temp(sbm_graph)

Field_Budget = 10.0
beta_factor = 1.5
Max_Iterations = 1000
step_size=50.0
block_control , mag_vals_block = two_block_class.block_level_control(beta_factor,Field_Budget,step_size=step_size)
full_control , mag_vals_full = two_block_class.full_graph_control(sbm_graph,beta_factor ,Field_Budget,Max_Iterations=Max_Iterations,step_size=step_size)

print("block control = {}".format(block_control))
print("\n\nFull control = {}".format(full_control))

print("Block control completed in {} iterations".format(len(mag_vals_block)-1))
print("Full control completed in {} iterations".format(len(mag_vals_full)-1))


#Plot Mags as a function of iterations to demonstrate convergence:
plt.figure(2,figsize=(8,8))
plt.plot(mag_vals_full,'bo')
plt.plot(mag_vals_full,'b')
plt.xlabel("Iterations")
plt.ylabel("Magnetization")
plt.savefig("Plots/mag_as_iterations_full_control")

plt.figure(1,figsize=(8,8))
plt.plot(mag_vals_block,'bo')
plt.plot(mag_vals_block,'b')
plt.xlabel("Iterations")
plt.ylabel("Magnetization")
plt.savefig("Plots/mag_as_iterations_block_control")


two_block_class.plot_full_control_on_network(sbm_graph,full_control,file_path='Plots/full_control_on_graph.jpg',label='(a)')
two_block_class.plot_control_vs_degrees(sbm_graph,full_control,block_control,file_path='Plots/control_vs_degrees')
two_block_class.plot_control_histogram(full_control,block_control,file_path='Plots/control_hist.jpg',label='(b)')
