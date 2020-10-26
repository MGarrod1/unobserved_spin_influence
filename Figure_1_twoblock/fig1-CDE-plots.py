"""

Plot the markups of the block level and
Full IIM controls as a function of the
Field Budget for different values of
beta.

Used for figures 1c, 1d & 1e.

M Garrod, Dec 2019.

"""

#Python modules:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
import glob

# Converter taken from: https://stackoverflow.com/questions/42755214/how-to-keep-numpy-array-when-saving-pandas-dataframe-to-csv
def from_np_array(array_string):
    array_string = ','.join(array_string.replace('[ ', '[').split())
    return np.array(ast.literal_eval(array_string))

def make_M_vals_plot_plot(mag_mark_data, beta_factor):

    as_h_vals_data = mag_mark_data.loc[mag_mark_data['beta_factor'] == beta_factor]
    M_IIM_Vals = list(as_h_vals_data['M(full)'])
    M_Block_Vals = list(as_h_vals_data['M(block)'])
    M_Unif_Vals = list(as_h_vals_data['M(uniform)'])
    Budget_Vals = list(as_h_vals_data['H'])

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.clf()
    plt.plot(Budget_Vals,M_IIM_Vals,'ro-',label='full graph')
    plt.plot(Budget_Vals, M_Block_Vals, 'bo-',label='block level')
    plt.plot(Budget_Vals, M_Unif_Vals, 'go-', label='uniform')
    plt.xlabel("H", fontsize=20)
    plt.ylabel("$\\Delta M_{MC}$", fontsize=20)
    plt.xscale('log')

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.legend(fontsize=15, loc='upper right')

    plt.savefig("Plots/Mag_as_H_beta_f_{}".format(beta_factor).replace('.', '-'), bbox_inches='tight')

def make_markup_plot(mag_mark_data, beta_factor,label=None):
    as_h_vals_data = mag_mark_data.loc[mag_mark_data['beta_factor'] == beta_factor]

    M_IIM_Vals = list(as_h_vals_data['M(full)'])
    M_IIM_SEs = list(as_h_vals_data['M(full)_SE'])

    M_Block_Vals = list(as_h_vals_data['M(block)'])
    M_Block_SEs = list(as_h_vals_data['M(block)_SE'])

    M_Unif_Vals = list(as_h_vals_data['M(uniform)'])
    M_Unif_SEs = list(as_h_vals_data['M(uniform)_SE'])

    Budget_Vals = list(as_h_vals_data['H'])

    Full_Markups = [(i - j) for i, j in zip(M_IIM_Vals, M_Unif_Vals)]
    # Standard error propagation formula on the SEs:
    Full_Markup_SEs = [(i ** 2 + j ** 2) ** 0.5 for i, j in zip(M_IIM_SEs, M_Unif_SEs)]

    Block_Markups = [(i - j) for i, j in zip(M_Block_Vals, M_Unif_Vals)]
    # Standard error propagation formula on the SEs:
    Block_Markup_SEs = [(i ** 2 + j ** 2) ** 0.5 for i, j in zip(M_Block_SEs, M_Unif_SEs)]

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.clf()
    plt.errorbar(Budget_Vals, Full_Markups, yerr=Full_Markup_SEs, fmt='bo', label="Full Graph")
    plt.plot(Budget_Vals, Full_Markups, 'b')

    plt.errorbar(Budget_Vals, Block_Markups, yerr=Block_Markup_SEs, fmt='rs', label="Block Level")
    plt.plot(Budget_Vals, Block_Markups, 'r')

    plt.xlabel("H", fontsize=20)
    plt.ylabel("$\\Delta M_{MC}$", fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xscale('log')
    plt.ylim(0.0, 0.055)

    if label is not None:
        plt.text(230, 0.045, label, fontsize=25)
    if label == '(e)\n$\\beta=1.2\\beta_c$ (hot)' :
        plt.legend(fontsize=15, loc='upper right')

    plt.savefig("Plots/makrup_beta_f_{}".format(beta_factor).replace('.', '-') +".jpg" , bbox_inches='tight')


def make_control_behaviour_plot(mag_mark_data,beta_factor) :

    Block_Size = 250

    as_h_vals_data = mag_mark_data.loc[mag_mark_data['beta_factor'] == beta_factor]
    block_controls = list(as_h_vals_data['block control'])
    full_controls = list(as_h_vals_data['full control'])
    Budget_Vals = list(as_h_vals_data['H'])

    #Convert block control on all N nodes to the total control on each block:
    blocker=[ [Block_Size*np.mean(block_controls[k][0:250]), Block_Size*np.mean(block_controls[k][250:])] for k in range(len(block_controls)) ]
    b1_frac_block = [i[0]/j for i, j in zip(blocker, Budget_Vals)]
    b1_frac_full = [np.sum(i[0:Block_Size]) / np.sum(i) for i in full_controls]

    plt.figure(figsize=(8, 8))
    plt.plot(Budget_Vals, b1_frac_full, 'bo-', label='MF-IIM (Full)')
    plt.plot(Budget_Vals, b1_frac_block, 'ro-', label='MF-IIM Block')

    # 50/50 as a guide for the eye
    plt.plot(Budget_Vals, 0.5 * np.ones(len(Budget_Vals)), 'k--', label='Uniform')
    plt.xscale('log')
    plt.ylim(0, 1)
    plt.xlabel("H", fontsize=20)
    plt.ylabel("Fraction on block 1", fontsize=20)
    plt.legend(loc='upper right', fontsize=15)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig("Plots/control_fractions_f_{}".format(beta_factor).replace('.', '-'), bbox_inches='tight')

mag_mark_data = pd.DataFrame()
file_names = glob.glob('Data/two_block_markup_data*.csv')
for fname in file_names:
    curr_dat = pd.read_csv(fname, converters={'block control': from_np_array, 'full control': from_np_array})
    mag_mark_data = mag_mark_data.append(curr_dat)


for beta_factor,label in zip([0.5,1.2,1.5],['(c)\n$\\beta=0.5\\beta_c$ (cold)','(d)\n$\\beta=1.2\\beta_c$ (critical)','(e)\n$\\beta=1.2\\beta_c$ (hot)']) :
    make_M_vals_plot_plot(mag_mark_data, beta_factor)
    make_markup_plot(mag_mark_data,beta_factor,label=label)
    make_control_behaviour_plot(mag_mark_data, beta_factor)
