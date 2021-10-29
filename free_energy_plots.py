# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 17:20:34 2021

@author: Lauren
"""

import matplotlib.pyplot as plt
import numpy as np
from sys import argv, exit

def main(argv):
    """
    Make free energy plots for a given solute.
    
    Input: python free_energy_plots.py solute solvent1 solvent2 ...
    For each solute-solvent combination, generates two plots--one with the free
    energy change between individual lambda steps and a second for the 
    cummulative free enegy change at each lambda step. 
    
    """
    # Check command line input, requires a solute and at least one solvent
    if len(argv) < 3:
        print("Usage: python free_energy_plots.py solute solvent1 solvent2 ...")
        exit(1)
        
    # Assign command line arguements
    solute = argv[1]
    solvents = []
    for solvent in argv[2:]:
        solvents.append(solvent)
        
    # Dictionary for solute name from abbreviation
    names = {"PF7" : r"Phenyl$-$CF$_2-$PF$_5$", \
             "NF7" : r"Naphthalene$-$CF$_2-$PF$_5$"}
    
    # Dictionary containing bar output as bar[solvent]["delta_Gs","delta_stds",
    # "cummulative", "lambdas", "l", "total_dg", "total_std"]
    data = free_energy_output(solute, solvents)
    
    # Subtract the energy related to the counterion
    ion_data = free_energy_output("NA", solvents)
    
    # Apply the finite size correction by adding the Wigner self-interaction energy term
    ion_data = finite_size_correction(solvents, ion_data)
    
    # Generate the relevant plots of the 
    free_energy_plots(solute, solvents, data, ion_data, names)
    
    return None

class Bar():
    """Store important outputs obtained from the "gmx bar ..." module."""
    
    def __init__(self, solvent):
        self.solvent = solvent
        self.dG = []
        self.dSTD = []
        self.cummulative = []
        self.lambdas = []
        self.total_dG = 0
        self.total_STD = 0
        self.size_correct = False

def free_energy_output(solute, solvents):
    """
    Extract the data from "bar.out" for each solute-solvent.
    
    Reads bennett acceptance ratio output for each solute-solvent combination
    and store it in data[solute] as a "Bar" object.

    Parameters
    ----------
    solute : str
        Name of the solute
    solvents : (str) list
        Names of the solvents 

    Returns
    -------
    bar : (Bar) dict
        Dictionary containing Bar objects with solvents as keys

    """
    # Grab bennett acceptance ratio output for each solvent and store in
    # data[solute] as a "Bar" object""
    files = []
    data = {}
    for solvent in solvents:
        path = "{0}_{1}/bar.out".format(solute,solvent)
        files.append(path)
    
    for filename, solv in zip(files, solvents):    
        with open(filename) as file: 
            bar_out = file.readlines()
        bar_data = Bar(solv)
        
        for i, line in enumerate(bar_out):
            if "Final results in kJ/mol:" in line:
                initial = i
            if "total" in line:
                total_dg = float(line.split("DG")[1].split()[0])
                total_std = float(line.split()[-1])
                break
        delta_Gs, delta_stds, cummulative, lambdas, l = [], [], [0], [], 0
        for line in bar_out[initial:]:
            if "point" in line:
                dg = line.split("DG")[1].split()[0]
                std = line.split()[-1]
                delta_Gs.append(float(dg))
                delta_stds.append(float(std))
                cummulative.append(cummulative[-1] + delta_Gs[-1])
                l += 1
                lambdas.append(l)
        bar_data.dG, bar_data.dSTD = delta_Gs, delta_stds
        bar_data.cummulative, bar_data.lambdas = cummulative, lambdas
        bar_data.total_dG, bar_data.total_std = total_dg, total_std
        data[solv] = bar_data
    
    return data

def finite_size_correction(solvents, data):
    """
    This correction is not negligible for solvents with low dielectric constants
    
    Includes the self-interaction between periodic images of the ion and the 
    uniform background charge as well as undersolvation. The volume during the
    first 5 lambda steps is entered in manually.

    """
    # dielectric of H2O: https://doi.org/10.1021/acs.jcim.8b00026
    # dielectric of DCM: https://doi.org/10.1021/ct200731v
    dielectrics = {"water" : 96.2, "DCM" : 10.1}
    cubic_constant = 2.837297
    # 1/(4*pi*(permittivity of vacuum))
    constant_frac = 1 / (4 * np.pi * 8.854187 * 10 ** (-12)) 
    elem_charge = 1.602176 * 10 ** (-19)
    charges = [1,0.75,0.50,0.25,0.0]
    volumes = {"water" : [31.8994,31.8988,31.9122,31.9035,31.9031],
               "DCM" : [31.1522,31.1884,31.2256,31.2713,31.2483]}
    
    for solvent in solvents:
    
        cummulative_correction = 0
        data[solvent].size_correct = True
        
        for i in range(4):
            length_ave = (volumes[solvent][i] + volumes[solvent][i+1]) / 2
            length = length_ave ** (1/3) * 10 ** (-9)
            original = data[solvent].dG[i]
            charge_diff = charges[i+1] ** 2 - charges[i] ** 2 
            charge_diff = charge_diff * (elem_charge ** 2) * 6.02214 * 10 ** (23)
            #charge = charges[i] ** 2 * (elem_charge ** 2) * 6.02214 * 10 ** (23)
            correction = constant_frac * cubic_constant * charge_diff / 2 / dielectrics[solvent] / length / 1000
            cummulative_correction -= correction
            data[solvent].dG[i] = original - correction
            data[solvent].cummulative[i+1] += cummulative_correction
        
        for k, j in enumerate(data[solvent].cummulative[5:]):
            data[solvent].cummulative[k+5] = j - cummulative_correction
            
        print("The finite size correction for Na in {0} is {1:.2} kJ/mol".format(\
               solvent, cummulative_correction))
        
    return data

def free_energy_plots(solute, solvents, data, data_ion, names):
    
    # Make figure with 2 subplots per solvent
    fig, axes = plt.subplots(len(solvents), 2, figsize=(12,6*len(solvents)), \
                             constrained_layout=True)
    
    # Set up variables
    font = {'color': 'black', 'weight': 'semibold', 'size': 16}
    solvents2x = [x for pair in zip(solvents, solvents) for x in pair]
    coulomb_axis = range(1,6,1)
    lj_axis = range(5,20,1)
    totals = {}
    for solvent in solvents:
        raw = data[solvent].cummulative[-1]
        ion = data_ion[solvent].cummulative[-1]
        totals[solvent] = -(raw - ion) 
        print("raw is {0}, ion is {1} and total is {2} in {3}".format(raw, ion, totals[solvent], solvent))
    deldelG = totals[solvents[1]] - totals[solvents[0]]
    logk = - deldelG / (298 * 8.314/1000 * 2.3026)
    if data_ion["water"].size_correct:
        correction = "with size correction"
    else:
        correction = "without corrections"
    #title = "{0}; log K = {1}".format(names[solute],np.round(logk,2))
    title = "log K = {0}; {1}".format(np.round(logk,2), correction)
    solv_color = {"water" : "skyblue", "DCM" : "thistle"}
    
    # For even plot_type: lambda vs. dG_lambda; 
    # otherwise: lambda vs dG(cummulative) 
    plot_type = 0
    
    for solvent, ax in zip(solvents2x, axes.flat):
        
        label = "solvent : {0}".format(solvent)
        # Access solvent data within the data dictionary
        bar = data[solvent]
        bar_ion = data_ion[solvent]
        
        # Plots for free energy change between lambda states    
        if plot_type % 2 == 0:
            
            # Color code area under curve wrt. Coulomb or LJ contribution
            dG_diff_coulomb = [i - j for i, j in zip(bar.dG[:5], bar_ion.dG[:5])] 
            dG_diff_lj = [i - j for i, j in zip(bar.dG[4:], bar_ion.dG[4:])]
            #dG_diff_coulomb = bar.dG[:5]
            #dG_diff_lj = bar.dG[4:]
            ax.fill_between(coulomb_axis, dG_diff_coulomb, \
                            color="plum", label="Coulomb")
            ax.fill_between(lj_axis, dG_diff_lj, \
                            color="paleturquoise", label="Lennard-Jones")
            
            dG_diff = bar.dG
            dSTD_sum = bar.dSTD
            dG_diff = [i - j for i, j in zip(bar.dG, bar_ion.dG)]
            dSTD_sum = [i + j for i, j in zip(bar.dSTD, bar_ion.dSTD)]
            # Plot curve between lambda steps with error bar  
            ax.plot(bar.lambdas, dG_diff, color="#292929")
            ax.scatter(bar.lambdas, dG_diff, color="red", marker="o", s=30)
            ax.errorbar(bar.lambdas, dG_diff, yerr=dSTD_sum, \
                        color="#292929", capsize=5, elinewidth=2, capthick=2,\
                        ecolor="red")
                
            # Label y-axis, add legend, label plot as "solute in solvent"
            ax.set_ylabel(r"$\Delta G_{\lambda}$ (kJ/mol)", fontdict=font)
            ax.legend(fontsize=14, loc=1)
            ax.text(0.96, 0.79, label, color="black", fontsize=14, \
                    transform=ax.transAxes, bbox=dict(facecolor=solv_color[solvent], \
                    edgecolor='dimgrey', boxstyle='round,pad=0.5', \
                    alpha=0.7), ha="right", va="bottom")
         
        # Plots for cummulative energy change wrt lambda
        else:
            
            # Plot the coulomb data
            #G_sum_coulomb = bar.cummulative[:5]
            G_sum_coulomb = [i - j for i, j in zip(bar.cummulative[:5], \
                                                   bar_ion.cummulative[:5])]
            ax.plot([0] + bar.lambdas[:4], G_sum_coulomb, \
                    color="mediumvioletred")
            ax.scatter([0] + bar.lambdas[:4], G_sum_coulomb, \
                       color="mediumvioletred", marker= "^", s=80, label='Coulomb')
                
            # Plot the LJ data
            #G_sum_lj = bar.cummulative[4:]
            G_sum_lj = [i - j for i, j in zip(bar.cummulative[4:], \
                                                   bar_ion.cummulative[4:])]
            ax.plot(bar.lambdas[3:], G_sum_lj, color="darkcyan")
            ax.scatter(bar.lambdas[3:], G_sum_lj, color="darkcyan", \
                       marker= "o", s=80, label='Lennard-Jones')
            
            # Plot dashed line at net energy change
            #total_G = bar.cummulative[-1]
            total_G = bar.cummulative[-1] - bar_ion.cummulative[-1]
            ax.plot(range(0,20,1), [total_G] * 20, \
                    linewidth=3, linestyle="dashed", color="#292929")
            
            # Plot settings
            ax.set_ylabel(r"Total $\Delta G$ (kJ/mol)", fontdict=font)
            ax.legend(fontsize=14, loc=4)
            
            # Label the net energy change on the plot
            delG = r"$\Delta G = ${0} $\pm$ {1} kJ$/$mol".format(str(np.round(total_G,1)),\
                   str(np.round(bar.total_std + bar_ion.total_std, 1))) 
            ax.annotate("", xy=(bar.lambdas[8], total_G), xycoords='data', size=20, \
                        xytext=(bar.lambdas[8], 0), \
                        arrowprops=dict(arrowstyle="<->"))
            ax.text(bar.lambdas[8], (total_G) / 2, delG, size=18, \
                    bbox=dict(boxstyle="round", fc="w", edgecolor='dimgrey', \
                    ec="0.5", alpha=0.7), ha="center", color="midnightblue")
            
            # Label plot as "solute"
            ax.text(0.96, 0.17, label, color='black', fontsize=14, \
                    va="bottom", transform=ax.transAxes, ha="right", \
                    bbox=dict(facecolor=solv_color[solvent], edgecolor='dimgrey', \
                    boxstyle='round,pad=0.5', alpha=0.9))
        
        # General plot settings for x-axis, ticks, border spines and grid
        ax.set_xticks(range(0,21,3))
        ax.set_xlim(0,19)
        ax.tick_params(axis='y', labelsize=16, direction='out', width=3, \
                       length=7, pad=10)
        ax.tick_params(axis='x', labelsize=16, direction='out', width=3, \
                       length=7, pad=10)
        ax.grid()
        for i in ["top","bottom","left","right"]:
            ax.spines[i].set_linewidth(2)
        if solvent in solvents[-1]:
            ax.set_xlabel(r"$\lambda$ states", fontdict=font, labelpad=10)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
            
        plot_type += 1
    
    plt.suptitle(title, fontsize=24)
    if data_ion["water"].size_correct:
        plt.savefig("free_energy_plots_{}_sized.png".format(solute))
    else:
        plt.savefig("free_energy_plots_{}.png".format(solute))
    
    return None

if __name__ ==  '__main__':
    main(argv)