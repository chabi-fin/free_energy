# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 16:08:35 2021

@author: Lauren
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# time given in ps, temp in K
temp_data = np.loadtxt("temperature_{}.xvg".format(sys.argv[1]), comments=["#","@"])
temperature1_data = np.loadtxt("temperature1_{}.xvg".format(sys.argv[1]), comments=["#","@"])
temperature2_data = np.loadtxt("temperature2_{}.xvg".format(sys.argv[1]), comments=["#","@"])
pressure1_data = np.loadtxt("pressure1_{}.xvg".format(sys.argv[1]), comments=["#","@"])
pressure2_data = np.loadtxt("pressure2_{}.xvg".format(sys.argv[1]), comments=["#","@"])
density1_data = np.loadtxt("density1_{}.xvg".format(sys.argv[1]), comments=["#","@"])
density2_data = np.loadtxt("density2_{}.xvg".format(sys.argv[1]), comments=["#","@"])

conditions = [temp_data, None, None, temperature1_data, temperature2_data, \
              pressure_data1, pressure_data2, density1_data, density2_data]

# font dictionary for plots
font = {'color': 'black', 'weight': 'semibold', 'size': 18}

# make a figure with two subplots
fig, axes = plt.subplots(3,3, figsize=(16,16), constrained_layout=True)
for condition, ax in zip(conditions, axes.flat):
    if condition == None:
    continue

    # Calculate and format statistics
    average = np.round(np.mean(condition[:,1]), decimals=1)
    std = np.round(np.std(condition[:,1]), decimals=1)
    stats = '\n'.join((
    r'$\mu=%.1f$' % (average, ),
    r'$\sigma=%.1f$' % (std, )))

    # Make individual plot
    ax.plot(condition[:,0],condition[:,1])
    ax.set_xlabel("Time (ps)", fontdict=font, labelpad=5)
    ax.tick_params(axis='y', labelsize=14, direction='in', width=2, \
                    length=5, pad=10)
    ax.tick_params(axis='x', labelsize=14, direction='in', width=2, \
                    length=5, pad=10)
    for i in ["top","bottom","left","right"]:
        ax.spines[i].set_linewidth(2)
    ax.grid(True)

    # Include stat summary in a text box
    ax.text(0.80, 0.15, stats, color='black', va='top', ha="left", transform=ax.transAxes,
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1'))

# add y-axis names and save figures
for i in range(0,9,3):
    axes[i].set_ylabel("Temperature (K)", fontdict=font, labelpad=5)
    axes[i+1].set_ylabel("Pressure (bar)", fontdict=font, labelpad=5)
    axes[i+2].set_ylabel("Density (g/L)", fontdict=font, labelpad=5)
plt.savefig("equilibriation.png")
