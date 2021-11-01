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

conditions = [temp_data, [None], [None], temperature1_data, pressure1_data, \
                density1_data, temperature2_data, pressure2_data, density2_data]

# font dictionary for plots
font = {'color': 'black', 'weight': 'semibold', 'size': 16}

# make a figure with two subplots
fig, axes = plt.subplots(3,3, figsize=(16,16), constrained_layout=True)
labels=[("NVT","red"),("NPT1","royalblue"),("NPT2","green")]
i=0
for condition, ax in zip(conditions, axes.flat):
    if None in condition:
        ax.axis("off")
        i += 1
        continue
    # Calculate and format statistics
    average = np.round(np.mean(condition[:,1]), decimals=1)
    std = np.round(np.std(condition[:,1]), decimals=1)
    stats = '\n'.join((
    r'$\mu=%.1f$' % (average, ),
    r'$\sigma=%.1f$' % (std, )))

    # Make individual plot
    ax.plot(condition[:,0],condition[:,1], label=labels[i//3][0], \
            color=labels[i//3][1])
    ax.set_xlabel("Time (ps)", fontdict=font, labelpad=5)
    ax.tick_params(axis='y', labelsize=14, direction='in', width=2, \
                    length=5, pad=10)
    ax.tick_params(axis='x', labelsize=14, direction='in', width=2, \
                    length=5, pad=10)
    for j in ["top","bottom","left","right"]:
        ax.spines[j].set_linewidth(2)
    ax.legend(loc=1)
    ax.grid(True)

    # Include stat summary in a text box
    ax.text(0.95, 0.10, stats, color='black', va='bottom', ha="right", transform=ax.transAxes,
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1'))
    i += 1

# add y-axis names and save figures
for i in range(0,3):
    axes[i,0].set_ylabel("Temperature (K)", fontdict=font, labelpad=5)
    if i > 0:
        axes[i,1].set_ylabel("Pressure (bar)", fontdict=font, labelpad=5)
        axes[i,2].set_ylabel("Density (g/L)", fontdict=font, labelpad=5)
plt.savefig("equilibriation.png")

