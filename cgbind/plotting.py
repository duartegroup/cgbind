import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from .log import logger


mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['axes.linewidth'] = 1


def plot_heatmap(ligand_names, substrate_names, ba_array, units_kcal_mol=True, units_kj_mol=False):
    """
    For a list of lists (/matrix) of binding affinities in the form

    # q1     q2    q3    q4    q5    q6
    [
     [5.0, 5.9, 7.6, 9.5, 5.5, 4.7],      # L-1
     [4.9, 5.9, 6.0, 7.3, 3.9, 4.9],      # L-2
     [3.9, 4.7, 5.7, 6.6, 4.1, 3.0],      # L-3
     [1.6, 2.5, 3.4, 2.1, 3.0, 1.2],      # L-4
     [-1.5, -2.3, 1.4, -0.8, -4.2, -1.8]  # L-5
     ]


    :param ligand_names: List of ligand names
    :param substrate_names:
    :param ba_array:
    :param units_kcal_mol:
    :param units_kj_mol:
    :return:
    """
    logger.info('Generating heatmap with mpl and saving as binding_affinities.png')

    fig, ax = plt.subplots()
    im = ax.imshow(-ba_array, cmap='RdYlGn')        # Negative here so green is the most negative ∆E (strongest binder)

    ax.set_xticks(np.arange(len(substrate_names)))
    ax.set_xticklabels(substrate_names, fontsize='large', fontweight='bold')

    ax.set_yticks(np.arange(len(ligand_names)))
    ax.set_yticklabels(ligand_names, fontsize='large', fontweight='bold')

    plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")

    for i in range(len(ligand_names)):
        for j in range(len(substrate_names)):
            ax.text(j, i, np.round(ba_array[i, j], 1), ha="center", va="center", color="k")

    cb = fig.colorbar(im)
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_visible(False)

    if units_kcal_mol:
        cb.set_label('$\Delta E$ / kcal mol$^{-1}$')
    if units_kj_mol:
        cb.set_label('$\Delta E$ / kJ mol$^{-1}$')

    plt.tight_layout()
    plt.savefig('binding_affinities.png')

    return 0
