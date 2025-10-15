""" Plots functions """

# librairies

from env import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from angmom import angmom_capture, g_mult_electrans

# plots design

plt.rcParams.update({'font.size': FONT_SIZE, 'font.family': 'times new roman'})
mpl.rcParams['mathtext.fontset']='stix'
mpl.rcParams['font.family']='STIXGeneral'
mpl.rcParams['axes.linewidth'] = 2

# plot gamma-rays multiplicity

def plot_g_mult(energies, g_mult, g_mult_err):
    """
    Plotting gamma-rays multiplicity vs. incident energy.

    Args:
        energies (narray): Incident neutron energies [MeV].
        g_mult (narray): Unfolded gamma-rays multiplicity measurements by SCONE.
        g_mult_err (narray): 1-sigma error on unfolded gamma-rays multiplicities.

    Returns:
        (str): Plot file direction.
    """

    # figure design

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.tight_layout()
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    plt.xticks(fontsize=FONT_SIZE)
    plt.yticks(fontsize=FONT_SIZE)
    plt.ylim(4.8,9.9)
    plt.xlabel(r'Incident neutron energy: $E_\mathrm{n}$ (MeV)',size=FONT_SIZE)
    plt.ylabel(r'Average $\gamma$-rays multiplicity: $\bar{n}_\gamma$',size=FONT_SIZE)

    # scone 

    plt.errorbar(energies, g_mult, 
                 xerr = [0.5 for i in energies], yerr= g_mult_err,
                 fmt='.', color='red', linewidth=2, markersize=2,
                 label="SCONE (thres. 200 keV)")

    # literature

    plt.errorbar(laborie_energies, laborie_mult, 
                 xerr = laborie_energies_err, yerr=laborie_mult_err, 
                 fmt='s', color='black', linewidth=3, markersize=7,
                 label="Laborie (thres. 190 keV)")
    
    plt.errorbar(qi_energies, qi_mult,
                 yerr=qi_mult_err,
                 fmt='o', color='green', linewidth=3, markersize=7,
                 label="Qi (thres. 200 keV)")
    
    # simulations

    plt.plot(cgmf_energies, cgmf_mult,
             linestyle='-', linewidth=3, color='blue',
             label="CGMF (thres. 200 keV)")
    
    plt.plot(gef_energies, gef_mult,
             linestyle='--', linewidth=3, color='saddlebrown',
             label="GEF (thres. 200 keV)")
    
    # save
    
    plt.legend(loc = 'upper left',frameon=False, title_fontsize=FONT_SIZE)
    dir = FIG_DIR/'g_mult.png'
    plt.savefig(dir, bbox_inches='tight')
    return dir


# plot angular momentum


def plot_angmom(energies, g_mult, g_mult_err):
    """
    Plotting gamma-rays multiplicity difference vs. angular momentum difference.

    Args:
        energies (narray): Incident neutron energies [MeV].
        ln_diff (narray): Incident angular momentum difference from 1 MeV.
        g_mult_err (narray): 1-sigma error on unfolded gamma-rays multiplicities.

    Returns:
        (str): Plot file direction.
    """

    # figure design

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.tight_layout()
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel(r'$\Delta{}J_0$ ($\hbar$)',size=FONT_SIZE)
    plt.ylabel(r'$\Delta{}\bar{n}_\gamma$',size=FONT_SIZE)
    ax.set_ylim(-0.5,6.5)
    ax.set_yticks([0,1,2,3,4,5,6])

    # scone

    j0, j0_err = angmom_capture(energies)
    dj0 = diff_init(j0)
    dg_mult = diff_init(g_mult)

    plt.errorbar(dj0, dg_mult,
                 xerr = j0_err, yerr = g_mult_err,
                 fmt='.', color='red', label="SCONE (thres. 200 keV)")

    # L = 0 assumption

    ymin = g_mult_electrans(j0, s = SNU_MAX, frag_ratio = 1.0, pole = 2)
    ymax = g_mult_electrans(j0, s = SNU_MIN, frag_ratio = 1.0, pole = 1)
    dymin = diff_init(ymin)
    dymax = diff_init(ymax)

    plt.fill_between(dj0, dymin, dymax, 
                     alpha=0.3, color='saddlebrown', 
                     label=r'$\Delta{}L=0$')

    # L from microscopic calculations (G. Scamps)

    ymin = g_mult_electrans(j0, s = SNU_MAX, frag_ratio = 0.3, pole = 2)
    ymax = g_mult_electrans(j0, s = SNU_MIN, frag_ratio = 0.3, pole = 1)
    dymin = diff_init(ymin)
    dymax = diff_init(ymax)

    plt.fill_between(dj0, dymin, dymax, 
                     alpha=0.3, color='blue', 
                     label=r'$\Delta{}L_\mathrm{micro}$')

    # legend and annotations
    
    plt.legend(loc = 'upper left', frameon=False, title_fontsize=FONT_SIZE)

    plt.annotate(r'Pure E2 and 1$\hbar$ / prompt neutron',
                 ha='center', va='bottom', xytext=(3.3,0.45), xy=(3.3,0.45), rotation=21)
    plt.annotate(r'Pure E1 and 0.5$\hbar$ / prompt neutron',
                 ha='center', va='bottom', xytext=(2.5,0.85), xy=(2.5,0.85), rotation=48)

    # save and return 

    dir = FIG_DIR/'angmom.png'
    plt.savefig(dir, bbox_inches='tight')
    return dir