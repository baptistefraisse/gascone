""" Main functions and run """


# librairies 


import os
from env import *
from unfolding import g_mult_unfolding
from plots import plot_g_mult, plot_angmom, plot_ab_fit
from angmom import angmom_capture, angmom_emission, angmom_frag, electric_trans

# run


if __name__ == "__main__":

    os.makedirs('figs', exist_ok=True)
    
    # SCONE raw measurements at 56us 

    energies, g_mult_raw_56us = scone_meas(filename = "238U_meas_mg_56us.csv")
    
    # SCONE raw measurements at 5.6us 

    _, g_mult_raw_5us6 = scone_meas(filename = "238U_meas_mg_56us.csv")

    # merging

    g_mult_raw = np.concatenate((g_mult_raw_56us[:4], g_mult_raw_5us6[4:30]))

    # A, B fit

    _ = plot_ab_fit(FILENAMES)

    # unfolding

    g_mult, g_mult_err = g_mult_unfolding(energies, g_mult_raw)

    # plot gamma-rays multiplicity vs. incident energy

    _ = plot_g_mult(energies[1:], g_mult[1:], g_mult_err[1:])

    # plot gamma-rays multiplicity vs. angular momentum

    _ = plot_angmom(energies[1:], g_mult[1:], g_mult_err[1:])

