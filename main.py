""" Main functions and run """


# librairies 


import os
from env import *
from unfolding import g_mult_unfolding
from plots import plot_g_mult, plot_angmom, plot_ab_fit


# run


if __name__ == "__main__":

    os.makedirs('outputs', exist_ok=True)
    
    # SCONE raw measurements at 56us 

    energies, g_mult_raw_56us, stat_err_56us = scone_meas(filename = "238U_meas_mg_56us.csv")
    
    # SCONE raw measurements at 5.6us 

    _, g_mult_raw_5us6, stat_err_5us6 = scone_meas(filename = "238U_meas_mg_5us6.csv")

    # merging

    g_mult_raw = np.concatenate((g_mult_raw_56us[:4], g_mult_raw_5us6[4:30]))
    stat_err_raw = np.concatenate((stat_err_56us[:4], stat_err_5us6[4:30]))

    # A, B fit on Geant4 simulations of SCONE

    _ = plot_ab_fit(FILENAMES)

    # unfolding and saving

    g_mult, syst_err, stat_err = g_mult_unfolding(energies, g_mult_raw, stat_err=stat_err_raw, out_name="g_mult.csv")

    # plot gamma-rays multiplicity vs. incident energy

    _ = plot_g_mult(energies[1:], g_mult[1:], syst_err[1:], stat_err[1:])

    # plot gamma-rays multiplicity vs. angular momentum

    _ = plot_angmom(energies[1:], g_mult[1:], syst_err[1:])

