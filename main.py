""" Main functions and run """


# librairies 


import os
from env import *
from unfolding import g_mult_unfolding
from plots import plot_g_mult, plot_angmom, plot_ab_fit, plot_gamma_energy
from calorimetry import compute_total_gamma_energy, GammaEnergyConfig

# run


if __name__ == "__main__":

    os.makedirs('outputs', exist_ok=True)
    
    # SCONE raw measurements at 56us 

    energies, g_mult_raw_56us, stat_err_56us = scone_meas(filename = "238U_meas_mg_56us.csv")
    
    # SCONE raw measurements at 5.6us 

    _, g_mult_raw_5us6, stat_err_5us6 = scone_meas(filename = "238U_meas_mg_5us6.csv")

    # merging @ 4 MeV

    merg = 3
    g_mult_raw = np.concatenate((g_mult_raw_56us[:merg], g_mult_raw_5us6[merg:30]))
    stat_err_raw = np.concatenate((stat_err_56us[:merg], stat_err_5us6[merg:30]))

    # A, B fit on Geant4 simulations of SCONE

    _ = plot_ab_fit(FILENAMES)

    # unfolding and saving

    g_mult, syst_err, stat_err = g_mult_unfolding(energies, g_mult_raw, stat_err=stat_err_raw, out_name="g_mult.csv")

    # plot gamma-rays multiplicity vs. incident energy

    _ = plot_g_mult(energies[1:], g_mult[1:], syst_err[1:], stat_err[1:])

    # plot gamma-rays multiplicity vs. angular momentum

    _ = plot_angmom(energies[1:], g_mult[1:], stat_err[1:])

    # gamma-rays calorimetry

    df_nubar = pd.read_csv(
        "./data/evaluations/238U_nubar_JEFF4.csv",
        sep=r"\s+|,|;",
        engine="python",
        comment="#",
        header=None,
    )

    df_nubar = df_nubar.iloc[:, :2].copy()
    df_nubar.columns = ["energy", "nubar"]
    df_nubar["energy"] = pd.to_numeric(df_nubar["energy"], errors="coerce")
    df_nubar["nubar"] = pd.to_numeric(df_nubar["nubar"], errors="coerce")
    df_nubar = df_nubar.dropna(subset=["energy", "nubar"])
    df_nubar = df_nubar.sort_values("energy")

    df_egamma = compute_total_gamma_energy(
        root_10ns=Path("./data/scone/Egtot_twoanodes_run19_20_21_22_24_dtfission_10ns_Q12.root"),
        root_30ns=Path("./data/scone/Egtot_twoanodes_run19_20_21_22_24_dtfission_30ns_Q12.root"),
        nubar_eval=df_nubar,
    )

    plot_gamma_energy(df_egamma)