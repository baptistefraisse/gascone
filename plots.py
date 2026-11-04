""" Plots functions """

# librairies

from env import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from angmom import angmom_capture, g_mult_electrans
from response import fission_gcasc_resp, g_mult_scone, fit_scone_gconst_multiple

# plots design

plt.rcParams.update({'font.size': FONT_SIZE, 'font.family': 'times new roman'})
mpl.rcParams['mathtext.fontset']='stix'
mpl.rcParams['font.family']='STIXGeneral'
mpl.rcParams['axes.linewidth'] = 2
plt.rcParams['pdf.fonttype'] = 42     # Texte en polices TrueType (pas converti en chemins)
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['lines.linewidth'] = 1.0 # Réduit l'épaisseur par défaut des courbes
plt.rcParams['axes.linewidth'] = 1.0  # Bordures d'axes plus fines

# plot gamma-rays multiplicity

def plot_g_mult(
    energies,
    g_mult,
    syst_err,
    stat_err,
    syst_AB,
    syst_C,
    syst_nubar,
    syst_pileup
):
    """
    Plot gamma-ray multiplicity vs incident neutron energy,
    with relative uncertainty contributions in a lower panel.
    """

    energies = np.asarray(energies, dtype=float)
    g_mult = np.asarray(g_mult, dtype=float)

    syst_err = np.asarray(syst_err, dtype=float)
    stat_err = np.asarray(stat_err, dtype=float)

    syst_AB = np.asarray(syst_AB, dtype=float)
    syst_C = np.asarray(syst_C, dtype=float)
    syst_nubar = np.asarray(syst_nubar, dtype=float)
    syst_pileup = np.asarray(syst_pileup, dtype=float)

    # ------------------------------------------------------------
    # Relative uncertainties [%]
    # ------------------------------------------------------------

    rel_stat = 100.0 * stat_err / g_mult

    rel_AB = 100.0 * syst_AB / g_mult
    rel_C = 100.0 * syst_C / g_mult
    rel_nubar = 100.0 * syst_nubar / g_mult
    rel_pileup = 100.0 * syst_pileup / g_mult

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------

    fig, (ax, ax_unc) = plt.subplots(
        2, 1,
        figsize=(10, 12),
        sharex=True,
        gridspec_kw={
            "height_ratios": [3, 1],
            "hspace": 0.05
        }
    )

    # ============================================================
    # TOP PANEL
    # ============================================================

    ax.set_ylim(4.8, 9.2)
    ax.set_yticks([5, 6, 7, 8, 9])

    ax.set_ylabel(
        r'Average $\gamma$-ray multiplicity $\bar{n}_\gamma$',
        size=FONT_SIZE
    )

    ax.tick_params(
        axis='both',
        which='major',
        size=10,
        width=2,
        direction='in',
        top=True,
        right=True,
        labelsize=FONT_SIZE
    )

    ax.tick_params(
        axis='both',
        which='minor',
        size=7,
        width=2,
        direction='in',
        top=True,
        right=True
    )

    # ------------------------------------------------------------
    # SCONE
    # ------------------------------------------------------------

    ax.errorbar(
        energies,
        g_mult,
        xerr=np.full_like(energies, 0.5),
        yerr=stat_err,
        fmt='.',
        color='red',
        linewidth=3,
        markersize=7,
        linestyle='none',
        label="SCONE (thresh. 200 keV)"
    )

    ax.fill_between(
        energies,
        g_mult - syst_err,
        g_mult + syst_err,
        color='red',
        alpha=0.15
    )

    # ------------------------------------------------------------
    # Literature
    # ------------------------------------------------------------

    ax.errorbar(
        laborie_energies,
        laborie_mult,
        xerr=laborie_energies_err,
        yerr=laborie_mult_err,
        fmt='s',
        color='black',
        linewidth=3,
        markersize=7,
        label="Laborie (thresh. 190 keV)"
    )

    ax.errorbar(
        qi_energies,
        qi_mult,
        yerr=qi_mult_err,
        fmt='o',
        color='green',
        linewidth=3,
        markersize=7,
        label="Qi (thresh. 200 keV)"
    )

    # ------------------------------------------------------------
    # Simulations
    # ------------------------------------------------------------

    ax.plot(
        cgmf_energies,
        cgmf_mult,
        linestyle='-',
        linewidth=3,
        color='blue',
        label="CGMF (thresh. 200 keV)"
    )

    ax.plot(
        gef_energies,
        gef_mult,
        linestyle='--',
        linewidth=3,
        color='saddlebrown',
        label="GEF (thresh. 200 keV)"
    )

    ax.legend(
        loc='upper left',
        frameon=False,
        fontsize=FONT_SIZE * 0.8
    )

    # ============================================================
    # BOTTOM PANEL: RELATIVE UNCERTAINTY CONTRIBUTIONS
    # ============================================================

    baseline = np.zeros_like(energies)

    ax_unc.fill_between(
        energies,
        -rel_AB,
        +rel_AB,
        color='#0072B2',
        alpha=0.2,
        linewidth=0,
        label=r'$\sigma_{A,B}$',
        zorder=1
    )

    ax_unc.fill_between(
        energies,
        -rel_C,
        +rel_C,
        color='#009E73',
        alpha=0.25,
        linewidth=0,
        label=r'$\sigma_C$',
        zorder=2
    )

    ax_unc.fill_between(
        energies,
        -rel_pileup,
        +rel_pileup,
        color='gray',
        alpha=0.3,
        linewidth=0,
        label=r'$\sigma_{\mathrm{pileup}}$',
        zorder=4
    )

    ax_unc.fill_between(
        energies,
        -rel_nubar,
        +rel_nubar,
        color='black',
        facecolor='none',
        edgecolor='black',
        #alpha=0.35,
        hatch='...',
        linewidth=0,
        label=r'$\sigma_{\bar{\nu}}$',
        zorder=3
    )

    ax_unc.errorbar(
        energies,
        baseline,
        xerr=np.full_like(energies, 0.),
        yerr=rel_stat,
        fmt='.',
        color='red',
        elinewidth=3,
        markersize=7,
        linestyle='none',
        label=r'$\sigma_{\mathrm{stat}}$',
        zorder=4
    )

    # ------------------------------------------------------------
    # Reference line
    # ------------------------------------------------------------

    ax_unc.axhline(
        0.0,
        color='black',
        linewidth=1.0,
        linestyle='--',
        zorder=0
    )

    # ------------------------------------------------------------
    # Axes
    # ------------------------------------------------------------

    ax_unc.set_ylabel(
        r'Uncertainties (%)',
        size=FONT_SIZE
    )

    ax_unc.set_xlabel(
        r'Incident neutron energy $E_\mathrm{n}$ (MeV)',
        size=FONT_SIZE
    )

    ax_unc.set_xticks([0, 10, 20, 30])

    ax_unc.tick_params(
        axis='both',
        which='major',
        size=10,
        width=2,
        direction='in',
        top=True,
        right=True,
        labelsize=FONT_SIZE
    )

    ax_unc.tick_params(
        axis='both',
        which='minor',
        size=7,
        width=2,
        direction='in',
        top=True,
        right=True
    )

    # ------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------

    ax_unc.legend(
        loc='upper center',
        frameon=False,
        fontsize=FONT_SIZE * 0.8,
        ncol=5,
        columnspacing=0.45,
        handletextpad=0.35,
        handlelength=1.2,
        borderaxespad=0.25
    )

    ax_unc.set_ylim(0,2.25)
    ax_unc.set_yticks([0,1,2])

    # ------------------------------------------------------------
    # Save
    # ------------------------------------------------------------

    dir = FIG_DIR / 'g_mult.pdf'

    plt.savefig(
        dir,
        format="pdf",
        bbox_inches='tight',
        dpi=300
    )

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
    plt.xticks(fontsize=FONT_SIZE)
    plt.yticks(fontsize=FONT_SIZE)
    plt.xlabel(r'$\Delta{}\bar{J}_0$ ($\hbar$)',size=FONT_SIZE)
    plt.ylabel(r'$\Delta{}\bar{n}_\gamma$',size=FONT_SIZE)
    ax.set_ylim(-0.5,6.5)
    ax.set_yticks([0,1,2,3,4,5,6])

    # SCONE

    j0, j0_err = angmom_capture(energies)
    dj0 = diff_init(j0)
    dg_mult = diff_init(g_mult)
    j0_err[0] = 0.0  # no error on reference energy
    
    # dg_mult_err = np.sqrt(np.array(g_mult_err)**2 + g_mult_err[0]**2)

    # plt.errorbar(dj0, dg_mult,
    #              xerr = j0_err, yerr = dg_mult_err,
    #              fmt='.', color='red', linewidth=3, markersize=7, linestyle='none',
    #              label="SCONE (thresh. 200 keV)")

    plt.errorbar(
        dj0,
        dg_mult,
        xerr=j0_err,
        yerr=g_mult_err,
        fmt='.',
        color='red',
        linewidth=3,
        markersize=7,
        linestyle='none',
        label="SCONE (thresh. 200 keV)"
    )
    
    # Literature data: Qi

    j0_qi, j0_qi_err = angmom_capture(qi_energies)
    dj0_qi = diff_init(j0_qi)
    dg_qi = diff_init(qi_mult)

    plt.errorbar(
        dj0_qi[1:],
        dg_qi[1:],
        xerr=None, # no energy error provided in Qi data
        yerr=[0.12], # stat only for this plot: difference study eliminates syst uncertainties
        fmt='o',
        color='green',
        linewidth=3,
        markersize=7,
        linestyle='none',
        label=r"Qi (thresh. 200 keV)"
    )

    # Literature data: Laborie

    j0_laborie, j0_laborie_err = angmom_capture(laborie_energies, de=[0.1,0.2,0.6])
    dj0_laborie = diff_init(j0_laborie)
    dj0_laborie_err = np.sqrt(j0_laborie_err**2 + j0_laborie_err[0]**2)
    dg_laborie = diff_init(laborie_mult)

    plt.errorbar(
        dj0_laborie[1:],
        dg_laborie[1:], 
        xerr=dj0_laborie_err[1:], # energy error laborie
        yerr=[0.12, 0.12], # stat only for this plot: difference study eliminates syst uncertainties
        fmt='s',
        color='black',
        linewidth=3,
        markersize=7,
        linestyle='none',
        label=r"Laborie (thresh. 190 keV)"
    )

    # L = 0 assumption

    ymin = g_mult_electrans(j0, s = SNU_MAX, frag_ratio = 1.0, pole = 2)
    ymax = g_mult_electrans(j0, s = SNU_MIN, frag_ratio = 1.0, pole = 1)
    dymin = diff_init(ymin)
    dymax = diff_init(ymax)

    # Light colored background
    plt.fill_between(
        dj0, dymin, dymax,
        color='saddlebrown',
        alpha=0.10,
        linewidth=0
    )

    # Hatching on top
    plt.fill_between(
        dj0, dymin, dymax,
        facecolor='none',
        hatch='...',
        edgecolor='saddlebrown',
        linewidth=0.,
        label='Full transfer to fragment spins'
    )

    # L from microscopic calculations (G. Scamps)

    ymin = g_mult_electrans(j0, s = SNU_MAX, frag_ratio = FRAG_MOM_MICRO, pole = 2)
    ymax = g_mult_electrans(j0, s = SNU_MIN, frag_ratio = FRAG_MOM_MICRO, pole = 1)
    dymin = diff_init(ymin)
    dymax = diff_init(ymax)

    plt.fill_between(dj0, dymin, dymax, 
                     alpha=0.3, color='blue', 
                     linewidth=0,
                     label='Microscopic partition\n(TDHF+BCS)')

    # legend and annotations
    
    plt.legend(loc = 'upper left', frameon=False, title_fontsize=FONT_SIZE)

    plt.annotate(r'Pure E2 and 0.5 $\hbar$ / prompt neutron',
                 ha='center', va='bottom', xytext=(4.,0.85), xy=(3.2,0.57), rotation=20.5, size=25)
    plt.annotate(r'Pure E1 and 0 $\hbar$ / prompt neutron',
                 ha='center', va='bottom', xytext=(3.06,1.3), xy=(2.45,0.85), rotation=47, size=25)

    # save and return 

    dir = FIG_DIR/'angmom.pdf'
    plt.savefig(dir, format="pdf", bbox_inches='tight', dpi=300)
    return dir


def parse_filename_label(filename):
    """
    Parse filenames like:
      - "Geant4_238U_GEF_10MeV.txt"    → "238U GEF (10 MeV)"
      - "Geant4_235U_CGMF_1MeV.txt"    → "235U CGMF (1 MeV)"
      - "Geant4_252Cf_FIFRELIN.txt"    → "252Cf FIFRELIN (spontaneous)"
      - "Geant4_252Cf.txt"             → "252Cf (spontaneous)"
    """
    base = filename.replace(".txt", "").replace("Geant4_", "")
    parts = base.split("_")

    # Cas 1 : seulement le noyau
    if len(parts) == 1:
        nucleus = parts[0]
        label = f"{nucleus} (spontaneous)"

    # Cas 2 : noyau + énergie
    elif len(parts) == 2 and parts[1].endswith("MeV"):
        nucleus, energy_raw = parts
        energy = energy_raw.replace("MeV", " MeV")
        label = f"{nucleus} ({energy})"

    # Cas 3 : noyau + modèle
    elif len(parts) == 2:
        model, nucleus = parts
        label = f"{nucleus} {model} (spontaneous)"

    # Cas 4 : noyau + modèle + énergie
    elif len(parts) == 3:
        model, nucleus, energy_raw = parts
        energy = energy_raw.replace("MeV", " MeV")
        label = f"{nucleus} {model} ({energy})"

    else:
        raise ValueError(f"Unexpected filename format: {filename}")

    return label



def plot_ab_fit(filenames, mult_range=None):
    """
    Plot the fit of A and B SCONE constants from multiple Geant4 simulations.
    """

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.tight_layout()
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    plt.xticks(fontsize=FONT_SIZE)
    plt.yticks(fontsize=FONT_SIZE)
    plt.xlabel(r"Emitted $\gamma$-rays multiplicity", size=FONT_SIZE)
    plt.ylabel("Average detected multiplicity (SCONE)", size=FONT_SIZE)
    ax.set_xlim(0.1, 17.5)
    ax.set_ylim(0.1, 13.5)
    ax.set_xticks([0, 4, 8, 12, 16])
    ax.set_yticks([2, 4, 6, 8, 10, 12])

    colors = plt.cm.tab10(np.linspace(0, 1, len(filenames)))
    markers = ["o", "s", "D", "^", "v", ">", "<", "p", "P", "X"]

    # Geant4 data
    handles_data, labels_data = [], []
    for i, fname in enumerate(filenames):
        label = parse_filename_label(fname)
        x, y = fission_gcasc_resp(fname, mult_range=mult_range)
        marker = markers[i % len(markers)]
        h = ax.scatter(x, y, color=colors[i], s=100, marker=marker, label=label)
        handles_data.append(h)
        labels_data.append(label)

    # Global fit
    a, b, da, db, _ = fit_scone_gconst_multiple(FILENAMES)
    print(a, b, da, db)
    xfit = np.linspace(0, 1.1 * max(x), 400)
    yfit = g_mult_scone(a, b, xfit)
    yfit_min = g_mult_scone(a - 3*da, b + 3*db, xfit)
    yfit_max = g_mult_scone(a + 3*da, b - 3*db, xfit)

    line_fit, = ax.plot(xfit, yfit, color="black", linewidth=3, label="Average fit")
    band_fit = ax.fill_between(xfit, yfit_min, yfit_max, color="gray", alpha=0.3, label=r"Fit uncertainty (3$\sigma$)")

    # --- Two separate legends ---
    order = sorted(range(len(labels_data)), key=lambda i: 0 if "252Cf" in labels_data[i] else 1)
    handles_data = [handles_data[i] for i in order]
    labels_data = [labels_data[i] for i in order]
    legend1 = ax.legend(handles=handles_data, labels=labels_data, loc='upper left', frameon=False)
    ax.add_artist(legend1)

    handles_fit = [line_fit, band_fit]
    labels_fit = ["Average fit", r"Fit uncertainty (3$\sigma$)"]
    ax.legend(handles_fit, labels_fit, loc='lower right', frameon=False)

    savename = "AB_fit.pdf"
    plt.savefig(FIG_DIR / savename, format="pdf", bbox_inches='tight', dpi=300)
    return savename


def plot_gamma_energy(df_egamma):
    """
    Plot total prompt gamma energy vs incident neutron energy.

    Args:
        df_egamma (pd.DataFrame): output of compute_total_gamma_energy()

    Returns:
        (str): Plot file direction.
    """

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.tight_layout()
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    plt.xticks(fontsize=FONT_SIZE)
    plt.yticks(fontsize=FONT_SIZE)
    plt.xlabel(r'Incident neutron energy: $E_\mathrm{n}$ (MeV)', size=FONT_SIZE)
    plt.ylabel(r'Total $\gamma$-ray energy: $E_{\gamma,\mathrm{tot}}$ (MeV)', size=FONT_SIZE)
    ax.set_ylim(4.5, 12.5)
    # ax.set_yticks([2, 4, 6, 8, 10, 12])
    ax.set_xlim(-1, 32)
    ax.set_xticks([0, 10, 20, 30])

    energies = df_egamma["energy"].to_numpy()
    egamma = df_egamma["egamma_tot"].to_numpy()
    stat_err = df_egamma["egamma_tot_stat"].to_numpy()
    sys_err = df_egamma["egamma_tot_sys"].to_numpy()

    # SCONE

    plt.errorbar(
        energies[1:31],
        egamma[1:31],
        xerr=[0.5 for _ in energies[1:31]],
        yerr=stat_err[1:31],
        fmt='.',
        color='red',
        linewidth=3,
        markersize=7,
        linestyle='none',
        label="SCONE"
    )

    # SCONE systematic uncertainty
    plt.fill_between(
        energies[1:31],
        egamma[1:31] - sys_err[1:31],
        egamma[1:31] + sys_err[1:31],
        color='red',
        alpha=0.15,
    )

    # Laborie

    plt.errorbar(
        laborie_energies,
        laborie_egtot,
        xerr=laborie_energies_err,
        yerr=laborie_egtot_err,
        fmt='s',
        color='black',
        linewidth=3,
        markersize=7,
        label="Laborie (thresh. 190 keV)"
    )

    # Qi

    plt.errorbar(
        qi_energies,
        qi_egtot,
        yerr=qi_egtot_err,
        fmt='o',
        color='green',
        linewidth=3,
        markersize=7,
        label="Qi (thresh. 200 keV)"
    )

    # simulations

    plt.plot(
        cgmf_energies_egtot,
        cgmf_egtot,
        linestyle='-',
        linewidth=3,
        color='blue',
        label="CGMF"
    )

    plt.plot(
        gef_energies_egtot,
        gef_egtot,
        linestyle='--',
        linewidth=3,
        color='saddlebrown',
        label="GEF"
    )

    plt.legend(loc='upper left', frameon=False, title_fontsize=FONT_SIZE)

    dir = FIG_DIR / 'egtot.pdf'
    plt.savefig(dir, format="pdf", bbox_inches='tight', dpi=300)
    return dir

def plot_spin_partition(eq_number=14, save_path='spin_partition.pdf'):

    x = np.arange(9, dtype=float)

    # Scission - dashed curves
    scission_l = np.array([
        0.000, 0.184, 0.368, 0.552, 0.736,
        0.920, 1.104, 1.288, 1.472
    ])

    scission_h = np.array([
        0.000, 0.041, 0.082, 0.123, 0.164,
        0.205, 0.246, 0.287, 0.328
    ])

    # Final - solid curves
    final_l = np.array([
        0.000, 0.253, 0.506, 0.759, 1.012,
        1.265, 1.518, 1.771, 2.024
    ])

    final_h = np.array([
        0.000, 0.037, 0.074, 0.111, 0.148,
        0.185, 0.222, 0.259, 0.296
    ])

    # Eq.(8) - green crosses
    eq8_l = np.array([
        0.000, 0.181, 0.364, 0.546, 0.728,
        0.916, 1.102, 1.283, 1.472
    ])

    # Eq.(8) - red circles
    eq8_h = np.array([
        0.000, 0.035, 0.073, 0.111, 0.151,
        0.191, 0.230, 0.268, 0.302
    ])

    fig, ax = plt.subplots(figsize=(8,4.5))

    green = "green"
    red = "red"

    # Four curves
    ax.plot(x, scission_l, "--", color=green,
            lw=2, label=r"Scission:  $\Delta J_l$")

    ax.plot(x, scission_h, "--", color=red,
            lw=2, label=r"Scission:  $\Delta J_h$")

    ax.plot(x, final_l, "-", color=green,
            lw=2, label=r"Final:  $\Delta J_l$")

    ax.plot(x, final_h, "-", color=red,
            lw=2, label=r"Final:  $\Delta J_h$")

    # Two scatter series
    ax.scatter(x, eq8_l, marker="+", color=green,
               s=55, linewidths=1.5, zorder=5,
               label=rf"Eq.({eq_number}):  $\Delta J_l$")

    ax.scatter(x, eq8_h, marker="o", color=red,
               s=55, zorder=6,
               label=rf"Eq.({eq_number}):  $\Delta J_h$")

    # Axes
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 2.5)

    ax.set_xticks(np.arange(0, 9, 1))
    ax.set_yticks(np.arange(0, 2.51, 0.5))

    ax.set_xlabel(r"$\Delta J_0$ ($\hbar$)", fontsize=18)
    ax.set_ylabel(r"$\Delta J_{h,l}$ ($\hbar$)", fontsize=18)

    ax.tick_params(axis="both", labelsize=14,
                   direction="out", length=5)

    ax.legend(loc="upper left", fontsize=13,
              frameon=False, title_fontsize=FONT_SIZE)

    plt.tight_layout()
    plt.savefig(save_path, format="pdf", bbox_inches='tight', dpi=300)
