""" Unfolding functions """

# librairies

import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from utils import *
from plots import *

# functions 

def gamma_unfolding(a, b, g_mult):
    """
    SCONE response function to pure fission gamma-rays.

    Args:
        a (float): First gamma-rays multiplicity SCONE constant.
        b (float): Second gamma-rays multiplicity SCONE constant.
        g_mult (float): Raw measured gamma-rays multiplicity.

    Returns:
        g_mult_corr (float): Unfolded gamma-rays multiplicity.
    """

    # vectorization 

    g_mult = np.asarray(g_mult, dtype=float)

    # SCONE response function 

    g_mult_corr = - b * np.log(1 - g_mult/a) # a*(1.0-np.exp(-g_mult/b))

    return g_mult_corr


def gamma_unfolding_uq(a, da, b, db, g_mult, g_mult_err=None):
    """
    SCONE response function to pure fission gamma-rays 
    with uncertainty quantification.

    Args:
        a (float): First gamma-rays multiplicity SCONE constant.
        da (float): 1-sigma uncertainty of a. 
        b (float): Second gamma-rays multiplicity SCONE constant.
        db (float): 1-sigma uncertainty of b. 
        g_mult (float): Raw measured gamma-rays multiplicity.
        g_mult_err (float): Statistical error of g_mult.

    Returns:
        g_mult_corr (float): Unfolded gamma-rays multiplicity.
        sigma (float): Sigma error of m_corr.
    """

    # average response correction

    g_mult = np.asarray(g_mult, dtype=float)
    g_mult_corr = gamma_unfolding(a, b, g_mult)

    # error propagation with partial derivatives

    df_da = - b * (g_mult / (a**2 * (1.0 - g_mult / a)))
    df_db = - np.log(1.0 - g_mult / a)
    df_dg =   b / (a * (1.0 - g_mult / a))

    var = (df_da * da)**2 + (df_db * db)**2
    if g_mult_err is not None:
        g_mult_err = np.asarray(g_mult_err, dtype=float)
        var += (df_dg * g_mult_err)**2

    return g_mult_corr, np.sqrt(var)


def neutron_contamination(c, dc, nubar, dnubar):
    """
    Neutron contamination in SCONE assemblies.

    Args:
        c (float): Neutron contamination constant of SCONE.
        dc (float): 1-sigma error of c.
        nubar (float or narray): Average neutron multiplicity of fission.
        dnubar (float or narray): 1-sigma error of nubar.

    Returns:
        (float or narray): Neutron contamination to the SCONE multplicity.
    """
    
    n_contam = c * nubar
    sigma_n_contam = np.sqrt((nubar * dc)**2 + (c * dnubar)**2)
    return n_contam, sigma_n_contam


def ng_pileup(e):
    """
    Neutron-gamma pile-up on SCONE assemblies.

    Args:
        e (float or narray): Incident neutron energy [MeV].

    Returns:
        (float or narray): Interative result for neutron-gamma pile-up multiplicity (fit vs. En to simplify).
    """
    pileup = 0.15 + 0.4*(e-1)/30
    return pileup

def g_mult_unfolding(energies, g_mult_raw, stat_err=None, out_name=None):
    """
    Gamma-rays multiplicity unfolded from SCONE measurements.

    The pile-up uncertainty is assumed to be 10% of the pile-up
    correction at each incident neutron energy.

    Args:
        energies (narray): Incident neutron energies [MeV].
        g_mult_raw (narray): Raw measurements of gamma-rays multiplicity by SCONE.
        stat_err (narray): 1-sigma statistical error of g_mult_raw.
        out_name (str): Name of the output CSV file.

    Returns:
        (narray): Unfolded gamma-rays multiplicities.
        (narray): 1-sigma total error on unfolded gamma-rays multiplicities.
        (narray): 1-sigma propagated statistical error.
    """

    # ----------------------------------------------------------
    # Neutron correction
    # ----------------------------------------------------------

    n_contam, n_contam_err = neutron_contamination(
        C, DC, nubar_jeff, nubar_jeff_err
    )

    # ----------------------------------------------------------
    # Neutron-gamma pile-up correction
    # ----------------------------------------------------------

    pileup = ng_pileup(energies)

    # Assume a 10% relative uncertainty
    pileup_err = 0.10 * np.abs(pileup)

    # Corrected detector multiplicity before gamma unfolding
    g_mult = g_mult_raw - n_contam + pileup

    # ----------------------------------------------------------
    # Uncertainty propagation before nonlinear unfolding
    # ----------------------------------------------------------

    err_terms = [
        n_contam_err,
        pileup_err
    ]

    if stat_err is not None:
        stat_err = np.asarray(stat_err, dtype=float)
        err_terms.append(stat_err)

    g_mult_err = np.sqrt(
        np.sum([term**2 for term in err_terms], axis=0)
    )

    # ----------------------------------------------------------
    # Nonlinear gamma unfolding
    # ----------------------------------------------------------

    g_mult_corr, g_mult_corr_err = gamma_unfolding_uq(
        A, DA,
        B, DB,
        g_mult,
        g_mult_err
    )

    # Statistical uncertainty propagated through the response
    df_dg = B / (
        A * (1.0 - g_mult / A)
    )

    if stat_err is not None:
        stat_err_corr = np.abs(df_dg) * stat_err
    else:
        stat_err_corr = None

    # ----------------------------------------------------------
    # Saving
    # ----------------------------------------------------------

    if out_name is not None:
        data = np.column_stack(
            [energies, g_mult_corr, g_mult_corr_err]
        )

        np.savetxt(
            OUT_DIR / out_name,
            data,
            fmt="%.3g %.3g %.3g",
            header="energy g_mult g_mult_err",
            comments=""
        )

    return g_mult_corr, g_mult_corr_err, stat_err_corr

# def g_mult_unfolding(energies, g_mult_raw, stat_err=None, out_name=None):
#     """
#     Gamma-rays multiplicity unfolded from SCONE measurements.

#     Args:
#         energies (narray): Incident neutron energies [MeV].
#         g_mult_raw (narray): Raw measurements of gamma-rays nultiplicity by SCONE.
#         stat_err (narray): 1-sigma statistical error of g_mult_raw.
#         out_name (str): Name of the output CSV file.

#     Returns:
#         (narray): Unfolded gamma-rays multiplicities.
#         (narray): 1-sigma error on unfolded gamma-rays multiplicities.
#     """

#     # neutron corrections

#     n_contam, n_contam_err = neutron_contamination(C, DC, nubar_jeff, nubar_jeff_err)
#     pileup = ng_pileup(energies)
#     g_mult = g_mult_raw - n_contam + pileup

#     # uncertainty propagation

#     err_terms = [n_contam_err]

#     if stat_err is not None:
#         err_terms.append(np.asarray(stat_err, dtype=float))

#     if len(err_terms) == 1:
#         g_mult_err = err_terms[0]

#     else:
#         g_mult_err = np.sqrt(np.sum([t**2 for t in err_terms], axis=0))

#     # unfolding

#     g_mult_corr, g_mult_corr_err = gamma_unfolding_uq(A, DA, B, DB, g_mult, g_mult_err)
#     stat_err_corr = gamma_unfolding(A, B, stat_err)

#     # saving csv

#     if out_name is not None:
#         data = np.column_stack([energies, g_mult_corr, g_mult_corr_err])
#         np.savetxt(
#             OUT_DIR/out_name,
#             data,
#             fmt="%.3g %.3g %.3g",
#             header="energy g_mult g_mult_err",
#             comments=""
#         )

#     return g_mult_corr, g_mult_corr_err, stat_err_corr

def delta_gamma_uncertainties(
    energies,
    g_mult,
    raw_stat_err,
    nubar,
    dnubar,
    A,
    B,
    C, DC,
    pcov_AB,
    reference_energy=2.0
):
    energies = np.asarray(energies, dtype=float)
    g_mult = np.asarray(g_mult, dtype=float)
    raw_stat_err = np.asarray(raw_stat_err, dtype=float)
    nubar = np.asarray(nubar, dtype=float)
    dnubar = np.asarray(dnubar, dtype=float)

    # Find reference-energy index
    iref = np.argmin(np.abs(energies - reference_energy))

    if not np.isclose(energies[iref], reference_energy):
        raise ValueError(
            f"Reference energy {reference_energy} MeV not found "
            f"(closest value: {energies[iref]} MeV)"
        )

    # ==========================================================
    # Derivatives of unfolded multiplicity
    # n_gamma = -B log(1 - g/A)
    # ==========================================================

    df_da = -B * g_mult / (
        A**2 * (1.0 - g_mult / A)
    )

    df_db = -np.log(
        1.0 - g_mult / A
    )

    df_dg = B / (
        A * (1.0 - g_mult / A)
    )

    # ==========================================================
    # Statistical uncertainty
    # ==========================================================

    stat = np.abs(df_dg) * raw_stat_err

    # Statistics are independent between energy bins
    dstat = np.sqrt(
        stat**2 +
        stat[iref]**2
    )

    # ==========================================================
    # A-B systematic uncertainty WITH covariance
    # ==========================================================

    D_A = df_da - df_da[iref]
    D_B = df_db - df_db[iref]

    var_A = pcov_AB[0, 0]
    var_B = pcov_AB[1, 1]
    cov_AB = pcov_AB[0, 1]

    var_AB = (
        D_A**2 * var_A
        + D_B**2 * var_B
        + 2.0 * D_A * D_B * cov_AB
    )

    var_AB = np.maximum(var_AB, 0.0)

    syst_AB = np.sqrt(var_AB)

    # ==========================================================
    # C systematic uncertainty
    # ==========================================================

    # g = g_raw - C*nubar + pileup
    # dg/dC = -nubar

    deriv_C = -df_dg * nubar

    D_C = deriv_C - deriv_C[iref]

    syst_C = np.abs(D_C) * DC

    # ==========================================================
    # nubar uncertainty
    # ==========================================================

    # Assumption: nubar uncertainties independent between bins

    syst_nubar = np.sqrt(
        (df_dg * C * dnubar)**2
        +
        (df_dg[iref] * C * dnubar[iref])**2
    )

    # ==========================================================
    # Pile-up uncertainty
    # ==========================================================

    pileup = ng_pileup(energies)

    # 10% relative uncertainty on the pile-up correction
    pileup_err = 0.10 * np.abs(pileup)

    # Since dg/d(pileup) = +1:
    #
    # dn_gamma/d(pileup) = df/dg
    #
    # Here we assume pile-up uncertainties are independent
    # between incident-energy bins.

    syst_pileup = np.sqrt(
        (df_dg * pileup_err)**2
        +
        (df_dg[iref] * pileup_err[iref])**2
    )

    # ==========================================================
    # Total residual systematic
    # ==========================================================

    syst_total = np.sqrt(
        syst_AB**2
        + syst_C**2
        + syst_nubar**2
        + syst_pileup**2
    )

    # ==========================================================
    # Total uncertainty
    # ==========================================================

    total = np.sqrt(
        dstat**2 +
        syst_total**2
    )

    # Reference point is Delta n_gamma = 0 exactly
    dstat[iref] = 0.0
    syst_AB[iref] = 0.0
    syst_C[iref] = 0.0
    syst_nubar[iref] = 0.0
    syst_pileup[iref] = 0.0
    syst_total[iref] = 0.0
    total[iref] = 0.0

    return (
        dstat,
        syst_total,
        total,
        syst_AB,
        syst_C,
        syst_nubar,
        syst_pileup
    )

# def delta_gamma_uncertainties(
#     energies,
#     g_mult,
#     raw_stat_err,
#     nubar,
#     dnubar,
#     A,
#     B,
#     C, DC,
#     pcov_AB,
#     reference_energy=2.0
# ):
#     energies = np.asarray(energies, dtype=float)
#     g_mult = np.asarray(g_mult, dtype=float)
#     raw_stat_err = np.asarray(raw_stat_err, dtype=float)
#     nubar = np.asarray(nubar, dtype=float)
#     dnubar = np.asarray(dnubar, dtype=float)

#     # Find reference-energy index
#     iref = np.argmin(np.abs(energies - reference_energy))

#     if not np.isclose(energies[iref], reference_energy):
#         raise ValueError(
#             f"Reference energy {reference_energy} MeV not found "
#             f"(closest value: {energies[iref]} MeV)"
#         )

#     # ==========================================================
#     # Derivatives of unfolded multiplicity
#     # n_gamma = -B log(1 - g/A)
#     # ==========================================================

#     df_da = -B * g_mult / (
#         A**2 * (1.0 - g_mult / A)
#     )

#     df_db = -np.log(
#         1.0 - g_mult / A
#     )

#     df_dg = B / (
#         A * (1.0 - g_mult / A)
#     )

#     # ==========================================================
#     # Statistical uncertainty
#     # ==========================================================

#     stat = np.abs(df_dg) * raw_stat_err

#     # Statistics are independent between energy bins
#     dstat = np.sqrt(
#         stat**2 +
#         stat[iref]**2
#     )

#     # ==========================================================
#     # A-B systematic uncertainty WITH covariance
#     # ==========================================================

#     # Derivatives of Delta n_gamma
#     D_A = df_da - df_da[iref]
#     D_B = df_db - df_db[iref]

#     var_A = pcov_AB[0, 0]
#     var_B = pcov_AB[1, 1]
#     cov_AB = pcov_AB[0, 1]

#     var_AB = (
#         D_A**2 * var_A
#         + D_B**2 * var_B
#         + 2.0 * D_A * D_B * cov_AB
#     )

#     # Numerical protection against tiny negative values
#     var_AB = np.maximum(var_AB, 0.0)

#     syst_AB = np.sqrt(var_AB)

#     # ==========================================================
#     # C systematic uncertainty
#     # ==========================================================

#     # g = g_raw - C*nubar + pileup
#     # therefore dg/dC = -nubar

#     deriv_C = -df_dg * nubar

#     D_C = deriv_C - deriv_C[iref]

#     syst_C = np.abs(D_C) * DC

#     # ==========================================================
#     # nubar uncertainty
#     # ==========================================================

#     # CURRENT assumption:
#     # nubar uncertainties independent between energy bins

#     syst_nubar = np.sqrt(
#         (df_dg * C * dnubar)**2
#         +
#         (df_dg[iref] * C * dnubar[iref])**2
#     )

#     # ==========================================================
#     # Total residual systematic
#     # ==========================================================

#     syst_total = np.sqrt(
#         syst_AB**2
#         + syst_C**2
#         + syst_nubar**2
#     )

#     # ==========================================================
#     # Total uncertainty
#     # ==========================================================

#     total = np.sqrt(
#         dstat**2 +
#         syst_total**2
#     )

#     # Reference point is Delta n_gamma = 0 exactly
#     dstat[iref] = 0.0
#     syst_AB[iref] = 0.0
#     syst_C[iref] = 0.0
#     syst_nubar[iref] = 0.0
#     syst_total[iref] = 0.0
#     total[iref] = 0.0

#     return (
#         dstat,
#         syst_total,
#         total,
#         syst_AB,
#         syst_C,
#         syst_nubar
#     )