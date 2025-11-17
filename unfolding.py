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
        (float or narray): Fitted neutron-gamma pile-up multiplicity.
    """
    pileup = 0.15 + 0.4*(e-1)/30
    return pileup


def g_mult_unfolding(energies, g_mult_raw, stat_err=None, out_name=None):
    """
    Gamma-rays multiplicity unfolded from SCONE measurements.

    Args:
        energies (narray): Incident neutron energies [MeV].
        g_mult_raw (narray): Raw measurements of gamma-rays nultiplicity by SCONE.
        stat_err (narray): 1-sigma statistical error of g_mult_raw.
        out_name (str): Name of the output CSV file.

    Returns:
        (narray): Unfolded gamma-rays multiplicities.
        (narray): 1-sigma error on unfolded gamma-rays multiplicities.
    """

    # neutron corrections

    n_contam, n_contam_err = neutron_contamination(C, DC, nubar_jeff, nubar_jeff_err)
    pileup = ng_pileup(energies)
    g_mult = g_mult_raw - n_contam + pileup

    # uncertainty propagation

    err_terms = [n_contam_err]

    if stat_err is not None:
        err_terms.append(np.asarray(stat_err, dtype=float))

    if len(err_terms) == 1:
        g_mult_err = err_terms[0]

    else:
        g_mult_err = np.sqrt(np.sum([t**2 for t in err_terms], axis=0))

    # unfolding

    g_mult_corr, g_mult_corr_err = gamma_unfolding_uq(A, DA, B, DB, g_mult, g_mult_err)
    stat_err_corr = gamma_unfolding(A, B, stat_err)

    # saving csv

    if out_name is not None:
        data = np.column_stack([energies, g_mult_corr, g_mult_corr_err])
        np.savetxt(
            OUT_DIR/out_name,
            data,
            fmt="%.3g %.3g %.3g",
            header="energy g_mult g_mult_err",
            comments=""
        )

    return g_mult_corr, g_mult_corr_err, stat_err_corr
