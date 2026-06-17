""" Utils, paths and reading functions """


# librairies


import numpy as np
import pandas as pd
from pathlib import Path


# paths


PROJECT_DIR = Path('.') 
OUT_DIR = PROJECT_DIR / 'outputs'
FIG_DIR = PROJECT_DIR / 'figs'
DATA_DIR = PROJECT_DIR / 'data'
SCONE_DIR = DATA_DIR / 'scone'
GEANT4_DIR = DATA_DIR / 'geant4'
EVAL_DIR = DATA_DIR / 'evaluations'
LIT_DIR = DATA_DIR / 'literature'
SIMU_DIR = DATA_DIR / 'simulations'

# general functions 


def diff_init(x):
    """
    Compute the difference of an array from its first value.

    Args:
        x (narray): Array of values.

    Returns:
        narray: Array of values difference from the first one.
    """
    x0 = x[0]
    return x - x0


# CSV data reader 


def csv_data_reader(filepath, observable="mult", reverse_asym_errors=False):
    df = pd.read_csv(filepath, sep=r"\s+")

    energy_col = next((c for c in df.columns if c == "energy"), None)
    energy_err_col = next(
        (c for c in df.columns if c in ["energy_error", "energy_err", "err_nrj"]),
        None,
    )

    energies = df[energy_col].to_numpy() if energy_col else None
    energies_err = df[energy_err_col].to_numpy() if energy_err_col else None

    value_col = next((c for c in df.columns if c == observable), None)

    if value_col is None:
        raise ValueError(f"Column '{observable}' not found in {filepath}")

    err_sup_col = next(
        (
            c for c in df.columns
            if c in [
                f"{observable}_error_sup",
                f"{observable}_err_up",
                f"err_{observable}_up",
            ]
        ),
        None,
    )

    err_inf_col = next(
        (
            c for c in df.columns
            if c in [
                f"{observable}_error_inf",
                f"{observable}_err_low",
                f"err_{observable}_low",
            ]
        ),
        None,
    )

    err_col = next(
        (
            c for c in df.columns
            if c in [
                f"{observable}_error",
                f"{observable}_err",
                f"err_{observable}",
            ]
        ),
        None,
    )

    values = df[value_col].to_numpy()

    if err_sup_col and err_inf_col:
        if reverse_asym_errors:
            values_err = [
                df[err_sup_col].to_numpy(),
                df[err_inf_col].to_numpy(),
            ]
        else:
            values_err = [
                df[err_inf_col].to_numpy(),
                df[err_sup_col].to_numpy(),
            ]
    elif err_col:
        values_err = df[err_col].to_numpy()
    else:
        values_err = None

    return energies, energies_err, values, values_err


# SCONE measurements reader


def scone_meas(filename = "238U_meas_mg_56us.csv"):
    """
    Read raw gamma-rays multiplicity distribution by SCONE.

    Args:
        filename (str): Name of the file from low-level analysis.

    Returns:
        narray: Incident neutron energies.
        narray: Raw gamma-rays multiplicity measurements.
    """
    
    # reading 

    data = np.loadtxt(SCONE_DIR/filename)

    # components extraction 

    X = data[:,0]
    Y = data[:,1]
    Z = data[:,2]

    # pandas treatment

    df = pd.DataFrame({'X':X,'Y':Y,'Z':Z})
    Z_pivot = df.pivot(index='Y',columns='X',values='Z')
    Z_normalized = Z_pivot.div(Z_pivot.sum(axis=0),axis=1)
    X_grid, Y_grid = np.meshgrid(Z_normalized.columns,Z_normalized.index)
    X_grid, Y_grid = np.clip(X_grid,1,30), np.clip(Y_grid,0,20)
    Z_grid = Z_normalized.values
    Z_grid = np.clip(Z_grid,0,0.2)

    # simple values

    y_vals = Z_normalized.index.values[:, None]
    p = Z_normalized.values 

    # statistical moments

    mean_Y = (y_vals*p).sum(axis=0)/Z_normalized.sum(axis=0)
    mean_Y2 = ((y_vals**2) * p).sum(axis=0)/Z_normalized.sum(axis=0)
    var_Y = mean_Y2 - mean_Y**2
    var_Y = np.where(var_Y < 0, 0.0, var_Y) # protection against negative artefacts

    with np.errstate(divide='ignore', invalid='ignore'):
        sigma_mean = np.sqrt(var_Y / Z_pivot.sum(axis=0).values)
        sigma_mean[~np.isfinite(sigma_mean)] = np.nan

    # final outputs

    multg_raw = np.array(mean_Y)[1:31]
    multg_err = np.array(sigma_mean)[1:31]
    energies = np.unique(X)[1:31]

    return energies, multg_raw, multg_err
