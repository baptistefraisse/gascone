""" SCONE response analysis """


# imports


from env import *
from scipy.optimize import fsolve, curve_fit


# functions


def g_mult_scone(a, b, g_mult_casc):
    """
    SCONE response function  to gamma-rays.

    Args:
        a (float): First SCONE gamma-rays constant.
        b (float): Sedond SCONE gamma-rays constant.
        g_mult_casc (narray): Gamma-rays multiplicity of a fission cascade.

    Returns:
        g_mult_scone (narray): Gamma-rays multiplicity measured by SCONE.
    """
    g_mult_scone = a*(1.0-np.exp(-g_mult_casc/b))
    return g_mult_scone


def g_mult_scone_inv(a, b, g_mult_scone):
    """
    Inverse of SCONE response function to gamma-rays.

    Args:
        a (float): First SCONE gamma-rays constant.
        b (float): Sedond SCONE gamma-rays constant.
        g_mult_scone (narray): Gamma-rays multiplicity measured by SCONE.

    Returns:
        g_mult_casc (narray): Gamma-rays multiplicity of a fission cascade.
    """
    g_mult_casc = fsolve(lambda x: g_mult_scone(a,b,x)-g_mult_scone,1.)
    return g_mult_casc


def fission_gcasc_resp(filename, mult_range=None):
    """
    Extract the average SCONE response to fission gamma-rays from GEANT4 simulations.

    Args:
        filename (str): _description_. Defaults to "data/geant4/Geant4_238U_1MeV.txt".
        mult_range (tuple or None): (min, max) range of emitted multiplicity to keep. Defaults to None.

    Returns:
        emitted_mult (narray): Mean emitted gamma-rays multiplicity by fission.
        detected_mult (narray): Mean detected gamma-rays multiplicity by SCONE.
    """
    # reading 

    simu = np.loadtxt(GEANT4_DIR/filename)

    # extracting full response

    X_simu = simu[:,0]
    Y_simu = simu[:,1]
    Z_simu = simu[:,2]
    df_simu = pd.DataFrame({'X':X_simu,'Y':Y_simu,'Z':Z_simu})
    Z_pivot = df_simu.pivot(index='Y',columns='X',values='Z')

    # delete empty columns

    col_sums = Z_pivot.sum(axis=0).to_numpy()
    good_cols = col_sums > 0
    Z_pivot = Z_pivot.loc[:, good_cols]

    # normalization

    Z_norm = Z_pivot.div(Z_pivot.sum(axis=0), axis=1)

    # average response

    y_vals = Z_norm.index.values[:, None]
    mean_Y = (y_vals * Z_norm.values).sum(axis=0) / Z_norm.values.sum(axis=0)
    emitted_mult = Z_norm.columns.values
    detected_mult = np.asarray(mean_Y)

    # security 

    finite = np.isfinite(emitted_mult) & np.isfinite(detected_mult)
    emitted_mult = emitted_mult[finite]
    detected_mult = detected_mult[finite]

    # cut extreme multiplicity values (anecdotic fission events)
    
    if mult_range is not None:
        low, high = mult_range
        mask = (emitted_mult >= low) & (emitted_mult <= high)
        emitted_mult = emitted_mult[mask]
        detected_mult = detected_mult[mask]

    return emitted_mult, detected_mult


def fit_scone_gconst(emitted_mult, detected_mult):
    """
    Extract the average SCONE response to fission gamma-rays from GEANT4 simulations.

    Args:
        emitted_mult (narray): Mean emitted gamma-rays multiplicity by fission.
        detected_mult (narray): Mean detected gamma-rays multiplicity by SCONE.

    Returns:
        a (float): First SCONE gamma-ray constant.
        b (float): Second SCONE gamma-ray constant.
        a_err (float): 1-sigma uncertainty on a.
        b_err (float): 1-sigma uncertainty on b.
    """

    popt, pcov = curve_fit(
        lambda x, a, b: g_mult_scone(a, b, x),
        emitted_mult,
        detected_mult,
        p0=[10, 20], # initial trial
        bounds=(0, np.inf) # boundaries
    )

    # extraction and error

    a, b = popt
    a_err, b_err = np.sqrt(np.diag(pcov))
    
    return a, b, a_err, b_err



def fit_scone_gconst_multiple(filenames, mult_range=None):
    """
    Fit the SCONE gamma-ray response constants A and B for several cascades.

    Args:
        filenames (list of str): List of Geant4 simulation filenames.
        mult_range (tuple or None): (min, max) range of emitted multiplicity to keep. Defaults to None.

    Returns:
        a (float): Mean fitted a constant.
        b (float): Mean fitted b constant.
        a_err (float): Standard deviation of fitted a.
        b_err (float): Standard deviation of fitted b.
    """

    xs, ys = [], []

    for fname in filenames:
        x, y = fission_gcasc_resp(fname, mult_range=mult_range)
        finite = np.isfinite(x) & np.isfinite(y)
        xs.append(x[finite])
        ys.append(y[finite])

    X = np.concatenate(xs)
    Y = np.concatenate(ys)

    popt, pcov = curve_fit(
        lambda t, A, B: g_mult_scone(A, B, t),
        X, Y,
        p0=[20, 30],
        bounds=(0, np.inf),
        maxfev=20000
    )
    a, b = popt
    a_err, b_err = np.sqrt(np.diag(pcov))

    return a, b, a_err, b_err
