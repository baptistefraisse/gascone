from pathlib import Path
import numpy as np
import pandas as pd
from env import *
from unfolding import g_mult_unfolding, delta_gamma_uncertainties, ng_pileup
from plots import plot_g_mult, plot_angmom, plot_ab_fit, plot_gamma_energy, plot_spin_partition
from calorimetry import compute_total_gamma_energy


USE_AB_COVARIANCE = True   # True: full fit covariance; False: conservative
N_BINS = 30
MERGE_INDEX = 3             # 56 us: first 3 bins; 5.6 us: remaining bins
REFERENCE_ENERGY = 2.0     # MeV for Delta n_gamma
PILEUP_REL_ERR = 0.10

OUTPUT_DIR = Path(OUT_DIR)
REPORT_ENERGIES = (5, 10, 15, 20, 25, 30)


def load_measurements():
    """Load aligned energy grids and combine acquisition windows."""
    e56, raw56, err56 = scone_meas(filename="238U_meas_mg_56us.csv")
    e5, raw5, err5 = scone_meas(filename="238U_meas_mg_5us6.csv")
    e56, e5 = np.asarray(e56, float), np.asarray(e5, float)
    if min(len(e56), len(e5)) < N_BINS:
        raise ValueError(f"Fewer than {N_BINS} energy bins in SCONE data")
    np.testing.assert_allclose(e56[:N_BINS], e5[:N_BINS])
    raw = np.r_[np.asarray(raw56[:MERGE_INDEX], float),
                np.asarray(raw5[MERGE_INDEX:N_BINS], float)]
    stat = np.r_[np.asarray(err56[:MERGE_INDEX], float),
                 np.asarray(err5[MERGE_INDEX:N_BINS], float)]
    return e56[:N_BINS], raw, stat


def validate_fit_covariance():
    """Validate covariance and derive the selected matrix."""
    pcov = np.asarray(PCOV_AB, float)
    if pcov.shape != (2, 2) or not np.allclose(pcov, pcov.T):
        raise ValueError("PCOV_AB must be a symmetric 2x2 matrix")
    if np.linalg.eigvalsh(pcov).min() < -1e-12:
        raise ValueError("PCOV_AB must be positive semidefinite")
    np.testing.assert_allclose(np.diag(pcov), [DA**2, DB**2], rtol=1e-3)

    # Omitting the cross term MUST be applied to both the absolute
    # and reference-subtracted uncertainty propagations.
    selected = pcov.copy() if USE_AB_COVARIANCE else np.diag(np.diag(pcov))
    return pcov, selected


def absolute_budget(energies, raw, stat_raw, n_gamma, selected_pcov):
    """1-sigma budget of absolute unfolded n_gamma(E)."""
    nubar = np.asarray(nubar_jeff, float)
    dnubar = np.asarray(nubar_jeff_err, float)
    arrays = (raw, stat_raw, n_gamma, nubar, dnubar)
    if any(np.shape(x) != np.shape(energies) for x in arrays):
        raise ValueError("All inputs must match the energy grid")

    pileup = np.asarray(ng_pileup(energies), float)
    g = raw - C * nubar + pileup
    if np.any((g < 0) | (g >= A)):
        raise ValueError("Corrected detector multiplicity must satisfy 0 <= g < A")

    # n_gamma = -B log(1-g/A); partial derivatives at fixed pile-up mean.
    d_a = -B * g / (A * (A - g))
    d_b = -np.log1p(-g / A)
    d_g = B / (A - g)
    np.testing.assert_allclose(n_gamma, -B * np.log1p(-g / A), rtol=1e-9)

    grad = np.column_stack([d_a, d_b])
    var_ab = np.einsum('ni,ij,nj->n', grad, selected_pcov, grad)
    if np.any(var_ab < -1e-10):
        raise ValueError("Negative A-B variance: verify fit covariance")

    parts = {
        'stat': np.abs(d_g * stat_raw),
        'ab': np.sqrt(np.maximum(var_ab, 0)),
        'c': np.abs(d_g * nubar * DC),
        'nubar': np.abs(d_g * C * dnubar),
        'pileup': np.abs(d_g * PILEUP_REL_ERR * pileup),
    }
    parts['syst'] = np.sqrt(sum(parts[k]**2 for k in ('ab', 'c', 'nubar', 'pileup')))
    parts['total'] = np.hypot(parts['stat'], parts['syst'])
    np.testing.assert_allclose(
        parts['syst']**2,
        sum(parts[k]**2 for k in ('ab', 'c', 'nubar', 'pileup')),
        rtol=1e-12, atol=1e-12,
    )
    return parts


def delta_budget(energies, raw, stat_raw, selected_pcov):
    """Use the SAME selected A-B treatment for Delta n_gamma(E)."""
    pileup = np.asarray(ng_pileup(energies), float)
    g = raw - C * np.asarray(nubar_jeff, float) + pileup
    dstat, dsyst, dtotal, dab, dc, dnubar, dpileup = delta_gamma_uncertainties(
        energies=energies,
        g_mult=g,
        raw_stat_err=stat_raw,
        nubar=nubar_jeff,
        dnubar=nubar_jeff_err,
        A=A, B=B, C=C, DC=DC,
        pcov_AB=selected_pcov,  # critical: same switch as Fig. 3
        reference_energy=REFERENCE_ENERGY,
    )
    return dict(stat=dstat, syst=dsyst, total=dtotal,
                ab=dab, c=dc, nubar=dnubar, pileup=dpileup)


def save_results(energies, n_gamma, budget, delta):
    """Save selected treatment with unambiguous columns."""
    np.savetxt(
        OUTPUT_DIR / 'g_mult.csv',
        np.column_stack((energies, n_gamma, budget['total'])),
        fmt='%.8g', header='energy g_mult g_mult_err', comments='',
    )
    pd.DataFrame({
        'energy_MeV': energies,
        'n_gamma': n_gamma,
        'stat': budget['stat'], 'syst': budget['syst'],
        'total': budget['total'], 'syst_AB': budget['ab'],
        'syst_C': budget['c'], 'syst_nubar': budget['nubar'],
        'syst_pileup': budget['pileup'],
        'delta_stat': delta['stat'], 'delta_syst': delta['syst'],
        'delta_total': delta['total'], 'delta_syst_AB': delta['ab'],
    }).to_csv(OUTPUT_DIR / 'g_mult_uncertainty_budget.csv', index=False)
    (OUTPUT_DIR / 'ab_covariance_mode.txt').write_text(
        'with_covariance\n' if USE_AB_COVARIANCE else 'without_covariance\n'
    )


def plot_calorimetry():
    df = pd.read_csv(
        './data/evaluations/238U_nubar_JEFF4.csv',
        sep=r'\s+|,|;', engine='python', comment='#', header=None,
    ).iloc[:, :2].copy()
    df.columns = ['energy', 'nubar']
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['energy', 'nubar']).sort_values('energy')
    df_egamma = compute_total_gamma_energy(
        root_10ns=Path('./data/scone/Egtot_twoanodes_run19_20_21_22_24_dtfission_10ns_Q12.root'),
        root_30ns=Path('./data/scone/Egtot_twoanodes_run19_20_21_22_24_dtfission_30ns_Q12.root'),
        nubar_eval=df,
    )
    plot_gamma_energy(df_egamma)


def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    energies, raw, stat_raw = load_measurements()
    plot_ab_fit(FILENAMES)
    pcov, selected_pcov = validate_fit_covariance()

    n_gamma, legacy_total, _ = g_mult_unfolding(
        energies, raw, stat_err=stat_raw, out_name=None,
    )

    budget = absolute_budget(energies, raw, stat_raw, n_gamma, selected_pcov)
    delta = delta_budget(energies, raw, stat_raw, selected_pcov)
    save_results(energies, n_gamma, budget, delta)

    sl = slice(1, None)  # Original plots exclude 1 MeV.
    # Both panels of Fig. 3 use the SAME selected budget.
    plot_g_mult(
        energies[sl], n_gamma[sl], budget['syst'][sl],
        budget['stat'][sl], budget['ab'][sl], budget['c'][sl],
        budget['nubar'][sl], budget['pileup'][sl],
    )
    # Angular momentum: delta from 2 MeV, with the SAME covariance choice.
    plot_angmom(energies[sl], n_gamma[sl], delta['total'][sl])
    plot_calorimetry()
    plot_spin_partition(eq_number=14, save_path='figs/spin_partition_Eq14.pdf')


if __name__ == '__main__':
    main()