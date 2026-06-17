from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import uproot


@dataclass
class GammaEnergyConfig:
    eta: float = 0.55
    alpha: float = 8.0

    rho10: float = 0.35
    rho30: float = 0.75
    kappa10: float = -577e-3
    kappa30: float = -653e-3

    sigma_eta: float = 0.0
    sigma_rho10: float = 0.02
    sigma_rho30: float = 0.05
    sigma_kappa10: float = 64e-3
    sigma_kappa30: float = 138e-3
    sigma_nubar: float = 0.05


def _profile_y_mean_and_std(root_path: Path, hist_name: str = "egtot"):
    """
    Equivalent ROOT:
        h2->ProjectionY("proj", i, i)->GetMean()
        h2->ProjectionY("proj", i, i)->GetStdDev()

    Returns one y-mean and y-std per x-bin.
    """
    with uproot.open(root_path) as f:
        h = f[hist_name]
        values, x_edges, y_edges = h.to_numpy()

    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    means = np.full(values.shape[0], np.nan)
    stds = np.full(values.shape[0], np.nan)
    counts = values.sum(axis=1)

    for i in range(values.shape[0]):
        w = values[i, :]
        norm = w.sum()

        if norm <= 0:
            continue

        mean = np.sum(w * y_centers) / norm
        var = np.sum(w * (y_centers - mean) ** 2) / norm

        means[i] = mean
        stds[i] = np.sqrt(var)

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

    return x_centers, means, stds, counts


def compute_total_gamma_energy(
    root_10ns: Path,
    root_30ns: Path,
    nubar_eval: pd.DataFrame,
    cfg: GammaEnergyConfig | None = None,
    hist_name: str = "egtot",
) -> pd.DataFrame:
    """
    Compute SCONE total gamma energy from the two coincidence windows.

    nubar_eval must contain:
        energy, nubar

    Returns:
        energy, egamma_tot, egamma_tot_stat, egamma_tot_sys,
        egamma_tot_min, egamma_tot_max
    """
    if cfg is None:
        cfg = GammaEnergyConfig()

    energy, mean_10, std_10, counts_10 = _profile_y_mean_and_std(root_10ns, hist_name)
    energy_30, mean_30, std_30, counts_30 = _profile_y_mean_and_std(root_30ns, hist_name)

    if not np.allclose(energy, energy_30):
        raise ValueError("The 10 ns and 30 ns histograms do not have the same x-binning.")

    nubar_interp = np.interp(
        energy,
        nubar_eval["energy"].to_numpy(),
        nubar_eval["nubar"].to_numpy(),
    )

    denom = cfg.eta * (cfg.rho30 - cfg.rho10)

    egamma = (
        cfg.rho30 * mean_10
        - cfg.rho10 * mean_30
        - (cfg.rho30 * cfg.kappa10 - cfg.rho10 * cfg.kappa30)
        * 0.5
        * nubar_interp
        / cfg.alpha
    ) / denom

    # Stat: équivalent raisonnable du C++, mais en utilisant le nombre réel d’entrées
    # du bin au lieu du vecteur stat codé en dur.
    stat = std_30 / np.sqrt(np.maximum(counts_30, 1.0))

    # Dérivées utilisées dans ton code C++
    erf_eta = -egamma / cfg.eta

    erf_rho30 = (
        mean_10
        - cfg.kappa10 * nubar_interp / cfg.alpha
        - cfg.eta * egamma
    ) / denom

    erf_rho10 = (
        -mean_30
        + cfg.kappa30 * nubar_interp / cfg.alpha
        + cfg.eta * egamma
    ) / denom

    erf_kappa30 = -cfg.rho30 * nubar_interp / (
        cfg.alpha * cfg.eta * (cfg.rho30 - cfg.rho10)
    )

    erf_kappa10 = cfg.rho10 * nubar_interp / (
        cfg.alpha * cfg.eta * (cfg.rho30 - cfg.rho10)
    )

    erf_nubar = (
        cfg.rho30 * cfg.kappa10 - cfg.rho10 * cfg.kappa30
    ) / (
        cfg.alpha * cfg.eta * (cfg.rho30 - cfg.rho10)
    )

    sys = np.sqrt(
        (erf_eta * cfg.sigma_eta) ** 2
        + (erf_rho10 * cfg.sigma_rho10) ** 2
        + (erf_rho30 * cfg.sigma_rho30) ** 2
        + (erf_kappa10 * cfg.sigma_kappa10) ** 2
        + (erf_kappa30 * cfg.sigma_kappa30) ** 2
        + (erf_nubar * cfg.sigma_nubar) ** 2
    )

    return pd.DataFrame(
        {
            "energy": energy,
            "egamma_tot": egamma,
            "egamma_tot_stat": stat,
            "egamma_tot_sys": sys,
            "egamma_tot_min": egamma - sys,
            "egamma_tot_max": egamma + sys,
            "mean_egamma_10ns": mean_10,
            "mean_egamma_30ns": mean_30,
            "counts_10ns": counts_10,
            "counts_30ns": counts_30,
        }
    )