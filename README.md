# **GASCONE : Gamma-rays Analysis for SCONE**  

## **Description**  

GASCONE is a high-level analysis module for SCONE [1]. It is intended to unfold and analyse $\gamma$-rays observables from fission experiments with SCONE.

## **Inputs**  

The main input is a CSV file produced by low-level analysis, containing the complete multplicity distribution of assemblies fired within a 50 ns coincidence-window after a fission chamber trigger and binned by incident neutron energy [2]. This input file is stored in ***data/scone/***.

External data are also used:

- In ***data/literature/***: existing measurements of fission $\gamma$-rays in the literature [3, 4].
- In ***data/evaluations/***: evaluated fission neutron multiplicities from JEFF-4.0 [5].
- In ***data/simulations/***: $\gamma$-rays from fission codes such as CGMF [6] and GEF [7].
- In ***data/geant4/***: Geant4 simulations of SCONE response to fission $\gamma$-rays.

## **Outputs**  

The main output is a CSV file of unfolded average $\gamma$-rays multiplicities saved in ***outputs/***. Some plots will also be generated: the fit of SCONE response to fission cascades, $\gamma$-rays multiplicities as a function of incident energy, $\gamma$-rays multiplicities as a function of the fissioning system angular momentum.

## **Installation and use**

Clone the repository, create a curated virtual environment and launch the main file.

```bash
git clone git@github.com:baptistefraisse/gascone.git
cd gascone
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 main.py
```

All the parameters the user would like to change are in ***env.py*** :

- The files of Geant4 simulations to SCONE response to fission cascade. 
- The SCONE's response constant to neutrons (C) and its uncertainty (DC).
- The range of angular momentum of prompt neutrons (SN_MIN, SN_MAX).
- The font size of plots (FONT_SIZE).

## **To do list**
 
- [short-term] Include G.Scamps simulation data and fit them directly in the code.
- [short-term] Include neutron cascades simulations and fit C directly in the code.
- [short-term] Include statistical errors.
- [long-term] Enable $\gamma$-rays analysis for spontaneous fission experiments.
- [long-term] Include the $\gamma$-rays total energy.

## **References**

[1] G. Belier _et al._, Nucl. Instrum. Methods Phys. Res. A **1072**, 170225 (2025)

[2] B.Fraisse, PhD thesis, Univ. Paris-Saclay (2024)

[3] L. Qi. PhD thesis, Univ. Paris-Saclay (2018).

[4] J.-M. Laborie et al. Phys. Rev. C **98**, 054604 (2018).

[5] Joint Evaluated Fission and Fusion Project. (2025). JEFF-4.0 Evaluated Data: neutron data [Data set]. OECD Nuclear Energy Agency. https://doi.org/10.82555/e9ajn-a3p20

[6] P. Talou et al. Comp. Phys. Comm. **269**, 108087 (2021).

[7] K.-H. Schmidt and B. Jurado et al. Nucl. Data Sheets **131**, 107–221 (2016).

[8] Article to come.

## **Contact**

fraisse@cua.edu