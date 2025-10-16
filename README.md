# **GASCONE : Gamma-rays Analysis for SCONE**  

## **Description**  

GASCONE is a high-level analysis module for SCONE [1]. It is intended to unfold and analyse $\gamma$-rays observables from neutron-induced fission experiments involving SCONE.

## **Inputs**  

GASCONE's main input is a CSV file from low-level analysis, containing the complete multplicity distribution of assemblies fired within a 50 ns coincidence-window opened by fission chamber triggers and binned by incident neutron energy [2]. The file is to be stored in ***data/scone/***.

GASCONE also uses external data.

- Existing measurements of fission $\gamma$-rays in the literature [3, 4] in ***data/scone/***.
- Average neutron-multiplicity from the evaluated library JEFF-4.0 [5] in ***data/evaluations/***.
- Fission codes simulations data from CGMF [6] and GEF [7] in ***data/simulations/***.
- Geant4 simulations of SCONE response to fission $\gamma$-rays cascades ***data/geant4/***.

## **Outputs**  

GASCONE's output is a CSV file of unfolded average $\gamma$-rays multiplicities saved in ***data/scone/***. Two plots will also be created in ***figs/***, including $\gamma$-rays multiplicities vs. incident energy and incident angular momentum [8].

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

- The fission $\gamma$-rays cascade files for fitting A and B (SCONE's response constants to $\gamma$-rays).
- The SCONE's response constant to neutrons (C) and its uncertainty (DC).
- The explored angular momentum range taken away by prompt neutrons (SN_MIN, SN_MAX).
- The font size for plots (FONT_SIZE).

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