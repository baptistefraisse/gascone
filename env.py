""" Global environment variables and data"""

# librairies 

from utils import *
from response import fit_scone_gconst_multiple

# plots font size

FONT_SIZE = 30

# physics constants 

SN = 0.5
SNU_MIN = 0.0
SNU_MAX = 0.5

# Geant4 files

FILENAMES = [
    "Geant4_CGMF_235U_10MeV.txt",
    "Geant4_GEF_238U_20MeV.txt",
    "Geant4_CGMF_238U_1MeV.txt",
    "Geant4_FIFRELIN_252Cf.txt"
    ]

# SCONE constants 

A, B, DA, DB = fit_scone_gconst_multiple(FILENAMES)
C, DC = 0.33, 0.01

# Microscopic calculations result for angular momentum transfer to fragments

FRAG_MOM_MICRO = 0.3

# evaluations : nubar JEFF-4.1

_, _, nubar_jeff, nubar_jeff_err = csv_data_reader(EVAL_DIR/"238U_nubar_JEFF4.csv")

# literature data : Qi

qi_energies, _, qi_mult, qi_mult_err = csv_data_reader(LIT_DIR/"Qi_200keV_data.csv")

# literature data: Laborie

laborie_energies, laborie_energies_err, laborie_mult, laborie_mult_err = csv_data_reader(LIT_DIR/"Laborie_190keV_data.csv")

# simulations: CGMF

cgmf_energies, _, cgmf_mult, _ = csv_data_reader(SIMU_DIR/"238U_CGMF_200keV.csv")

# simulations: GEF

gef_energies, _, gef_mult, _ = csv_data_reader(SIMU_DIR/"238U_GEF_200keV.csv")