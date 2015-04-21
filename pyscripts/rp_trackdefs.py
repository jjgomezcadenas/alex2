"""
trackdefs.py

Function definitions and constants for tracking.

"""
import sys
import numpy as np
import scipy.integrate as integrate
import random as rd
import os
from math import *
import logging 

# -------------------------------------------------------------------------------------------------------
# Frequently modified parameters:
# -------------------------------------------------------------------------------------------------------

rev_trk = False;  # Set to true to fit reverse tracks and false to fit forward tracks (for kftrackfit.py only)

# Output directories
dat_outdir = "/data4/NEXT/users/jrenner/kalmanfilter/alex2/build/alexMain/out";    # output directory for tracks

run_name = "nmagbb2";   # name assigned to this run; will be used in naming the output files
num_tracks = 10;     # number of tracks to generate and/or fit

chi2_low = 0.0; #1.0e-5;     # lower chi2 boundary for certain chi2 analyses
chi2_outlier = 500.;     # upper chi2 boundary for chi2 profile and other analyses
cfxy_low = 1.0e-5;     # similar to chi2_low but for for cfxy
cfxy_outlier = 10000.; # similar to chi2_outlier but for cfxy

sigma_xm = 0.1;       # x measurement error for fitting to toyMC tracks
sigma_ym = 0.1;       # y measurement error for fitting to toyMC tracks

nfits = 1;

# -------------------------------------------------------------------------------------------------------
# Less frequently modified parameters:
# -------------------------------------------------------------------------------------------------------

prof_name = "magseHtest01_1mm_5sp";  # name of run used to generate the profiles (only relevant if running fitprof.py)
E_0 = 2.447+0.511;         # initial energy in MeV
MS0 = 13.6;        # multiple scattering parameter (should be 13.6)
eslice = 0.05;     # slice energy in MeV
E_tol = 1e-3;      # energy tolerance value (energy is considered to be 0 if less than this value)

chi2_lim = 40.;     # chi2 limit: a new segment is started if chi2 jumps above this limit during fit

# Parameters for the chi2 vs. k/N profiles.
nbins_kon = 30;
kon_min = 0.05;
kon_max = 0.95;

# Plot options.
plt_tracks = True;
plt_prediction = True;
plt_filtered = True;
plt_smearedHits = True;  # if False, will plot TrueHits
plt_chi2 = False;
plt_units = "mm";

plt_show = False;
plt_print = True;
plt_3dprint = True;

me = 0.511;       # electron rest mass in MeV

# Physics constants.
pc_rho0 = 2.6867774e19;   # density of ideal gas at T=0C, P=1 atm in cm^(-3)
pc_m_Xe = 131.293;        # mass of xenon in amu
pc_NA = 6.02214179e23;    # Avogadro constant

# Segment statistical analysis parameters.
stat_nseg = 4;    # number of segments used to construct the forward-reverse discriminating param
stat_efac = 0.5;  # fraction of total track that is considered to be in the first part of the track
