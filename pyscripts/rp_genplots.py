"""
genplots.py

Plots the relevant statistical quantities for the forward and reverse tracks.

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from math import *
from rp_trackdefs import *

from abc import ABCMeta, abstractmethod
import logging 

# Ensure the track directory has been created.
if(not os.path.isdir(dat_outdir)): 
    print "ERROR: track directory {0} not available".format(dat_outdir);
    sys.exit();

plt_base = "{0}/{1}/plt".format(dat_outdir,run_name);

# Create the profiles directory if it does not exist.
if(not os.path.isdir("{0}/prof".format(plt_base))):
    os.mkdir("{0}/prof".format(plt_base)); 
    print "Creating profiles directory {0}/prof ...".format(plt_base);

# Plot options.
plt_show = False;
plt_print = True;

# Extract data from the tracks files for num_tracks tracks.
splot_fchi2 = []; splot_rchi2 = [];

prof_fkon = []; prof_fchi2 = [];
prof_rkon = []; prof_rchi2 = [];
for ntrk in range(num_tracks):
    
    print "-- Chi2 analysis for track {0}\n".format(ntrk);
        
    # Read in the forward fit.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    fittbl = np.loadtxt("{0}/{1}/trk/{2}_f{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
    fit_chi2f = fittbl[:,11];
 
    # Read in the reverse fit.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    rfittbl = np.loadtxt("{0}/{1}/trk/{2}_r{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
    rfit_chi2f = rfittbl[:,11];

    # Record the k/N and chi2 and cfxy values for points in the forward track
    #  with valid chi2.
    nn = 0;
    fchi2avg = 0.; navg = 0;
    ntotpts = len(fit_chi2f);
    while(nn < ntotpts):
        chi2 = fit_chi2f[nn];
        if(chi2 > 0. and chi2 < chi2_outlier):
            prof_fkon.append(1.0*nn/ntotpts);
            prof_fchi2.append(chi2);        # note: multiplication by 2 is a temporary fix
            fchi2avg += chi2; navg += 1
        nn += 1;
    if(fchi2avg/navg < 0.5): print "forward fchi2avg = {0}".format(fchi2avg/navg);
 
    # Record the k/N and chi2 average for points in the reverse track
    #  with valid chi2.
    rnn = 0;
    rchi2avg = 0.; rnavg = 0;
    rntotpts = len(rfit_chi2f);
    while(rnn < rntotpts):
        rchi2 = rfit_chi2f[rnn];
        if(rchi2 > 0. and rchi2 < chi2_outlier):
            prof_rkon.append(1.0*rnn/rntotpts);
            prof_rchi2.append(rchi2);         # note: multiplication by 2 is a temporary fix
            rchi2avg += rchi2; rnavg += 1
        rnn += 1;
    if(rchi2avg/rnavg < 0.5): print "reverse chi2 average is {0}".format(np.mean(rchi2avg/rnavg))

# ---------------------------------------------------------------------------
# Create the chi2 vs. k/N profiles.
prof_kon = [];
fprof_nvals = []; fprof_chi2 = []; fprof_sigma = [];
rprof_nvals = []; rprof_chi2 = []; rprof_sigma = [];

for nn in range(nbins_kon):
    prof_kon.append(1.0*nn/nbins_kon);

    fprof_nvals.append(0); fprof_chi2.append(0.); fprof_sigma.append(0.);
    rprof_nvals.append(0); rprof_chi2.append(0.); rprof_sigma.append(0.);

for kon,fchi2 in zip(prof_fkon,prof_fchi2):
    bb = int(kon*nbins_kon);
    fprof_nvals[bb] += 1;
    fprof_chi2[bb] += fchi2;
    fprof_sigma[bb] += fchi2**2;

for kon,rchi2 in zip(prof_rkon,prof_rchi2):
    bb = int(kon*nbins_kon);
    rprof_nvals[bb] += 1;
    rprof_chi2[bb] += rchi2;
    rprof_sigma[bb] += rchi2**2;

# Normalize.
for bb in range(nbins_kon):
    if(fprof_nvals[bb] > 1):
        NN = fprof_nvals[bb];
        mu = fprof_chi2[bb]/fprof_nvals[bb];
        fprof_chi2[bb] = mu;
        fprof_sigma[bb] = sqrt((fprof_sigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);
    if(rprof_nvals[bb] > 1):
        NN = rprof_nvals[bb];
        mu = rprof_chi2[bb]/rprof_nvals[bb];
        rprof_chi2[bb] = mu;
        rprof_sigma[bb] = sqrt((rprof_sigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);

# Print out the profiles.
f_fprof = open("{0}/prof/prof_{1}_forward.dat".format(plt_base,run_name),"w")
f_fprof.write("# (k/N) (chi2) (sigma) (N)\n")
for kon,chi2,sigma,NN in zip(prof_kon,fprof_chi2,fprof_sigma,fprof_nvals):
    f_fprof.write("{0} {1} {2} {3}\n".format(kon,chi2,sigma,NN));
f_fprof.close()

f_rprof = open("{0}/prof/prof_{1}_reverse.dat".format(plt_base,run_name),"w")
f_rprof.write("# (k/N) (chi2) (sigma) (N)\n")
for kon,chi2,sigma,NN in zip(prof_kon,rprof_chi2,rprof_sigma,rprof_nvals):
    f_rprof.write("{0} {1} {2} {3}\n".format(kon,chi2,sigma,NN));
f_rprof.close()

# ---------------------------------------------------------------------------

# Plot and print the chi2 profiles.
fig = plt.figure(11);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, fprof_chi2, '-', color='black', label='Forward fit');
plt.plot(prof_kon, rprof_chi2, '--', color='black', label='Reverse fit');
#plt.title("chi2 profile: $\chi^2$ limit = {0}".format(chi2_outlier));
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("Fraction of track fit");
plt.ylabel("$\chi^{2}$ local average");
plt.savefig("{0}/fit_profiles.pdf".format(plt_base), bbox_inches='tight');
