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

# Record the number of tracks that were actually processed.
nprocessed = 0;

# Extract data from the tracks files for num_tracks tracks.
splot_fchi2 = []; splot_rchi2 = [];

failfrac_f = []; failfrac_r = [];

prof_fkon = []; prof_fchi2 = [];
prof_rkon = []; prof_rchi2 = [];
prof_fekon = []; prof_fe = []; 
prof_rekon = []; prof_re = [];
all_fchi2 = []; all_rchi2 = [];
for ntrk in range(num_tracks):

    if(not os.path.isfile("{0}/{1}/trk/{2}_f{3}.dat".format(dat_outdir,run_name,run_name,ntrk))):
      continue;
    
    print "-- Chi2 analysis for track {0}\n".format(ntrk);
        
    # Read in the forward fit.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    fittbl = np.loadtxt("{0}/{1}/trk/{2}_f{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
    fit_chi2f = fittbl[:,11];
    fit_e = fittbl[:,12];
 
    # Read in the reverse fit.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    rfittbl = np.loadtxt("{0}/{1}/trk/{2}_r{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
    rfit_chi2f = rfittbl[:,11];
    rfit_e = rfittbl[:,12];

    # Record the k/N and chi2 and cfxy values for points in the forward track
    #  with valid chi2.  Also note the number of failed fits.
    nn = 0;
    fchi2avg = 0.; navg = 0;
    ntotpts = len(fit_chi2f);
    fnfail = 0;
    while(nn < ntotpts):
        chi2 = fit_chi2f[nn];
        if(chi2 > 0. and chi2 < 200.):
            all_fchi2.append(chi2);
        if(chi2 > 0. and chi2 < chi2_outlier):
            prof_fkon.append(1.0*nn/ntotpts);
            prof_fchi2.append(chi2);        # note: multiplication by 2 is a temporary fix
            fchi2avg += chi2; navg += 1
        if(chi2 < 0.): fnfail += 1;
        prof_fekon.append(1.0*nn/ntotpts);
        prof_fe.append(fit_e[nn]);
        nn += 1;
    if(fchi2avg/navg < 0.5): print "forward fchi2avg = {0}".format(fchi2avg/navg);
    #print "FWD: Found {0} failures / {1} points".format(fnfail,ntotpts);
    if(fnfail > 0): failfrac_f.append(1.0*fnfail/ntotpts); 

    # Record the k/N and chi2 average for points in the reverse track
    #  with valid chi2.
    rnn = 0;
    rchi2avg = 0.; rnavg = 0;
    rntotpts = len(rfit_chi2f);
    rnfail = 0;
    while(rnn < rntotpts):
        rchi2 = rfit_chi2f[rnn];
        if(rchi2 > 0. and rchi2 < 200.):
            all_rchi2.append(rchi2);
        if(rchi2 > 0. and rchi2 < chi2_outlier):
            prof_rkon.append(1.0*rnn/rntotpts);
            prof_rchi2.append(rchi2);         # note: multiplication by 2 is a temporary fix
            rchi2avg += rchi2; rnavg += 1
        if(rchi2 < 0.): rnfail += 1;
        prof_rekon.append(1.0*rnn/rntotpts);
        prof_re.append(rfit_e[rnn]);
        rnn += 1;
    if(rchi2avg/rnavg < 0.5): print "reverse chi2 average is {0}".format(np.mean(rchi2avg/rnavg))
    #print "REV: Found {0} failures / {1} points".format(rnfail,rntotpts);
    if(rnfail > 0): failfrac_r.append(1.0*rnfail/rntotpts);

    nprocessed += 1;

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
# Create the energy vs. k/N profiles.
prof_ekon = [];
fprof_envals = []; fprof_e = []; fprof_esigma = [];
rprof_envals = []; rprof_e = []; rprof_esigma = [];

for nn in range(nbins_kon):
    prof_ekon.append(1.0*nn/nbins_kon);

    fprof_envals.append(0); fprof_e.append(0.); fprof_esigma.append(0.);
    rprof_envals.append(0); rprof_e.append(0.); rprof_esigma.append(0.);

for kon,fe in zip(prof_fekon,prof_fe):
    bb = int(kon*nbins_kon);
    fprof_envals[bb] += 1;
    fprof_e[bb] += fe;
    fprof_esigma[bb] += fe**2;

for kon,re in zip(prof_rekon,prof_re):
    bb = int(kon*nbins_kon);
    rprof_envals[bb] += 1;
    rprof_e[bb] += re;
    rprof_esigma[bb] += re**2;

# Normalize.
for bb in range(nbins_kon):
    if(fprof_envals[bb] > 1):
        NN = fprof_envals[bb];
        mu = fprof_e[bb]/fprof_envals[bb];
        fprof_e[bb] = mu;
        fprof_esigma[bb] = sqrt((fprof_esigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);
    if(rprof_envals[bb] > 1):
        NN = rprof_envals[bb];
        mu = rprof_e[bb]/rprof_envals[bb];
        rprof_e[bb] = mu;
        rprof_esigma[bb] = sqrt((rprof_esigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);

# ---------------------------------------------------------------------------

# Plot and print the chi2 profiles.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, fprof_chi2, '-', color='black', label='Forward fit');
plt.plot(prof_kon, rprof_chi2, '--', color='black', label='Reverse fit');
#plt.title("chi2 profile: $\chi^2$ limit = {0}".format(chi2_outlier));
lnd = plt.legend(loc=1,frameon=False,handletextpad=0); # loc=4 is bottom-right of plot
plt.xlabel("Fraction of track fit");
plt.ylabel("$\chi^{2}$ local average");
plt.savefig("{0}/fit_profiles.pdf".format(plt_base), bbox_inches='tight');

# Plot the energy profiles.
fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_ekon, fprof_e, '-', color='black', label='Forward fit');
plt.plot(prof_ekon, rprof_e, '--', color='black', label='Reverse fit');
#plt.title("chi2 profile: $\chi^2$ limit = {0}".format(chi2_outlier));
lnd = plt.legend(loc=4,frameon=False,handletextpad=0);
plt.xlabel("Fraction of track fit");
plt.ylabel("q/p (MeV$^{-1}$)");
plt.savefig("{0}/energy_profiles.pdf".format(plt_base), bbox_inches='tight');

# ---------------------------------------------------------------------------
# Histogram all chi2 values

fig = plt.figure(3);
ax2 = fig.add_subplot(111);
chi2fn, chi2fbins, chi2fpatches = plt.hist(all_fchi2, 1000, normed=0, histtype='step',color='blue',label='Forward fit');
chi2rn, chi2rbins, chi2rpatches = plt.hist(all_rchi2, 1000, normed=0, histtype='step',color='red',label='Reverse fit');
ax2.set_xlabel("$\chi^{2}$");
ax2.set_ylabel("Counts/bin");
ax2.set_xlim(0,10);
#plt.yscale("log");
plt.savefig("{0}/chi2_histograms.pdf".format(plt_base), bbox_inches='tight');

# ---------------------------------------------------------------------------
# Histogram all failure fractions.

if(len(failfrac_f) == 0 or len(failfrac_r) == 0):
    print "No failure fraction histograms will be made, as {0} forward fits and {1} reverse fits had node failures.".format(len(failfrac_f),len(failfrac_r));
else:
    fig = plt.figure(4);
    ax2 = fig.add_subplot(111);
    chi2fn, chi2fbins, chi2fpatches = plt.hist(failfrac_f, 40, normed=0, histtype='step',color='blue',label='Forward fit');
    chi2rn, chi2rbins, chi2rpatches = plt.hist(failfrac_r, 40, normed=0, histtype='step',color='red',label='Reverse fit');
    lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
    ax2.set_xlabel("Fraction of failed fits");
    ax2.set_ylabel("Counts/bin");
    #ax2.set_xlim(0,1);
    plt.savefig("{0}/failed_fits_fraction.pdf".format(plt_base), bbox_inches='tight');
    print "Note: {0} forward fits and {1} reverse fits had node failures.".format(len(failfrac_f),len(failfrac_r));

print "{0} tracks in total were processed".format(nprocessed);
