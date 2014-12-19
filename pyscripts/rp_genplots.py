"""
rp_genplots.py

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

# Create the plt directory if it does not exist.
if(not os.path.isdir("{0}".format(plt_base))):
    os.mkdir("{0}".format(plt_base));
    print "Creating plt directory {0} ...".format(plt_base);

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
posfrac_f = []; posfrac_r = [];   # fraction of total points with sign > 0

lseg_f = []; lseg_r = [];
lseg1_f = []; lseg1_r = [];
lseg2_f = []; lseg2_r = [];
lseg3_f = []; lseg3_r = [];
lseg_qoverp_f = []; lseg_qoverp_r = [];

prof_fkon = []; prof_fchi2 = [];
prof_rkon = []; prof_rchi2 = [];
prof_fekon = []; prof_fe = []; 
prof_rekon = []; prof_re = [];
prof_fskon = []; prof_fsgn = [];
prof_rskon = []; prof_rsgn = [];
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
    fit_fqoverp = fittbl[:,13];
    fit_brk = fittbl[:,15];

    # Read in the reverse fit.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    rfittbl = np.loadtxt("{0}/{1}/trk/{2}_r{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
    rfit_chi2f = rfittbl[:,11];
    rfit_e = rfittbl[:,12];
    rfit_fqoverp = rfittbl[:,13];
    rfit_brk = rfittbl[:,15];

    # -----------------------------------------------------------------------------------------
    # Determine the length of the segments between breaks.
    # -----------------------------------------------------------------------------------------
    
    # Forward fit
    nn = 0; ntotpts = len(fit_brk);
    lmax1 = -1; lmax2 = -1; lmax3 = -1;
    lastbrk = 0;
    while(nn < ntotpts):
        if(fit_brk[nn] == 1):
            seglen = nn-lastbrk;
            lastbrk = nn;

            # Update the maximum values.
            if(seglen > lmax1):
              lmax3 = lmax2; lmax2 = lmax1; lmax1 = seglen;
            elif(seglen > lmax2):
              lmax3 = lmax2; lmax2 = seglen;
            elif(seglen > lmax3):
              lmax3 = seglen;

            # Add the segment length to the inclusive histogram.
            lseg_f.append(seglen);

        nn += 1;

    # Add the top 3 segment lengths to histograms.
    if(lmax1 > 0): lseg1_f.append(lmax1);
    if(lmax2 > 0): lseg2_f.append(lmax2);
    if(lmax3 > 0): lseg3_f.append(lmax3);

    # Reverse fit
    rnn = 0; rntotpts = len(rfit_brk);
    lmax1 = -1; lmax2 = -1; lmax3 = -1;
    lastbrk = 0;
    while(rnn < rntotpts):
        if(rfit_brk[rnn] == 1):
            seglen = rnn-lastbrk;
            lastbrk = rnn;

            # Update the maximum values.
            if(seglen > lmax1):
              lmax3 = lmax2; lmax2 = lmax1; lmax1 = seglen;
            elif(seglen > lmax2):
              lmax3 = lmax2; lmax2 = seglen;
            elif(seglen > lmax3):
              lmax3 = seglen;

            # Add the segment length to the inclusive histogram.
            lseg_r.append(seglen);

        rnn += 1;

    # Add the top 3 segment lengths to histograms.
    if(lmax1 > 0): lseg1_r.append(lmax1);
    if(lmax2 > 0): lseg2_r.append(lmax2);
    if(lmax3 > 0): lseg3_r.append(lmax3);

    # Record the k/N and chi2 and cfxy values for points in the forward track
    #  with valid chi2.  Also note the number of failed fits.
    nn = 0;
    fchi2avg = 0.; navg = 0;
    ntotpts = len(fit_chi2f);
    fnfail = 0; fnpos = 0;
    while(nn < ntotpts):
        chi2 = fit_chi2f[nn];
        if(chi2 > 0. and chi2 < 200.):
            all_fchi2.append(chi2);
        if(chi2 > 0. and chi2 < chi2_outlier):
            prof_fkon.append(1.0*nn/ntotpts);
            prof_fchi2.append(chi2);        # note: multiplication by 2 is a temporary fix
            fchi2avg += chi2; navg += 1;
        #if(chi2 < 0.): fnfail += 1;
        if(fit_brk[nn] == 1): fnfail += 1;
        if(chi2 >= 0.):
          prof_fekon.append(1.0*nn/ntotpts);
          prof_fe.append(fit_e[nn]);
          prof_fskon.append(1.0*nn/ntotpts);
          #prof_fsgn.append(fit_fqoverp[nn]);
          if(fit_fqoverp[nn] > 0):
              prof_fsgn.append(1);
              fnpos += 1;
          else:
              prof_fsgn.append(-1);
        nn += 1;
    if(fchi2avg/navg < 0.5): print "forward fchi2avg = {0}".format(fchi2avg/navg);
    #print "FWD: Found {0} failures / {1} points".format(fnfail,ntotpts);
    if(fnfail > 0): failfrac_f.append(1.0*fnfail/ntotpts);
    posfrac_f.append(1.0*fnpos/ntotpts);

    # Record the k/N and chi2 average for points in the reverse track
    #  with valid chi2.
    rnn = 0;
    rchi2avg = 0.; rnavg = 0;
    rntotpts = len(rfit_chi2f);
    rnfail = 0; rnpos = 0;
    while(rnn < rntotpts):
        rchi2 = rfit_chi2f[rnn];
        if(rchi2 > 0. and rchi2 < 200.):
            all_rchi2.append(rchi2);
        if(rchi2 > 0. and rchi2 < chi2_outlier):
            prof_rkon.append(1.0*rnn/rntotpts);
            prof_rchi2.append(rchi2);         # note: multiplication by 2 is a temporary fix
            rchi2avg += rchi2; rnavg += 1
        #if(rchi2 < 0.): rnfail += 1;
        if(rfit_brk[rnn] == 1): rnfail += 1;
        if(rchi2 > 0.):
          prof_rekon.append(1.0*rnn/rntotpts);
          prof_re.append(rfit_e[rnn]);
          prof_rskon.append(1.0*rnn/rntotpts);
          #prof_rsgn.append(rfit_fqoverp[rnn]);
          if(rfit_fqoverp[rnn] > 0):
              prof_rsgn.append(1);
              rnpos += 1;
          else:
              prof_rsgn.append(-1);
        rnn += 1;
    if(rchi2avg/rnavg < 0.5): print "reverse chi2 average is {0}".format(np.mean(rchi2avg/rnavg))
    #print "REV: Found {0} failures / {1} points".format(rnfail,rntotpts);
    if(rnfail > 0): failfrac_r.append(1.0*rnfail/rntotpts);
    posfrac_r.append(1.0*rnpos/rntotpts);

    # Read in the largest segment fits if they exist.
    if(os.path.isfile("{0}/{1}/trk/{2}_flsegf{3}.dat".format(dat_outdir,run_name,run_name,ntrk))):

        # Read in the forward fit for the longest segment.
        lf_fittbl = np.loadtxt("{0}/{1}/trk/{2}_flsegf{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
        lf_fqoverp = lf_fittbl[:,13];

        lseg_qoverp_f.append(np.mean(lf_fqoverp));

        # Read in the reverse fit for the longest segment.
        lr_fittbl = np.loadtxt("{0}/{1}/trk/{2}_flsegr{3}.dat".format(dat_outdir,run_name,run_name,ntrk));
        lr_fqoverp = lr_fittbl[:,13];

        lseg_qoverp_r.append(np.mean(lr_fqoverp));
    else:
        print "--- No largest segment file found for track {0}".format(ntrk);

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
# Create the sign vs. k/N profiles.
prof_skon = [];
fprof_sgnvals = []; fprof_sgn = []; fprof_ssigma = [];
rprof_sgnvals = []; rprof_sgn = []; rprof_ssigma = [];

for nn in range(nbins_kon):
    prof_skon.append(1.0*nn/nbins_kon);

    fprof_sgnvals.append(0); fprof_sgn.append(0.); fprof_ssigma.append(0.);
    rprof_sgnvals.append(0); rprof_sgn.append(0.); rprof_ssigma.append(0.);

for kon,fsgn in zip(prof_fskon,prof_fsgn):
    bb = int(kon*nbins_kon);
    fprof_sgnvals[bb] += 1;
    fprof_sgn[bb] += fsgn;
    fprof_ssigma[bb] += fsgn**2;

for kon,rsgn in zip(prof_rskon,prof_rsgn):
    bb = int(kon*nbins_kon);
    rprof_sgnvals[bb] += 1;
    rprof_sgn[bb] += rsgn;
    rprof_ssigma[bb] += rsgn**2;

# Normalize.
for bb in range(nbins_kon):
    if(fprof_sgnvals[bb] > 1):
        NN = fprof_sgnvals[bb];
        mu = fprof_sgn[bb]/fprof_sgnvals[bb];
        fprof_sgn[bb] = mu;
        fprof_ssigma[bb] = sqrt((fprof_ssigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);
    if(rprof_sgnvals[bb] > 1):
        NN = rprof_sgnvals[bb];
        mu = rprof_sgn[bb]/rprof_sgnvals[bb];
        rprof_sgn[bb] = mu;
        rprof_ssigma[bb] = sqrt((rprof_ssigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);

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

# Plot the sign profiles.
fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_skon, fprof_sgn, '-', color='black', label='Forward fit');
plt.plot(prof_skon, rprof_sgn, '--', color='black', label='Reverse fit');
#plt.title("chi2 profile: $\chi^2$ limit = {0}".format(chi2_outlier));
lnd = plt.legend(loc=4,frameon=False,handletextpad=0);
plt.xlabel("Fraction of track fit");
plt.ylabel("Average sign");
plt.savefig("{0}/sign_profiles.pdf".format(plt_base), bbox_inches='tight');

# ---------------------------------------------------------------------------
# Histogram all chi2 values

fig = plt.figure(4);
ax2 = fig.add_subplot(111);
chi2fn, chi2fbins, chi2fpatches = plt.hist(all_fchi2, 1000, normed=0, histtype='step',color='blue',label='Forward fit');
chi2rn, chi2rbins, chi2rpatches = plt.hist(all_rchi2, 1000, normed=0, histtype='step',color='red',label='Reverse fit');
ax2.set_xlabel("$\chi^{2}$");
ax2.set_ylabel("Counts/bin");
ax2.set_xlim(0,10);
#plt.yscale("log");
plt.savefig("{0}/chi2_histograms.pdf".format(plt_base), bbox_inches='tight');

# ---------------------------------------------------------------------------
# Draw the segment length histograms.
if(len(lseg_f) > 0 and len(lseg_r) > 0):
  fig = plt.figure(5);
  ax4 = fig.add_subplot(111);
  seglenfn, seglenfbins, seglenfpatches = plt.hist(lseg_f, 300, normed=0, histtype='step',color='blue',label='Forward fit');
  seglenrn, seglenrbins, seglenrpatches = plt.hist(lseg_r, 300, normed=0, histtype='step',color='red',label='Reverse fit');
  lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
  ax4.set_xlabel("Segment length (inclusive)");
  ax4.set_ylabel("Counts/bin");
  ax4.set_xlim(0,300);
  #plt.yscale("log");
  plt.savefig("{0}/lseg_hist.pdf".format(plt_base), bbox_inches='tight');

if(len(lseg1_f) > 0 and len(lseg1_r) > 0):
  fig = plt.figure(6);
  ax5 = fig.add_subplot(111);
  seglenfn, seglenfbins, seglenfpatches = plt.hist(lseg1_f, 300, normed=0, histtype='step',color='blue',label='Forward fit');
  seglenrn, seglenrbins, seglenrpatches = plt.hist(lseg1_r, 300, normed=0, histtype='step',color='red',label='Reverse fit');
  lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
  ax5.set_xlabel("Segment length (longest segment)");
  ax5.set_ylabel("Counts/bin");
  ax5.set_xlim(0,300);
  #plt.yscale("log");
  plt.savefig("{0}/lseg1_hist.pdf".format(plt_base), bbox_inches='tight');

if(len(lseg2_f) > 0 and len(lseg2_r) > 0):
  fig = plt.figure(7);
  ax6 = fig.add_subplot(111);
  seglenfn, seglenfbins, seglenfpatches = plt.hist(lseg2_f, 300, normed=0, histtype='step',color='blue',label='Forward fit');
  seglenrn, seglenrbins, seglenrpatches = plt.hist(lseg2_r, 300, normed=0, histtype='step',color='red',label='Reverse fit');
  lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
  ax6.set_xlabel("Segment length (2nd longest segment)");
  ax6.set_ylabel("Counts/bin");
  ax6.set_xlim(0,300);
  #plt.yscale("log");
  plt.savefig("{0}/lseg2_hist.pdf".format(plt_base), bbox_inches='tight'); 

if(len(lseg3_f) > 0 and len(lseg3_r) > 0):
  fig = plt.figure(8);
  ax7 = fig.add_subplot(111);
  seglenfn, seglenfbins, seglenfpatches = plt.hist(lseg3_f, 300, normed=0, histtype='step',color='blue',label='Forward fit');
  seglenrn, seglenrbins, seglenrpatches = plt.hist(lseg3_r, 300, normed=0, histtype='step',color='red',label='Reverse fit');
  lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
  ax7.set_xlabel("Segment length (3rd longest segment)");
  ax7.set_ylabel("Counts/bin");
  ax7.set_xlim(0,300);
  #plt.yscale("log");
  plt.savefig("{0}/lseg3_hist.pdf".format(plt_base), bbox_inches='tight');

# ---------------------------------------------------------------------------
# Histogram all failure fractions.

if(len(failfrac_f) == 0 or len(failfrac_r) == 0):
    print "No failure fraction histograms will be made, as {0} forward fits and {1} reverse fits had node failures.".format(len(failfrac_f),len(failfrac_r));
else:
    fig = plt.figure(9);
    ax2 = fig.add_subplot(111);
    chi2fn, chi2fbins, chi2fpatches = plt.hist(failfrac_f, 20, normed=0, histtype='step',color='blue',label='Forward fit');
    chi2rn, chi2rbins, chi2rpatches = plt.hist(failfrac_r, 20, normed=0, histtype='step',color='red',label='Reverse fit');
    lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
    ax2.set_xlabel("Breakpoint fraction (breakpt nodes/total nodes)");
    ax2.set_ylabel("Counts/bin");
    #ax2.set_xlim(0,1);
    plt.savefig("{0}/failed_fits_fraction.pdf".format(plt_base), bbox_inches='tight');
    print "Note: {0} forward fits and {1} reverse fits had node failures.".format(len(failfrac_f),len(failfrac_r));

# ---------------------------------------------------------------------------
# Histogram all positive fractions.

if(len(posfrac_f) == 0 or len(posfrac_r) == 0):
    print "No positive fraction histograms will be made, as {0} forward fits and {1}.".format(len(posfrac_f),len(posfrac_r));
else:
    fig = plt.figure(10);
    ax2 = fig.add_subplot(111);
    posfracfn, posfracbins, posfracpatches = plt.hist(posfrac_f, 20, normed=0, histtype='step',color='blue',label='Forward fit');
    posfracrn, posfracrbins, posfracrpatches = plt.hist(posfrac_r, 20, normed=0, histtype='step',color='red',label='Reverse fit');
    lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
    ax2.set_xlabel("Positive fraction (nodes with positive q/p / total nodes)");
    ax2.set_ylabel("Counts/bin");
    #ax2.set_xlim(0,1);
    plt.savefig("{0}/positive_fraction.pdf".format(plt_base), bbox_inches='tight');
    #print "Note: {0} forward fits and {1} reverse fits had node failures.".format(len(posfrac_f),len(posfrac_r));

# ---------------------------------------------------------------------------
# Histogram all average q/p values for the longest segments.

fig = plt.figure(11);
ax2 = fig.add_subplot(111);
lsegqopfn, lsegqopfbins, lsegqopfpatches = plt.hist(lseg_qoverp_f, 30, normed=0, histtype='step',color='blue',label='Forward fit');
lsegqoprn, lsegqoprbins, lsegqoprpatches = plt.hist(lseg_qoverp_r, 30, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
ax2.set_xlabel("Average q/p");
ax2.set_ylabel("Counts/bin");
#ax2.set_xlim(0,1);
plt.savefig("{0}/lseg_avg_qoverp.pdf".format(plt_base), bbox_inches='tight');
print "Note (avg. q/p histogram): {0} entries in forward histogram and {1} entries in reverse histogram.".format(len(lseg_qoverp_f),len(lseg_qoverp_r));

print "{0} tracks in total were processed".format(nprocessed);
