"""
rp_fitprof.py

Uses the previously determined average chi2 profiles to fit on a 
track-by-track basis.

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from math import *
from rp_trackdefs import *
from scipy.interpolate import interp1d

from abc import ABCMeta, abstractmethod
import logging 

fnb_trk = "{0}/{1}/trk".format(dat_outdir,run_name);
plt_base = "{0}/{1}/plt".format(dat_outdir,run_name);

# Ensure the track directory exists.
if(not os.path.isdir(fnb_trk)):
    print "ERROR: track directory {0} not available".format(dat_outdir);
    sys.exit();

# Ensure the profiles directory exists.
if(not os.path.isdir("{0}/prof".format(plt_base))):
    print "ERROR: profiles directory not available";
    sys.exit();

# Plot options.
plt_show = False;
plt_print = True;

# Read in and interpolate the profiles.
fproftbl = np.loadtxt("{0}/{1}/plt/prof/prof_{2}_forward.dat".format(dat_outdir,run_name,prof_name));
fprof = interp1d(fproftbl[:,0],fproftbl[:,1],kind='linear');
fprof_sigma = interp1d(fproftbl[:,0],fproftbl[:,2],kind='linear');
rproftbl = np.loadtxt("{0}/{1}/plt/prof/prof_{2}_reverse.dat".format(dat_outdir,run_name,prof_name));
rprof = interp1d(rproftbl[:,0],rproftbl[:,1],kind='linear');
rprof_sigma = interp1d(rproftbl[:,0],rproftbl[:,2],kind='linear');

#kon_vals = [];
#fprof_vals = [];
#rprof_vals = [];
#for ii in np.arange(0.,0.9,0.01):
#    kon_vals.append(ii);
#    fprof_vals.append(fprof(ii))
#    rprof_vals.append(rprof(ii))
#    
#fig = plt.figure(2);
#fig.set_figheight(5.0);
#fig.set_figwidth(7.5);
#plt.plot(kon_vals,fprof_vals,'-',color='blue',markersize=1.0,label='Forward fit');
#plt.plot(kon_vals,rprof_vals,'-',color='blue',markersize=1.0,label='Reverse fit');
#lnd = plt.legend(loc=4,frameon=False,handletextpad=0);
#plt.title("Profile chi2 comparisons ($\chi^2$ limit = {0}, k/N from 0.1 to 0.9)".format(chi2_outlier));
#plt.xlabel("k/N");
#plt.ylabel("$\chi^{2}$");
#plt.show();

# Run the profile analysis for each track.
splot_fchi2F = []; splot_fchi2R = [];
splot_rchi2F = []; splot_rchi2R = [];
for ntrk in range(num_tracks):

    # Skip tracks that were not written.
    if(not os.path.isfile("{0}/{1}_f{2}.dat".format(fnb_trk,run_name,ntrk))):
      continue;
    
    print "-- Profile analysis for track {0}\n".format(ntrk);

    # Read in the forward track.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    trktbl = np.loadtxt("{0}/{1}_f{2}.dat".format(fnb_trk,run_name,ntrk));
    trk_node = trktbl[:,0];
    trk_xM = trktbl[:,1];
    trk_yM = trktbl[:,2];
    trk_zM = trktbl[:,3];
    trk_xP = trktbl[:,4];
    trk_yP = trktbl[:,5];
    trk_zP = trktbl[:,6];
    trk_chi2P = trktbl[:,7];
    trk_xF = trktbl[:,8];
    trk_yF = trktbl[:,9];
    trk_zF = trktbl[:,10];
    trk_chi2F = trktbl[:,11];

    # Read in the reverse track.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    rtrktbl = np.loadtxt("{0}/{1}_r{2}.dat".format(fnb_trk,run_name,ntrk));
    rtrk_node = rtrktbl[:,0];
    rtrk_xM = rtrktbl[:,1];
    rtrk_yM = rtrktbl[:,2];
    rtrk_zM = rtrktbl[:,3];
    rtrk_xP = rtrktbl[:,4];
    rtrk_yP = rtrktbl[:,5];
    rtrk_zP = rtrktbl[:,6];
    rtrk_chi2P = rtrktbl[:,7];
    rtrk_xF = rtrktbl[:,8];
    rtrk_yF = rtrktbl[:,9];
    rtrk_zF = rtrktbl[:,10];
    rtrk_chi2F = rtrktbl[:,11];
    
    # -----------------------------------------------------------------------
    # Forward fit:
        
    # Compute the fchi2F and fchi2R.
    nn = 0; ntotpts = len(trk_node);
    chi2F = 0.; chi2R = 0.; ndof = 0;
    for nn in range(ntotpts):
        kon = 1.0*nn/ntotpts;
        chi2 = trk_chi2F[nn];
        
        if(chi2 >= chi2_low and chi2 < chi2_outlier and kon > kon_min and kon < kon_max):
            chi2F += (chi2-fprof(kon))**2/fprof_sigma(kon)**2;
            chi2R += (chi2-rprof(kon))**2/rprof_sigma(kon)**2;
            ndof += 1;
    splot_fchi2F.append(chi2F/ndof);
    splot_fchi2R.append(chi2R/ndof);

    # -----------------------------------------------------------------------
    # Reverse fit:
        
    # Compute the rchi2F and rchi2R.
    nn = 0; ntotpts = len(rtrk_node);
    chi2F = 0.; chi2R = 0.; ndof = 0;
    for nn in range(ntotpts):
        kon = 1.0*nn/ntotpts;
        chi2 = rtrk_chi2F[nn];
        if(chi2 >= chi2_low and chi2 < chi2_outlier and kon > kon_min and kon < kon_max):
            chi2F += (chi2-fprof(kon))**2/fprof_sigma(kon)**2;
            chi2R += (chi2-rprof(kon))**2/rprof_sigma(kon)**2;
            ndof += 1;
    splot_rchi2F.append(chi2F/ndof);
    splot_rchi2R.append(chi2R/ndof);
    
# ---------------------------------------------------------------------------
# Compute the chi2 ratios
splot_fratio = []; splot_rratio = []; 
splot_fdiff = []; splot_rdiff = []; splot_ratiodiff = [];
for fchi2F,fchi2R,rchi2F,rchi2R in zip(splot_fchi2F,splot_fchi2R,splot_rchi2F,splot_rchi2R): 
    splot_fratio.append(fchi2F/fchi2R);
    splot_rratio.append(rchi2F/rchi2R);
    splot_fdiff.append(fchi2F-rchi2F);
    splot_rdiff.append(fchi2R-rchi2R);
    splot_ratiodiff.append(fchi2F/fchi2R - rchi2F/rchi2R);

# ---------------------------------------------------------------------------
# Determine the efficiency vs. 1-b for several cuts.
eff_vals = []; bgr_vals = [];
cut_vals = [];
c1max = max(splot_fratio);
c2max = max(splot_rratio);
cmax = max(c1max,c2max);
nrvals = len(splot_fratio);
for cval in np.arange(0,cmax,cmax/50.):
    
    # Determine the number of signal and background events passing the cuts.
    nsig = 0; nbg = 0;
    for fr,rr in zip(splot_fratio,splot_rratio):
        if(fr < cval): nsig += 1;
        if(rr < cval): nbg += 1;
    
    # Add the fractions of signal and background events to the list.
    eff_vals.append(1.0*nsig/nrvals);
    bgr_vals.append(1.0-1.0*nbg/nrvals);
        
    # Add the cut to the cut values list.
    cut_vals.append(cval);

#print "\n\n";
#print "# eff (1-b)\n";
#for eff,bgr in zip(eff_vals,bgr_vals):
#    print "{0} {1}".format(eff,bgr);
#print "\n\n";
    
# Plot the quantities in a scatter plot.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_fchi2F,splot_fchi2R,'.',color='blue',markersize=1.0,label='Forward fit');
plt.plot(splot_rchi2F,splot_rchi2R,'.',color='red',markersize=1.0,label='Reverse fit');
lnd = plt.legend(loc=2,frameon=False,handletextpad=0);
plt.title("Profile chi2 comparisons ($\chi^2$ limit = {0})".format(chi2_outlier));
plt.xlabel("$\chi^{2}_F$");
plt.ylabel("$\chi^{2}_R$");
plt.axis([0, 75000., 0, 75000.]);
plt.savefig("{0}/prof_compare.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_fchi2F,splot_fchi2R,'.',color='blue',markersize=1.0,label='Forward fit');
lnd = plt.legend(loc=2,frameon=False,handletextpad=0);
plt.title("Forward fit: Profile chi2 comparisons ($\chi^2$ limit = {0})".format(chi2_outlier));
plt.xlabel("$\chi^{2}_F$");
plt.ylabel("$\chi^{2}_R$");
plt.axis([0, 75000., 0, 75000.]);
plt.savefig("{0}/prof_compare_forward.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_rchi2F,splot_rchi2R,'.',color='red',markersize=1.0,label='Reverse fit');
lnd = plt.legend(loc=2,frameon=False,handletextpad=0);
plt.title("Reverse fit: Profile chi2 comparisons ($\chi^2$ limit = {0})".format(chi2_outlier));
plt.xlabel("$\chi^{2}_F$");
plt.ylabel("$\chi^{2}_R$");
plt.axis([0, 75000., 0, 75000.]);
plt.savefig("{0}/prof_compare_reverse.pdf".format(plt_base), bbox_inches='tight');

# Histogram the likelihood ratios.
fig = plt.figure(4);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(splot_fratio, 50, normed=0, histtype='step',color='blue',label='Forward fit');
krn, krbins, krpatches = plt.hist(splot_rratio, 50, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("ratio $\chi^{2}_{F}/\chi^{2}_{R}$");
plt.ylabel("Counts/bin");
plt.savefig("{0}/likelihood_ratios.pdf".format(plt_base), bbox_inches='tight');

# Histogram the differences.
fig = plt.figure(5);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(splot_ratiodiff, 50, normed=0, histtype='step',color='blue',label='Difference');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("difference in $\chi^{2}_{F}/\chi^{2}_{R}$");
plt.ylabel("Counts/bin");
plt.savefig("{0}/lratio_diff.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(6);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(splot_fdiff, 50, normed=0, histtype='step',color='blue',label='chi2F difference');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("difference in $\chi^{2}_{F}/\chi^{2}_{R}$");
plt.ylabel("Counts/bin");
plt.savefig("{0}/chi2F_diff.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(7);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(bgr_vals,eff_vals,'-',color='blue',markersize=1.0);
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("Reverse fit rejection (1-b)");
plt.ylabel("Forward fit identification efficiency ($\epsilon$)");
plt.savefig("{0}/signal_eff_vs_bgr.pdf".format(plt_base), bbox_inches='tight');
