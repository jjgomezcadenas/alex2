"""
rp_lvsh.py

Creates histograms comparing the straight-line and helix models.

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

base_name = "magse";
Bfield = "02";

fnb_trkL = "{0}/{1}{2}L/trk".format(dat_outdir,base_name,Bfield);
fnb_trkH = "{0}/{1}{2}H/trk".format(dat_outdir,base_name,Bfield);
plt_base = "{0}/{1}{2}L/plt".format(dat_outdir,base_name,Bfield);

dlim = 100;

# Ensure the track directory exists.
if(not os.path.isdir(fnb_trkL) or not os.path.isdir(fnb_trkH)):
    print "ERROR: track directory {0} or {1} not available".format(fnb_trkL,fnb_trkH);
    sys.exit();

# Ensure the plot directory exists.
if(not os.path.isdir(plt_base)):
    print "Creating plot directory {0}".format(plt_base);
    os.mkdir(plt_base);

# Plot options.
plt_show = False;
plt_print = True;

# Run the profile analysis for each track.
prof_konlist = [];
dlprof_xlist = [];
dhprof_xlist = [];
lhprof_xlist = [];
lhfprof_xlist = [];
for ntrk in range(num_tracks):

    if(not os.path.isfile("{0}/{1}{2}L/trk/{3}{4}L_f{5}.dat".format(dat_outdir,base_name,Bfield,base_name,Bfield,ntrk)) or not os.path.isfile("{0}/{1}{2}H/trk/{3}{4}H_f{5}.dat".format(dat_outdir,base_name,Bfield,base_name,Bfield,ntrk))):
      continue;
    
    print "-- Analysis for track {0}\n".format(ntrk);

    # Read in the straight-line model track.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    if(rev_trk):
        ltrktbl = np.loadtxt("{0}/{1}{2}L/trk/{3}{4}L_r{5}.dat".format(dat_outdir,base_name,Bfield,base_name,Bfield,ntrk));
    else:
        ltrktbl = np.loadtxt("{0}/{1}{2}L/trk/{3}{4}L_f{5}.dat".format(dat_outdir,base_name,Bfield,base_name,Bfield,ntrk));
    ltrk_node = ltrktbl[:,0];
    ltrk_xM = ltrktbl[:,1];
    ltrk_yM = ltrktbl[:,2];
    ltrk_zM = ltrktbl[:,3];
    ltrk_xP = ltrktbl[:,4];
    ltrk_yP = ltrktbl[:,5];
    ltrk_zP = ltrktbl[:,6];
    ltrk_chi2P = ltrktbl[:,7];
    ltrk_xF = ltrktbl[:,8];
    ltrk_yF = ltrktbl[:,9];
    ltrk_zF = ltrktbl[:,10];
    ltrk_chi2F = ltrktbl[:,11];

    # Read in the helix model track.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    if(rev_trk):
        htrktbl = np.loadtxt("{0}/{1}{2}H/trk/{3}{4}H_r{5}.dat".format(dat_outdir,base_name,Bfield,base_name,Bfield,ntrk));
    else:
        htrktbl = np.loadtxt("{0}/{1}{2}H/trk/{3}{4}H_f{5}.dat".format(dat_outdir,base_name,Bfield,base_name,Bfield,ntrk));
    htrk_node = htrktbl[:,0];
    htrk_xM = htrktbl[:,1];
    htrk_yM = htrktbl[:,2];
    htrk_zM = htrktbl[:,3];
    htrk_xP = htrktbl[:,4];
    htrk_yP = htrktbl[:,5];
    htrk_zP = htrktbl[:,6];
    htrk_chi2P = htrktbl[:,7];
    htrk_xF = htrktbl[:,8];
    htrk_yF = htrktbl[:,9];
    htrk_zF = htrktbl[:,10];
    htrk_chi2F = htrktbl[:,11];

    # Ensure the same number of measurements in each run.    
    nlvals = len(ltrk_node);
    nhvals = len(htrk_node);
    if(nlvals != nhvals):
        print "ERROR: number of measurements in L and H runs do not match for {0}".format(ntrk);
        sys.exit();

    # -----------------------------------------------------------------------
    # Record the k/N and relevant values.
    nn = 0;
    while(nn < nlvals):

        dl_dx = ltrk_xM[nn] - ltrk_xP[nn];
        dl_dy = ltrk_yM[nn] - ltrk_yP[nn];
        dl_dz = ltrk_zM[nn] - ltrk_zP[nn];
        dl = sqrt(dl_dx**2 + dl_dy**2 + dl_dz**2);

        dh_dx = htrk_xM[nn] - htrk_xP[nn];
        dh_dy = htrk_yM[nn] - htrk_yP[nn];
        dh_dz = htrk_zM[nn] - htrk_zP[nn];
        dh = sqrt(dh_dx**2 + dh_dy**2 + dh_dz**2);

        lh_dx = ltrk_xP[nn] - htrk_xP[nn];
        lh_dy = ltrk_yP[nn] - htrk_yP[nn];
        lh_dz = ltrk_zP[nn] - htrk_zP[nn];
        lh = sqrt(lh_dx**2 + lh_dy**2 + lh_dz**2);

        lhf_dx = ltrk_xF[nn] - htrk_xF[nn];
        lhf_dy = ltrk_yF[nn] - htrk_yF[nn];
        lhf_dz = ltrk_zF[nn] - htrk_zF[nn];
        lhf = sqrt(lhf_dx**2 + lhf_dy**2 + lhf_dz**2);
       
        if(dl < dlim and dh < dlim and lh < dlim):
            prof_konlist.append(1.0*nn/nlvals);
            dlprof_xlist.append(dl);
            dhprof_xlist.append(dh);
            lhprof_xlist.append(lh);
            lhfprof_xlist.append(lhf);

        nn += 1;

# ---------------------------------------------------------------------------
# Create the profiles vs. k/N.
prof_kon = []; prof_nvals = [];
dlprof_x = []; dlprof_sigma = [];
dhprof_x = []; dhprof_sigma = [];
lhprof_x = []; lhprof_sigma = [];
lhfprof_x = []; lhfprof_sigma = [];

for nn in range(nbins_kon):
    prof_kon.append(1.0*nn/nbins_kon);
    prof_nvals.append(0);

    dlprof_x.append(0.); dlprof_sigma.append(0.);
    dhprof_x.append(0.); dhprof_sigma.append(0.);
    lhprof_x.append(0.); lhprof_sigma.append(0.);
    lhfprof_x.append(0.); lhfprof_sigma.append(0.);

for kon,dl,dh,lh,lhf in zip(prof_konlist,dlprof_xlist,dhprof_xlist,lhprof_xlist,lhfprof_xlist):
    bb = int(kon*nbins_kon);
    prof_nvals[bb] += 1;
    dlprof_x[bb] += dl; dlprof_sigma[bb] += dl**2;
    dhprof_x[bb] += dh; dhprof_sigma[bb] += dh**2;
    lhprof_x[bb] += lh; lhprof_sigma[bb] += lh**2;
    lhfprof_x[bb] += lhf; lhfprof_sigma[bb] += lhf**2;

# Normalize.
for bb in range(nbins_kon):
    if(prof_nvals[bb] > 1):
        NN = prof_nvals[bb];
        mu_dl = dlprof_x[bb]/prof_nvals[bb];
        mu_dh = dhprof_x[bb]/prof_nvals[bb];
        mu_lh = lhprof_x[bb]/prof_nvals[bb];
        mu_lhf = lhfprof_x[bb]/prof_nvals[bb];
        dlprof_x[bb] = mu_dl;
        dhprof_x[bb] = mu_dh;
        lhprof_x[bb] = mu_lh;
        lhfprof_x[bb] = mu_lhf; 
        dlprof_sigma[bb] = sqrt((dlprof_sigma[bb])/(NN-1) - NN*mu_dl**2/(NN-1))/sqrt(NN);
        dhprof_sigma[bb] = sqrt((dhprof_sigma[bb])/(NN-1) - NN*mu_dh**2/(NN-1))/sqrt(NN);
        lhprof_sigma[bb] = sqrt((lhprof_sigma[bb])/(NN-1) - NN*mu_lh**2/(NN-1))/sqrt(NN);
        lhfprof_sigma[bb] = sqrt((lhfprof_sigma[bb])/(NN-1) - NN*mu_lhf**2/(NN-1))/sqrt(NN);
 
# Histogram the dl, dh, and lh values.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(dlprof_xlist, 50, normed=0, histtype='step',color='blue',label='|data - linear|');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("|pos(data) - pos(straight line)|");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_dl.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(dhprof_xlist, 50, normed=0, histtype='step',color='blue',label='|data - helix|');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("|pos(data) - pos(helix)|");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_dh.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(lhprof_xlist, 50, normed=0, histtype='step',color='blue',label='|linear - helix|');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("|pos(straight line) - pos(helix)|");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_lh.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(4);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(lhfprof_xlist, 50, normed=0, histtype='step',color='blue',label='|linear - helix|');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("filtered |pos(straight line) - pos(helix)|");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_lhf.pdf".format(plt_base), bbox_inches='tight');

# Draw the profiles.
fig = plt.figure(5);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, dlprof_x, '-', color='black', label='|data - linear|');
plt.title("|data - linear| profile");
lnd = plt.legend(loc=1,frameon=False,handletextpad=0); # loc=4 is bottom-right of plot
plt.xlabel("Fraction of track fit");
plt.ylabel("|pos(data) - pos(straight line)| average");
plt.savefig("{0}/prof_dl.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(6);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, dhprof_x, '-', color='black', label='|data - helix|');
plt.title("|data - helix| profile");
lnd = plt.legend(loc=1,frameon=False,handletextpad=0); # loc=4 is bottom-right of plot
plt.xlabel("Fraction of track fit");
plt.ylabel("|pos(data) - pos(helix)| average");
plt.savefig("{0}/prof_dh.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(7);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, lhprof_x, '-', color='black', label='|linear - helix|');
plt.title("|straight line - helix| profile");
lnd = plt.legend(loc=1,frameon=False,handletextpad=0); # loc=4 is bottom-right of plot
plt.xlabel("Fraction of track fit");
plt.ylabel("|pos(straight line) - pos(helix)| average");
plt.savefig("{0}/prof_lh.pdf".format(plt_base), bbox_inches='tight');

fig = plt.figure(8);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, lhfprof_x, '-', color='black', label='|linear - helix|');
plt.title("filtered |straight line - helix| profile");
lnd = plt.legend(loc=1,frameon=False,handletextpad=0); # loc=4 is bottom-right of plot
plt.xlabel("Fraction of track fit");
plt.ylabel("filtered |pos(straight line) - pos(helix)| average");
plt.savefig("{0}/prof_lhf.pdf".format(plt_base), bbox_inches='tight');
