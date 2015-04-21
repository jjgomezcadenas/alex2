"""
rp_sevsbb.py

Single-electron vs. 0vbb events

Note: here the convention taken is: "passing" the cut means the event
is single-electron-like

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

# Get the track names as arguments.
args = sys.argv;
if(args[1] != "" and args[2] != ""):
    trk_name1 = args[1];
    trk_name2 = args[2];
else:
    trk_name1 = "magseHtest01_1mm_5sp"; # "magse05H_le_5sp_lseg_pfit";  # the single-electron event run
    trk_name2 = "magbbHtest01_1mm_5sp"; # "magbb05H_le_5sp_lseg_pfit";  # the double-beta event run

# Plot options.
plt_show = False;
plt_print = True;

#trk_outdir = "/Users/jrenner/IFIC/pylosk/tracks";
#trk_outdir = "/Users/jrenner/IFIC/software/alex2/build/alexMain/out";
trk_outdir = "/data4/NEXT/users/jrenner/kalmanfilter/alex2/build/alexMain/out";

plt_base1 = "{0}/{1}/plt_curv".format(trk_outdir,trk_name1);
plt_base2 = "{0}/{1}/plt_curv".format(trk_outdir,trk_name2);
apply_lpf = True;

fn_plt_sigvsb = "{0}/sigvsb_all.pdf".format(plt_base1);
fn_plt_curv = "{0}/scurv_diff_means.pdf".format(plt_base1);

# Number of points in combined curve.
ncpts = 50;

# Read in the mean signed curvature values.
scurvtbl1 = np.loadtxt("{0}/scurv_means.dat".format(plt_base1));
scurvtbl2 = np.loadtxt("{0}/scurv_means.dat".format(plt_base2));

l_scurv1 = scurvtbl1[:,1];
l_scurv2 = scurvtbl2[:,1];
    
# Cut ranges.
minCurv = min(min(l_scurv1), min(l_scurv2));
maxCurv = max(max(l_scurv1), max(l_scurv1));
curv_cuts = np.arange(minCurv, maxCurv, (maxCurv-minCurv)/100); # cut on curvature asymmetry values

# Determine which events passed the curvature cut.
# 0: would not pass the curvature cut
# 1: passes the curvature cut with currently selected orientation
# 2: passes the curvature cut if the sign of the asymmetry factor is flipped
t_cut_curv1 = []; t_cut_curv2 = [];
for mdiff_c in curv_cuts:
    
    cut_curv1 = []; cut_curv2 = [];
    
    for mdiff in l_scurv1:
        if(mdiff > mdiff_c): cut_curv1.append(1);
        else: cut_curv1.append(0);
    t_cut_curv1.append(cut_curv1);
    
    for mdiff in l_scurv2:
        if(mdiff > mdiff_c): cut_curv2.append(1);
        else: cut_curv2.append(0);
    t_cut_curv2.append(cut_curv2);
    
# Determine the signal (bb) acceptance and background (se) rejection.
nevts1 = len(l_scurv1);
nevts2 = len(l_scurv2);
sgacc = []; bgrej = [];
sgacc_curv = []; bgrej_curv = [];
i_ccut = 0;
fm = open("{0}/sevsb_tbl.dat".format(plt_base1),"w");
fm.write("# (curv_cut) (se_curv_fpass) (bb_curv_fpass)\n");
for mdiff_c in curv_cuts:

    cut_curv1 = t_cut_curv1[i_ccut];
    cut_curv2 = t_cut_curv2[i_ccut];
    
        
    # Determine the number of signal and background events.
    f1cv = 0.; f2cv = 0.;  # fraction of events passing curv cut
    for pcv1 in cut_curv1:
            
        # Curvature cut
        if(pcv1 > 0): f1cv += 1;
            
    for pcv2 in cut_curv2:
            
        # Curvature cut
        if(pcv2 > 0): f2cv += 1;
        
    # Convert to fractions.
    f1cv /= nevts1;
    f2cv /= nevts2;
        
    # Print to file.
    fm.write("{0} {1} {2}\n".format(mdiff_c,f1cv,f2cv));
        
    # Fill one array of fse vs. fbb for the curvature cuts.
    sgacc_curv.append(1 - f2cv);
    bgrej_curv.append(f1cv);
    
    i_ccut += 1; 

fm.close();

# Output the combined signal acceptance vs. background rejection.
fm = open("{0}/sevsb_combined.dat".format(plt_base1),"w");
fm.write("# (signal_eff) (bg_rejection)\n");
mdist90 = 1.0; s90 = 0.0; b90 = 0.0;       # record the signal acceptance rate for 90% background rejection
mdist80 = 1.0; s80 = 0.0; b80 = 0.0;       # record the signal acceptance rate for 80% background rejection
#for sv,bg in zip(lc_svals,lc_bvals):
for sv,bg in zip(sgacc_curv,bgrej_curv):
    fm.write("{0} {1}\n".format(sv,bg));
    if(abs(bg-0.9) < mdist90):
        mdist90 = abs(bg-0.9);
        s90 = sv;
        b90 = bg;
    if(abs(bg-0.8) < mdist80):
        mdist80 = abs(bg-0.8);
        s80 = sv;
        b80 = bg;
print "*** Curvature only: For bg rejection of {0} we have signal efficiency of {1}".format(b90,s90);
print "*** Curvature only: For bg rejection of {0} we have signal efficiency of {1}".format(b80,s80);
fm.close();

# Make the efficiency vs. background plot.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(bgrej_curv,sgacc_curv,'-.',color='black',label='Curvature');
#plt.plot(lc_svals,lc_bvals,'-',color='black',label='Combined');
lnd = plt.legend(loc=3,frameon=False,handletextpad=0);
plt.xlabel("Background rejection (1-b)");
plt.ylabel("Signal acceptance efficiency ($\epsilon$)");
plt.savefig(fn_plt_sigvsb, bbox_inches='tight');
plt.close();

fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
fmeann, fmeanbins, fmeanpatches = plt.hist(l_scurv1, 40, lw=2, normed=0, histtype='step',color='red',label='Single e-');
rmeann, rmeanbins, rmeanpatches = plt.hist(l_scurv2, 40, normed=0, lw=2, histtype='step',color='blue',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("Curvature asymmetry factor $\phi_C$");
plt.ylabel("Counts/bin");
plt.savefig(fn_plt_curv, bbox_inches='tight');
plt.close();
