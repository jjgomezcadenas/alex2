# -*- coding: utf-8 -*-
"""
rp_cfmeans.py

Compares the mean curvatures across two runs.

@author: josh
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import *

# Get the track names as arguments.
args = sys.argv;
if(args[1] != "" and args[2] != ""):
    trk_name1 = args[1];
    trk_name2 = args[2];
else:
    trk_name1 = "magseHtest01_1mm_5sp";   # the single-electron event run
    trk_name2 = "magbbHtest01_1mm_5sp";   # the double-beta event run

# Input variables
#trk_outdir = "/Users/jrenner/IFIC/pylosk/tracks";
#trk_outdir = "/Users/jrenner/IFIC/software/alex2/build/alexMain/out";
trk_outdir = "/data4/NEXT/users/jrenner/kalmanfilter/alex2/build/alexMain/out";

plt_base1 = "{0}/{1}".format(trk_outdir,trk_name1);
plt_base2 = "{0}/{1}".format(trk_outdir,trk_name2);
apply_lpf = True;

meandiff_cuts = np.arange(-0.03,0.03,0.06/100);

# Read in the mean signed curvature values.
if(apply_lpf):
    scurvtbl1 = np.loadtxt("{0}/plt_fwd_flt/scurv_means.dat".format(plt_base1));
    scurvtbl2 = np.loadtxt("{0}/plt_fwd_flt/scurv_means.dat".format(plt_base2));
else:
    scurvtbl1 = np.loadtxt("{0}/plt_fwd/scurv_means.dat".format(plt_base1));
    scurvtbl2 = np.loadtxt("{0}/plt_fwd/scurv_means.dat".format(plt_base2));
    
l_scurv1 = scurvtbl1[:,1];
l_time1 = scurvtbl1[:,6];

l_scurv2 = scurvtbl2[:,1];
l_time2 = scurvtbl2[:,6];

# Prepare the plotted quantities.
snp_diff1 = [];
for npos,nneg in zip(l_npos1_f,l_nneg1_f):
    snp_diff1.append(npos-nneg);
    
snp_diff2 = [];
for npos,nneg in zip(l_npos2_f,l_nneg2_f):
    snp_diff2.append(npos-nneg);

# Create the histograms.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
fmeann, fmeanbins, fmeanpatches = plt.hist(l_scurv1, 40, linestyle='--', normed=0, histtype='step',color='black',label='Single e-');
rmeann, rmeanbins, rmeanpatches = plt.hist(l_scurv2, 40, linestyle='--', normed=0, histtype='step',color='black',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("Curvature Asymmetry Factor");
plt.ylabel("Counts/bin");

# Print the plot.
if(apply_lpf):
    fn_plt = "{0}/scurv_diff_means_lpf.pdf".format(plt_base1);
else:
    fn_plt = "{0}/scurv_diff_means.pdf".format(plt_base1);

plt.savefig(fn_plt, bbox_inches='tight');
plt.close();

fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
fnpn1, fnpbins1, fnppatches1 = plt.hist(snp_diff1, 40, normed=0, histtype='step',color='black',label='Single e-');
fnpn2, fnpbins2, fnppatches2 = plt.hist(snp_diff2, 40, normed=0, histtype='step',color='black',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("(Positive - Negative) Track Segments");
plt.ylabel("Counts/bin");
    
# Print the plot.
if(apply_lpf):
    fn_plt = "{0}/npnn_diff_lpf.pdf".format(plt_base1);
else:
    fn_plt = "{0}/npnn_diff.pdf".format(plt_base1);
    
plt.savefig(fn_plt, bbox_inches='tight');
plt.close();

# Time histograms
fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
fnpn1, fnpbins1, fnppatches1 = plt.hist(l_time1, 800, linestyle='--', normed=0, histtype='step',color='black',label='Single e-');
fnpn2, fnpbins2, fnppatches2 = plt.hist(l_time2, 400, normed=0, histtype='step',color='black',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
axes = plt.gca()
axes.set_xlim([0.1e-9,1.0e-8])
plt.xlabel("Time Elapsed");
plt.ylabel("Counts/bin");

# Print the plot.
if(apply_lpf):
    fn_plt = "{0}/timehist_lpf.pdf".format(plt_base1);
else:
    fn_plt = "{0}/timehist.pdf".format(plt_base1);
    
plt.savefig(fn_plt, bbox_inches='tight');
plt.close();
