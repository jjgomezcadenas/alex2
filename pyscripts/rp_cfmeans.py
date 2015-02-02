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
    trk_name1 = "magseHtest01_1mm_5sp"; # "magse05H_le_5sp_lseg_pfit";  # the single-electron event run
    trk_name2 = "magbbHtest01_1mm_5sp"; # "magbb05H_le_5sp_lseg_pfit";  # the double-beta event run

# Input variables
#trk_outdir = "/Users/jrenner/IFIC/pylosk/tracks";
trk_outdir = "/Users/jrenner/IFIC/software/alex2/build/alexMain/out";

plt_base1 = "{0}/{1}".format(trk_outdir,trk_name1);
plt_base2 = "{0}/{1}".format(trk_outdir,trk_name2);
apply_lpf = True;

meandiff_cuts = np.arange(-0.03,0.03,0.06/100);

# Read in the mean signed curvature values.
if(apply_lpf):
    scurvtbl1_f = np.loadtxt("{0}/plt_fwd_flt/scurv_means.dat".format(plt_base1));
#    scurvtbl1_r = np.loadtxt("{0}/plt_rev_flt/scurv_means.dat".format(plt_base1));
    scurvtbl2_f = np.loadtxt("{0}/plt_fwd_flt/scurv_means.dat".format(plt_base2));
#    scurvtbl2_r = np.loadtxt("{0}/plt_rev_flt/scurv_means.dat".format(plt_base2));
else:
    scurvtbl1_f = np.loadtxt("{0}/plt_fwd/scurv_means.dat".format(plt_base1));
#    scurvtbl1_r = np.loadtxt("{0}/plt_rev/scurv_means.dat".format(plt_base1));
    scurvtbl2_f = np.loadtxt("{0}/plt_fwd/scurv_means.dat".format(plt_base2));
#    scurvtbl2_r = np.loadtxt("{0}/plt_rev/scurv_means.dat".format(plt_base2));
    
l_scurv1_f = scurvtbl1_f[:,1];
l_pfrac1_f = scurvtbl1_f[:,2];
l_npos1_f = scurvtbl1_f[:,4];
l_nneg1_f = scurvtbl1_f[:,5];
l_time1_f = scurvtbl1_f[:,6];
#l_scurv1_r = scurvtbl1_r[:,1];
#l_pfrac1_r = scurvtbl1_r[:,2];
#l_npos1_r = scurvtbl1_r[:,4];
#l_nneg1_r = scurvtbl1_r[:,5];

l_scurv2_f = scurvtbl2_f[:,1];
l_pfrac2_f = scurvtbl2_f[:,2];
l_npos2_f = scurvtbl2_f[:,4];
l_nneg2_f = scurvtbl2_f[:,5];
l_time2_f = scurvtbl2_f[:,6];
#l_scurv2_r = scurvtbl2_r[:,1];
#l_pfrac2_r = scurvtbl2_r[:,2];
#l_npos2_r = scurvtbl2_r[:,4];
#l_nneg2_r = scurvtbl2_r[:,5];

# Prepare the plotted quantities.
#scurv_diff1 = [];
#for scurvf,scurvr in zip(l_scurv1_f,l_scurv1_r):
#    scurv_diff1.append(scurvf - scurvr);
#
#scurv_diff2 = [];
#for scurvf,scurvr in zip(l_scurv2_f,l_scurv2_r):
#    scurv_diff2.append(scurvf - scurvr);
    
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
fmeann, fmeanbins, fmeanpatches = plt.hist(l_scurv1_f, 40, normed=0, histtype='step',color='blue',label='Single e-');
rmeann, rmeanbins, rmeanpatches = plt.hist(l_scurv2_f, 40, normed=0, histtype='step',color='red',label='0vbb');
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
fnpn1, fnpbins1, fnppatches1 = plt.hist(snp_diff1, 40, normed=0, histtype='step',color='blue',label='Single e-');
fnpn2, fnpbins2, fnppatches2 = plt.hist(snp_diff2, 40, normed=0, histtype='step',color='red',label='0vbb');
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

# Determine the signal (bb) acceptance and background (se) rejection.
sgacc = []; bgrej = []; 
for mdiff_c in meandiff_cuts:

    btemp = 0;
    for mdiff in l_scurv1_f:
        if(mdiff < mdiff_c): btemp += 1;
    btemp = 1.0-1.0*btemp/len(l_scurv1_f);
    
    stemp = 0;
    for mdiff in l_scurv2_f:
        if(mdiff < mdiff_c): stemp += 1;
    stemp = 1.0*stemp/len(l_scurv2_f); 
    
    sgacc.append(stemp);
    bgrej.append(btemp);
    
    print "Stemp = {0}, btemp = {1}".format(stemp,btemp);

fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(sgacc,bgrej,'-',color='black');
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("Background rejection (r)");
plt.ylabel("Single-electron acceptance efficiency ($\epsilon$)");
    
# Print the plot.
if(apply_lpf):
    fn_plt = "{0}/sigvsbg_lpf.pdf".format(plt_base1);
else:
    fn_plt = "{0}/sigvsbg.pdf".format(plt_base1);
    
plt.savefig(fn_plt, bbox_inches='tight');
plt.close();

# Time histograms
fig = plt.figure(4);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
fnpn1, fnpbins1, fnppatches1 = plt.hist(l_time1_f, 800, normed=0, histtype='step',color='blue',label='Single e-');
fnpn2, fnpbins2, fnppatches2 = plt.hist(l_time2_f, 400, normed=0, histtype='step',color='red',label='0vbb');
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