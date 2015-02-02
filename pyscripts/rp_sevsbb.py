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
trk_outdir = "/Users/jrenner/IFIC/software/alex2/build/alexMain/out";

plt_base1 = "{0}/{1}".format(trk_outdir,trk_name1);
plt_base2 = "{0}/{1}".format(trk_outdir,trk_name2);
apply_lpf = True;

fn_plt_sigvsb = "{0}/sigvsb_all.pdf".format(plt_base1);
fn_plt_chi2rdiff = "{0}/chi2rdiff.pdf".format(plt_base1);
fn_plt_curvchi2 = "{0}/curvchi2_scatter.pdf".format(plt_base1);
fn_plt_curvchi2s = "{0}/curvchi2_s_scatter.pdf".format(plt_base1);
fn_plt_curvchi2b = "{0}/curvchi2_b_scatter.pdf".format(plt_base1);

# Number of points in combined curve.
ncpts = 50;

# Read in the mean signed curvature values.
if(apply_lpf):
    scurvtbl1_f = np.loadtxt("{0}/plt_fwd_flt/scurv_means.dat".format(plt_base1));
    scurvtbl2_f = np.loadtxt("{0}/plt_fwd_flt/scurv_means.dat".format(plt_base2));
else:
    scurvtbl1_f = np.loadtxt("{0}/plt_fwd/scurv_means.dat".format(plt_base1));
    scurvtbl2_f = np.loadtxt("{0}/plt_fwd/scurv_means.dat".format(plt_base2));

l_scurv1_f = scurvtbl1_f[:,1];
l_pfrac1_f = scurvtbl1_f[:,2];
l_npos1_f = scurvtbl1_f[:,4];
l_nneg1_f = scurvtbl1_f[:,5];

l_scurv2_f = scurvtbl2_f[:,1];
l_pfrac2_f = scurvtbl2_f[:,2];
l_npos2_f = scurvtbl2_f[:,4];
l_nneg2_f = scurvtbl2_f[:,5];
    
# Read in the computed chi2 values.
chi2tbl1 = np.loadtxt("{0}/plt/chi2values.dat".format(plt_base1));
chi2tbl2 = np.loadtxt("{0}/plt/chi2values.dat".format(plt_base2));

l_fchi2F1 = chi2tbl1[:,0];
l_fchi2R1 = chi2tbl1[:,1];
l_rchi2F1 = chi2tbl1[:,2];
l_rchi2R1 = chi2tbl1[:,3];

l_fchi2F2 = chi2tbl2[:,0];
l_fchi2R2 = chi2tbl2[:,1];
l_rchi2F2 = chi2tbl2[:,2];
l_rchi2R2 = chi2tbl2[:,3];

# Create the quantities of interest for the chi2 cut.
l_chi2rf1 = []; l_chi2rr1 = []; l_chi2rd1 = [];
for fcF,fcR,rcF,rcR in zip(l_fchi2F1,l_fchi2R1,l_rchi2F1,l_rchi2R1):
    l_chi2rf1.append(fcF/fcR);
    l_chi2rr1.append(rcF/rcR);
    l_chi2rd1.append(rcF/rcR - fcF/fcR);
l_chi2rf2 = []; l_chi2rr2 = []; l_chi2rd2 = [];
for fcF,fcR,rcF,rcR in zip(l_fchi2F2,l_fchi2R2,l_rchi2F2,l_rchi2R2):
    l_chi2rf2.append(fcF/fcR);
    l_chi2rr2.append(rcF/rcR);
    l_chi2rd2.append(rcF/rcR - fcF/fcR);
    
# Cut ranges.
chi2_cuts = np.arange(-3.,4.,7.0/100);     # cut on difference in chi2 ratios
curv_cuts = np.arange(-0.008,0.008,0.016/100); # cut on curvature asymmetry values

# Determine if a forward fit was found.
# 0: no forward fits
# 1: forward fit was found to be the forward-most fit
# 2: the reverse fit was found to be the forward-most fit
t_cut_chi21 = []; t_cut_chi22 = [];
for chi2r_c in chi2_cuts:
    
    cut_chi21 = []; cut_chi22 = [];
    
    for chi2rd in l_chi2rd1:
        
        if(chi2rd > chi2r_c): cut_chi21.append(1);
        else: cut_chi21.append(0);
        
    for chi2rd in l_chi2rd2:
        
        if(chi2rd > chi2r_c): cut_chi22.append(1);
        else: cut_chi22.append(0);
            
    t_cut_chi21.append(cut_chi21);
    t_cut_chi22.append(cut_chi22);

# Determine which events passed the curvature cut.
# 0: would not pass the curvature cut
# 1: passes the curvature cut with currently selected orientation
# 2: passes the curvature cut if the sign of the asymmetry factor is flipped
t_cut_curv1 = []; t_cut_curv2 = [];
for mdiff_c in curv_cuts:
    
    cut_curv1 = []; cut_curv2 = [];
    
    for mdiff in l_scurv1_f:
        if(mdiff > mdiff_c): cut_curv1.append(1);
        else: cut_curv1.append(0);
    t_cut_curv1.append(cut_curv1);
    
    for mdiff in l_scurv2_f:
        if(mdiff > mdiff_c): cut_curv2.append(1);
        else: cut_curv2.append(0);
    t_cut_curv2.append(cut_curv2);
    
# Create the arrays to make the combined signal vs. background curve.
lc_bvals = []; lc_svals = [];
for b in range(ncpts):
    lc_svals.append(1.0*b/ncpts);
    lc_bvals.append(0);

# Determine the signal (bb) acceptance and background (se) rejection.
nevts1 = len(chi2tbl1);
nevts2 = len(chi2tbl2);
sgacc = []; bgrej = [];
sgacc_curv = []; bgrej_curv = [];
sgacc_chi2 = []; bgrej_chi2 = [];
i_ccut = 0;
fm = open("{0}/sevsb_tbl.dat".format(plt_base1),"w");
fm.write("# (chi2_cut) (curv_cut) (se_chi2_fpass) (se_curv_fpass) (bb_chi2_fpass) (bb_curv_fpass) (se_both_pass) (bb_both_pass)\n");
for mdiff_c in curv_cuts:

    cut_curv1 = t_cut_curv1[i_ccut];
    cut_curv2 = t_cut_curv2[i_ccut];
    
    i_c2cut = 0;
    for chi2rd_c in chi2_cuts:

        cut_chi21 = t_cut_chi21[i_c2cut];
        cut_chi22 = t_cut_chi22[i_c2cut];
        
        # Determine the number of signal and background events.
        f1c2 = 0.; f2c2 = 0.;  # fraction of events passing chi2 cut
        f1cv = 0.; f2cv = 0.;  # fraction of events passing curv cut
        f1a = 0.; f2a = 0.;    # fraction of events passing both cuts
        for pc21,pcv1 in zip(cut_chi21,cut_curv1):
            
            # chi2 cut
            if(pc21 > 0): f1c2 += 1;
            
            # Curvature cut
            if(pcv1 > 0): f1cv += 1;
            
            # Record if the event has passed all cuts.
            if(pc21 > 0 and pcv1 > 0): f1a += 1;
            
        for pc22,pcv2 in zip(cut_chi22,cut_curv2):
            
            # chi2 cut
            if(pc22 > 0): f2c2 += 1;
            
            # Curvature cut
            if(pcv2 > 0): f2cv += 1;
            
            # Record if the event has passed all cuts.
            if(pc22 > 0 and pcv2 > 0): f2a += 1;
        
        # Convert to fractions.
        f1c2 /= nevts1;
        f1cv /= nevts1;
        f2c2 /= nevts2;
        f2cv /= nevts2;
        f1a /= nevts1;
        f2a /= nevts2;
        
        # Print to file.
        fm.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(chi2rd_c,mdiff_c,f1c2,f1cv,f2c2,f2cv,f1a,f2a));
        
        # Update the combined curve.
        b = int(f1a*ncpts);
        if(b == ncpts): b = ncpts-1;
        fbb_check = lc_bvals[b];
        if(fbb_check < (1 - f2a)): lc_bvals[b] = 1 - f2a;
        
        # Fill one array of fse vs. fbb for the chi2 cuts.
        if(i_ccut == 0):
            sgacc_chi2.append(1 - f2c2);
            bgrej_chi2.append(f1c2);
        
        i_c2cut += 1;
    
    # Fill one array of fse vs. fbb for the curvature cuts.
    sgacc_curv.append(1 - f2cv);
    bgrej_curv.append(f1cv);
    
    i_ccut += 1; 
        
#        fbg_chi2 = 1.0*sum(cut_chi21)/nevts1;
#        fsg_chi2 = 1.0 - 1.0*sum(cut_chi22)/nevts2;
#        
#        fbg_curv = 1.0*sum(cut_curv1)/nevts1;
#        fsg_curv = 1.0 - 1.0*sum(cut_curv2)/nevts2;

fm.close();

# Make the efficiency vs. background plot.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(bgrej_chi2,sgacc_chi2,'-',color='blue',label='Mult. scattering');
plt.plot(bgrej_curv,sgacc_curv,'-',color='green',label='Curvature');
plt.plot(lc_svals,lc_bvals,'-',color='black',label='Combined');
lnd = plt.legend(loc=3,frameon=False,handletextpad=0);
plt.xlabel("Background rejection (1-b)");
plt.ylabel("Signal acceptance efficiency ($\epsilon$)");
plt.savefig(fn_plt_sigvsb, bbox_inches='tight');
plt.close();

fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
ffracn, ffracbins, ffracpatches = plt.hist(l_chi2rd1, 40, normed=0, histtype='step',color='blue',label='Single e-');
rfracn, rfracbins, rfracpatches = plt.hist(l_chi2rd2, 40, normed=0, histtype='step',color='red',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("$(\chi^2_F/\chi^2_R)_{\mathrm{reverse}} - (\chi^2_F/\chi^2_R)_{\mathrm{forward}}$");
plt.ylabel("Counts/bin");
plt.savefig(fn_plt_chi2rdiff, bbox_inches='tight');
plt.close();

fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(l_scurv1_f, l_chi2rd1,'.',color='blue',label='Single electrons');
plt.plot(l_scurv2_f, l_chi2rd2,'.',color='green',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("curvature asymmetry factor");
plt.ylabel("chi2 difference");
plt.savefig(fn_plt_curvchi2, bbox_inches='tight');
plt.close();

fig = plt.figure(4);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(l_scurv1_f, l_chi2rd1,'.',color='blue',label='Single electrons');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("curvature asymmetry factor");
plt.ylabel("chi2 difference");
plt.savefig(fn_plt_curvchi2s, bbox_inches='tight');
plt.close();

fig = plt.figure(5);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(l_scurv2_f, l_chi2rd2,'.',color='green',label='0vbb');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("curvature asymmetry factor");
plt.ylabel("chi2 difference");
plt.savefig(fn_plt_curvchi2b, bbox_inches='tight');
plt.close();