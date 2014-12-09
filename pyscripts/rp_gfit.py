"""
rp_gfit.py

Attempt at a global fit to a track.

"""
import sys
import numpy as np
import scipy.integrate as integrate
import random as rd
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from math import *
from scipy.interpolate import interp1d
import scipy.signal as signal

from abc import ABCMeta, abstractmethod

# Input variables
trk_outdir = "/Users/jrenner/IFIC/pylosk/tracks";
trk_name = "magse07_ms0";
num_tracks = 5;
rev_trk = False;
apply_lpf = True;
plt_drawfilter = True;
plt_drawtrk = True;
output_means = False;

if(rev_trk):
    plt_base = "{0}/{1}/plt_rev".format(trk_outdir,trk_name);
else:
    plt_base = "{0}/{1}/plt_fwd".format(trk_outdir,trk_name);
    
if(apply_lpf):
    plt_base = "{0}_flt".format(plt_base);

# Create the plot directory.
if(not os.path.isdir(plt_base)):
    os.mkdir(plt_base);
    print "Creating plot directory {0}...".format(plt_base);

# Remove outliers in a list.
def remove_outliers(lst):
    
    #print "\n --------------- Removing outliers --------------------"; 

    nloop = 0;
    mui = np.mean(lst); sgi = np.std(lst);
    mu = mui; sg = sgi;
    while(nloop == 0 or sg < 0.99*sgi):
    
        nn = 0;
        while(nn < len(lst)):
        
            val = lst[nn];
            if(abs(val-mu) > 2*sg):
                if(nn == 0 or nn == len(lst)-1):
                    lst[nn] = 0;
                else:
                    lst[nn] = (lst[nn-1] + lst[nn+1])/2.;
            
            nn += 1;
    
        mui = mu; sgi = sg;
        mu = np.mean(lst); sg = np.std(lst);
        #print "Mean lst is {0}, stdev lst is {1}".format(mu,sg);

        nloop += 1;

# Design the Butterworth filter
#[od, wn] = signal.buttord([0.05,0.15], [0.08,0.12], 3, 40);
#[bc, ac] = signal.butter(od,wn,'bandstop');
[od, wn] = signal.buttord([0.08,0.12], [0.05,0.15], 3, 40);
[bc, ac] = signal.butter(od,wn,'bandpass');
wfreq, hfreq = signal.freqz(bc, ac);
wfreq_norm = [];
for ww in wfreq:
    wfreq_norm.append(ww/pi);
if(plt_drawfilter):
    
    fig = plt.figure(1);
    fig.set_figheight(5.0);
    fig.set_figwidth(7.5);

    plt.plot(wfreq_norm, abs(hfreq)); #20*np.log10(abs(hfreq)));
    plt.xlabel("Frequency (half cycles/s)");
    plt.ylabel("Frequency Response");
    plt.xscale('log');
    
    plt.savefig("{0}/filter_freq_response.pdf".format(plt_base), bbox_inches='tight');

# Record the values of the curvature.
l_scurv_mean = [];
for trk_num in range(num_tracks):
    
    print "Processing track {0}".format(trk_num);

    # Read in the track.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    if(rev_trk):
        trktbl = np.loadtxt("/Users/jrenner/IFIC/pylosk/tracks/{0}/rev/{1}_{2}.dat".format(trk_name,trk_name,trk_num));
    else:
        trktbl = np.loadtxt("/Users/jrenner/IFIC/pylosk/tracks/{0}/{1}_{2}.dat".format(trk_name,trk_name,trk_num));
    #trktbl = np.loadtxt("data/magse07H_lowerr_b30_f0.dat");
    #trktbl = np.loadtxt("data/mag2e00H_f0.dat");
    #trk_node = trktbl[:,0];
    trk_xM = trktbl[:,0]; #trktbl[:,1];
    trk_yM = trktbl[:,1]; #trktbl[:,2];
    trk_zM = trktbl[:,2]; #trktbl[:,3];
    #trk_xP = trktbl[:,4];
    #trk_yP = trktbl[:,5];
    #trk_zP = trktbl[:,6];
    #trk_chi2P = trktbl[:,7];
    #trk_xF = trktbl[:,8];
    #trk_yF = trktbl[:,9];
    #trk_zF = trktbl[:,10];
    #trk_chi2F = trktbl[:,11];
    
    # Compute the derivatives.
    dxdn = np.diff(trk_xM);
    dydn = np.diff(trk_yM);
    dzdn = np.diff(trk_zM);
    
    d2xdn2 = np.append(np.diff(dxdn),dxdn[-2]-dxdn[-1]);
    d2ydn2 = np.append(np.diff(dydn),dydn[-2]-dydn[-1]);
    d2zdn2 = np.append(np.diff(dzdn),dzdn[-2]-dzdn[-1]);
    
    mu_dxdn = np.mean(dxdn); sg_dxdn = np.std(dxdn);
    mu_dydn = np.mean(dydn); sg_dydn = np.std(dydn);
    mu_dzdn = np.mean(dzdn); sg_dzdn = np.std(dzdn);
    print "Mean dxdn is {0}, stdev dxdn is {1}".format(mu_dxdn,sg_dxdn);
    
    # Remove outliers.
    #remove_outliers(dxdn);
    #remove_outliers(dzdn);
    
    # Calculate first derivatives.
    zz = []; nn = []; dxdz = []; dydz = [];
    n = 0;
    print "z length is {0}, derivative lengths are {1}".format(len(trk_zM[0:-1]),len(dxdn));
    for zval,ddx,ddy,ddz in zip(trk_zM[0:-1],dxdn,dydn,dzdn):
        nn.append(n);
        zz.append(zval);
        if(ddz == 0):
            dxdz.append(0.);
            dydz.append(0.);
        else:
            dxdz.append(ddx/ddz);
            dydz.append(ddy/ddz);
        n += 1;
        
    remove_outliers(dxdz);
    remove_outliers(dydz);
    #print "Len = {0}, {1}".format(len(dxdn),len(d2xdn2));
    
    # Calculate the difference in derivatives.
    diff_dxdy = [];
    for dx,dy in zip(dxdz,dydz):
        diff_dxdy.append(dx-dy);
    
    # Calculate seconds derivatives.
    d2xdz2 = []; d2ydz2 = []; rcurv = []; scurv = [];
    for dxz,dyz,ddz,d2dx,d2dy,d2dz in zip(dxdz,dydz,dzdn,d2xdn2,d2ydn2,d2zdn2):
    
        # 2nd derivatives
        if(ddz == 0):
            d2xz = 0.;
            d2yz = 0.;
        else:
            d2xz = (d2dx - d2dz*dxz)/(ddz**2);
            d2yz = (d2dy - d2dz*dyz)/(ddz**2);
        d2xdz2.append(d2xz);
        d2ydz2.append(d2yz);
    
    # Remove outliers from the second derivatives.
    remove_outliers(d2xdz2);
    remove_outliers(d2ydz2);
    
    # Apply the filters.
    if(apply_lpf):
        f_dxdz = signal.lfilter(bc,ac,dxdz);
        f_dydz = signal.lfilter(bc,ac,dydz);
        f_d2xdz2 = signal.lfilter(bc,ac,d2xdz2);
        f_d2ydz2 = signal.lfilter(bc,ac,d2ydz2);
    else:
        f_dxdz = dxdz;
        f_dydz = dydz;
        f_d2xdz2 = d2xdz2;
        f_d2ydz2 = d2ydz2;
    
    for dxz,dyz,d2xz,d2yz in zip(f_dxdz,f_dydz,f_d2xdz2,f_d2ydz2):
    
        # Radius of curvature and signed curvature.
        if(dxz*d2yz - dyz*d2xz == 0):
            sc = 0.;        
            rc = 0.;
        else:
            sc = (dxz*d2yz - dyz*d2xz)/(dxz**2 + dyz**2)**1.5;
            rc = abs(1.0/sc);
        rcurv.append(rc);
        scurv.append(sc);
        
    # Remove outliers.
    remove_outliers(d2xdz2);
    remove_outliers(d2ydz2);
    remove_outliers(rcurv);
    remove_outliers(scurv);
    
    # Calculate the mean curvature.
    scurv_mean = np.mean(scurv);
    if(trk_zM[-1] < trk_zM[0]): scurv_mean *= -1;    # flip the sign for the opposite z-direction
    print "Mean curvature is {0}".format(scurv_mean);
    l_scurv_mean.append(scurv_mean);
    
    # Compute FFTs
    ffdxdz = np.fft.fft(dxdz);
    fnfreqs = np.fft.fftfreq(len(nn));
    rffdxdz = np.real(ffdxdz);
    #print "Number of samples = {0}".format(len(dxdz));
    
    # Make the plot.
    if(plt_drawtrk):
    
        fig = plt.figure(2);
        fig.set_figheight(15.0);
        fig.set_figwidth(10.0);
        
        ax1 = fig.add_subplot(321);
        ax1.plot(nn,dxdz,'.',color='green')
        ax1.plot(nn,f_dxdz,'-',color='red');
        ax1.set_xlabel("hit number n");
        ax1.set_ylabel("dx/dz");
        
        ax2 = fig.add_subplot(322);
        ax2.plot(fnfreqs,ffdxdz,'-',color='black')
        ax2.set_xlabel("frequency");
        ax2.set_ylabel("FFT(dxdz)");
        
        ax3 = fig.add_subplot(323);
        ax3.plot(nn,dydz,'.',color='green')
        ax3.plot(nn,f_dydz,'-',color='red');
        ax3.set_xlabel("hit number n");
        ax3.set_ylabel("dy/dz");
        
        ax4 = fig.add_subplot(324);
        ax4.plot(nn,scurv,'-',color='green')
        ax4.set_xlabel("hit number n");
        ax4.set_ylabel("signed curvature");
        
        # Create the 3D track plot.
        ax5 = fig.add_subplot(325, projection='3d');
        ax5.plot(trk_xM,trk_yM,trk_zM,'.-');
        ax5.set_xlabel("x (mm)");
        ax5.set_ylabel("y (mm)");
        ax5.set_zlabel("z (mm)");
        
        # Create the x-y projection.
        ax6 = fig.add_subplot(326);
        ax6.plot(trk_xM,trk_yM,'.-');
        ax6.set_xlabel("x (mm)");
        ax6.set_ylabel("y (mm)");
        
        plt.savefig("{0}/plt_signals_{1}.pdf".format(plt_base,trk_num), bbox_inches='tight');
        
        plt.close();
        #plt.show();

#ax3 = fig.add_subplot(313);
#ax3.plot(wfreq_norm, 20*np.log10(abs(hfreq)));

# Output the list of mean curvature values.
if(output_means):
    fm = open("{0}/scurv_means.dat".format(plt_base),"w");
    fm.write("# (trk) (scurv_avg)\n");
    ntrk = 0;
    for scurv in l_scurv_mean:
        fm.write("{0} {1}\n".format(ntrk,scurv));
        ntrk += 1;
    fm.close();