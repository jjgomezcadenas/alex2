"""
rp_trackplot.py

Plots tracks and their Kalman filter fit
generated from Alex2 and Recpack.

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
from rp_trackdefs import *

from abc import ABCMeta, abstractmethod

# Two main directories: the fit and plot directories.
fnb_trk = "{0}/{1}/trk".format(dat_outdir,run_name);
fnb_plt = "{0}/{1}/plt".format(dat_outdir,run_name);

if(rev_trk):
  trk_name = "{0}_r".format(run_name);
else:
  trk_name = "{0}_f".format(run_name);
    
if(not os.path.isdir(fnb_trk)): 
    print "ERROR: attempting to plot original hits and track directory {0} not available".format(fnb_trk);
    sys.exit();

if(not os.path.isdir(fnb_plt)):
    os.mkdir(fnb_plt);
    print "Creating plot directory {0}...".format(fnb_plt);

# Keep a running list of the values of chi2.
chi2_totlist = [];

# Create num_tracks tracks.
for ntrk in range(num_tracks):
    
    print "-- Plotting track {0}\n".format(ntrk);

    # chi2 lists for this track only
    chi2_list = [];
    chi2_k = [];

    # Read in the track.
    # xM yM zM xP yP zP chi2P xF yF zF chi2F
    trktbl = np.loadtxt("{0}/{1}{2}.dat".format(fnb_trk,trk_name,ntrk));
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

    # Remove outliers from predict and fit lists.
    plt_xP = []; plt_yP = []; plt_zP = []; plt_chi2P = [];
    plt_xF = []; plt_yF = []; plt_zF = []; plt_chi2F = [];   
    for xP,yP,zP,chi2P,xF,yF,zF,chi2F in zip(trk_xP,trk_yP,trk_zP,trk_chi2P,trk_xF,trk_yF,trk_zF,trk_chi2F):

      if(not (xP < -1.0e8 or yP < -1.0e8 or zP < -1.0e8 or chi2P < -1.0e8)):
        plt_xP.append(xP); plt_yP.append(yP); plt_zP.append(zP); plt_chi2P.append(chi2P);
      if(not (xF < -1.0e8 or yF < -1.0e8 or zF < -1.0e8 or chi2F < -1.0e8)):
        plt_xF.append(xF); plt_yF.append(yF); plt_zF.append(zF); plt_chi2F.append(chi2F);
 
    # Construct the chi2 lists.
    for k,chi2 in zip(trk_node,trk_chi2F):

        # Add the chi2 if it is not an outlier.
        if(abs(chi2) < chi2_outlier):
            chi2_list.append(chi2);
            chi2_k.append(k);
        
        # Add to the total chi2 list regardless of the value.
        chi2_totlist.append(chi2);
    
    # Plot the 3-D track plot with projections and chi2.
    if(plt_tracks):

        # Make the plot.
        fig = plt.figure(1);
        fig.set_figheight(5.0);
        fig.set_figwidth(5.0);

        # Create the 3D track plot.
        ax1 = fig.add_subplot(111, projection='3d');
        ax1.plot(trk_xM,trk_yM,trk_zM,'.',color='black');
        if(plt_filtered): ax1.plot(plt_xF,plt_yF,plt_zF,'-',color='black');
        if(plt_prediction): ax1.plot(plt_xP,plt_yP,plt_zP,'.',color='blue');
        xst, xen = ax1.get_xlim(); ax1.xaxis.set_ticks(np.arange(xst,xen,(xen-xst)/5.));
        yst, yen = ax1.get_ylim(); ax1.yaxis.set_ticks(np.arange(yst,yen,(yen-yst)/5.));
        zst, zen = ax1.get_zlim(); ax1.zaxis.set_ticks(np.arange(zst,zen,(zen-zst)/5.));
        ax1.set_xlabel("x ({0})".format(plt_units),fontsize=10.);
        ax1.set_ylabel("y ({0})".format(plt_units),fontsize=10.);
        ax1.set_zlabel("z ({0})".format(plt_units),fontsize=10.);
        lb_x = ax1.get_xticklabels();
        lb_y = ax1.get_yticklabels();
        lb_z = ax1.get_zticklabels();
        for lb in (lb_x + lb_y + lb_z):
            lb.set_fontsize(7.);

        # Print the 3D plot.
        if(plt_3dprint):
            fn_plt = "{0}/plt3d_{1}_{2}.pdf".format(fnb_plt,trk_name,ntrk);
            plt.savefig(fn_plt, bbox_inches='tight');
        
        # Make the plot.
        fig = plt.figure(2);
        fig.set_figheight(15.0);
        fig.set_figwidth(10.0);
        
        # Create the 3D track plot.
        ax1 = fig.add_subplot(321, projection='3d');
        ax1.plot(trk_xM,trk_yM,trk_zM,'.',color='black');
        if(plt_filtered): ax1.plot(plt_xF,plt_yF,plt_zF,'-',color='black');
        if(plt_prediction): ax1.plot(plt_xP,plt_yP,plt_zP,'.',color='blue');
        xst, xen = ax1.get_xlim(); ax1.xaxis.set_ticks(np.arange(xst,xen,(xen-xst)/5.));
        yst, yen = ax1.get_ylim(); ax1.yaxis.set_ticks(np.arange(yst,yen,(yen-yst)/5.));
        zst, zen = ax1.get_zlim(); ax1.zaxis.set_ticks(np.arange(zst,zen,(zen-zst)/5.)); 
        ax1.set_xlabel("x ({0})".format(plt_units));
        ax1.set_ylabel("y ({0})".format(plt_units));
        ax1.set_zlabel("z ({0})".format(plt_units));
        lb_x = ax1.get_xticklabels();
        lb_y = ax1.get_yticklabels();
        lb_z = ax1.get_zticklabels();
        for lb in (lb_x + lb_y + lb_z):
            lb.set_fontsize(8.);
        
        # Create the x-y projection.
        ax2 = fig.add_subplot(322);
        ax2.plot(trk_xM,trk_yM,'.',color='black');
        if(plt_filtered): ax2.plot(plt_xF,plt_yF,'-',color='black');
        if(plt_prediction): ax2.plot(plt_xP,plt_yP,'.',color='blue');
        ax2.set_xlabel("x ({0})".format(plt_units));
        ax2.set_ylabel("y ({0})".format(plt_units));
        
        # Create the x-z projection.
        ax3 = fig.add_subplot(323);
        ax3.plot(trk_xM,trk_zM,'.',color='black');
        if(plt_filtered): ax3.plot(plt_xF,plt_zF,'-',color='black');
        if(plt_prediction): ax3.plot(plt_xP,plt_zP,'.',color='blue');
        ax3.set_xlabel("x ({0})".format(plt_units));
        ax3.set_ylabel("z ({0})".format(plt_units));    
        
        # Create the y-z projection.
        ax4 = fig.add_subplot(324);
        ax4.plot(trk_yM,trk_zM,'.',color='black');
        if(plt_filtered): ax4.plot(plt_yF,plt_zF,'-',color='black');
        if(plt_prediction): ax4.plot(plt_yP,plt_zP,'.',color='blue');
        ax4.set_xlabel("y ({0})".format(plt_units));
        ax4.set_ylabel("z ({0})".format(plt_units));
        
        # Plot the chi2.
        ax5 = fig.add_subplot(325);
        ax5.plot(chi2_k,chi2_list,color='black');
        ax5.set_xlabel("k");
        ax5.set_ylabel("$\chi^{2}$");
        #ax5.set_yscale("log");
        ax5.set_title("Avg. chi2f = {0}".format(np.mean(plt_chi2F)));

        # Show and/or print the plot.
        if(plt_print):
            fn_plt = "{0}/plt_{1}_{2}.pdf".format(fnb_plt,trk_name,ntrk);
            plt.savefig(fn_plt, bbox_inches='tight');
        if(plt_show):
            plt.show();
            
        plt.close();
            
# Plot the chi2 histogram.
if(plt_chi2):
    fig = plt.figure(4);
    fig.set_figheight(5.0);
    fig.set_figwidth(7.5);
    chi2n, chi2bins, chi2patches = plt.hist(chi2_list, 100, normed=0, histtype='step',color='blue',label='Forward fit');
    plt.xlabel("$\chi^{2}$");
    plt.ylabel("Counts/bin");
    
    # Show and/or print the plot.
    if(plt_print):
        fn_plt = "{0}/chi2_{1}.pdf".format(fnb_plt,trk_name);
        plt.savefig(fn_plt, bbox_inches='tight');
    if(plt_show):
        plt.show();
        
    plt.close();
