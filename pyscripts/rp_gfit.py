"""
rp_gfit.py

Attempt at a global fit to a track.

"""
import sys
import numpy as np
import scipy.integrate as itg
import random as rd
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from math import *
from scipy.interpolate import interp1d
import scipy.signal as signal

from abc import ABCMeta, abstractmethod

# Get the arguments, if any.
args = sys.argv;
if(args[1] != ""):
    trk_name = args[1];
else:
    trk_name = "magseHtest01_1mm_5sp";

# Input variables
#trk_outdir = "/Users/jrenner/IFIC/pylosk/tracks";
trk_outdir = "/Users/jrenner/IFIC/software/alex2/build/alexMain/out";
num_tracks = 5;
rev_trk = False;
apply_lpf = True;
plt_drawfilter = False;
plt_drawtrk = True;
output_means = False;

brk_win = 5;
#brk_ns = 2;
#Ttrack = -1.;

# Constants
KEinit = 2.447;   # initial electron kinetic energy (for a single-electron background event)
Pgas = 10.;       # gas pressure in atm
Tgas = 293.15;    # gas temperature in Kelvin
Bfield = 0.5;     # applied magnetic field in Tesla

pc_rho0 = 2.6867774e19;   # density of ideal gas at T=0C, P=1 atm in cm^(-3)
pc_m_Xe = 131.293;        # mass of xenon in amu
pc_NA = 6.02214179e23;    # Avogadro constant
pc_eC = 1.602176487e-19;  # electron charge in C
pc_me = 9.10938215e-31;   # electron mass in kg
pc_clight = 2.99792458e8; # speed of light in m/s

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
            if(abs(val-mu) > 5*sg):
                if(nn == 0 or nn == len(lst)-1):
                    lst[nn] = 0;
                else:
                    lst[nn] = (lst[nn-1] + lst[nn+1])/2.;
            
            nn += 1;
    
        mui = mu; sgi = sg;
        mu = np.mean(lst); sg = np.std(lst);
        #print "Mean lst is {0}, stdev lst is {1}".format(mu,sg);

        nloop += 1;
        
# Find the vertex given the curvature profile.
def calc_vertex(cprof):
    
    nvert = 0;
    npts = len(cprof);
    sum_max = 0;
    
    # Attempt to maximize the integral of:
    #  (-profile) from 0 to n (+profile) from n+1 to npts-1 (inclusive)
    for nn in range(npts):
        
        # Calculate the vertex sum.
        lsum = 0; rsum = 0;
        if(nn != 0):
            lsum = sum(cprof[0:nn]);
        if(nn < npts):
            rsum = sum(cprof[nn:npts]);
        cs = rsum - lsum;
        
        # Update the vertex if the sum
        if(nn == 0 or cs > sum_max):
            nvert = nn;
            sum_max = cs;
        
        #print "-- Point nn = {0} gives sum = {1}".format(nn,cs);
            
    # Return the vertex.
    return (1.0*nvert/npts);

# Find the breakpoints (potential vertices) given a list of second derivatives.
def find_brkpts(d2vals):
    
    # Find the mean and sigma of the given list.
    mu = np.mean(d2vals);
    sgm = sqrt(np.var(d2vals));
    
    # Subtract the mean.
    d2vals_mu0 = [];
    for vv in d2vals: d2vals_mu0.append(vv);#-mu);
        
    # Determine the breakpoints.
    nn = 0;
    brkpts = [];
    while(nn < len(d2vals_mu0)):
        
        # If we've reached a potential breakpoint, only identify one point
        #  in the peak.
        if(abs(d2vals_mu0[nn]) > 5.): #brk_ns*sgm):
            
            npk = nn;
            npk_val = d2vals_mu0[nn];
            nn += 1;
            
            # Ensure only one point for every brk_win points.
            while(nn < len(d2vals_mu0) and (nn-npk) < brk_win):
                
                if(abs(d2vals_mu0[nn]) > abs(npk_val)):
                    npk = nn;
                    npk_val = d2vals_mu0[nn];
        
                nn += 1;
            
            # Add the breakpoint.
            brkpts.append(npk);
            
        else:
            nn += 1;


    # Return the list of breakpoints.
    return brkpts;
    
# Analyze the breakpoints and curvature profile to determine whether this
#  track belongs to a single-electron.
def analyze_se(brkprof,curvprof):
    
    # Single-electron flag to be returned.
    se = 0;
    
    # Segment lists
    l_seg_i = [];
    l_seg_len = [];
    l_seg_sgn = [];

    # Loop quantities
    ibrk = 0;
    iseg = 0;
    sgn_sum = 0;
    nn = 0;
    
    # Perform the segment analysis.
    while(nn < len(curvprof)):
        
        sgn_sum += curvprof[nn];
        
        if(nn == len(curvprof)-1 or (ibrk < len(brkprof) and nn == brkprof[ibrk])):

            # Calculate the quantities of interest.
            slen = nn-iseg;
            ssgn = 1;
            if(slen > 0 and sgn_sum/slen < 0): ssgn = -1;
            
            # Append to the segment lists.
            l_seg_i.append(iseg);
            l_seg_len.append(slen);
            l_seg_sgn.append(ssgn);

            # Update loop quantities.
            ibrk += 1;
            iseg = nn;
            sgn_sum = 0;
        
        nn += 1;

    print "\n\n*** Segment analysis ***";
    print "Initial index:";
    print l_seg_i;
    print "Length:";
    print l_seg_len;
    print "Sign:";
    print l_seg_sgn;    
    
    # Determine whether or not this is a single-electron event.
    npseg = 0; nnseg = 0;
    for ii,ll,ss in zip(l_seg_i, l_seg_len, l_seg_sgn):
        
        if(ll >= 0 and ss >= 0 and 1.0*(ii+ll)/len(curvprof) > 0.5): npseg += ll;
        if(ll >= 0 and ss < 0 and 1.0*ii/len(curvprof) <= 0.5): nnseg += ll;
    
    print "{0} positive segments and {1} negative segments".format(npseg,nnseg);
    return npseg, nnseg;

# Calculate the time required to create the track.
def calc_track_time():

    # Prepare the stopping power table.    
    xesp_tbl = np.loadtxt("data/xe_estopping_power_NIST.dat");
    rho = pc_rho0*(Pgas/(Tgas/273.15))*(pc_m_Xe/pc_NA);
    e_vals = xesp_tbl[:,0];
    dEdx_vals = xesp_tbl[:,1];
    e_vals = np.insert(e_vals,0,0.0);
    dEdx_vals = np.insert(dEdx_vals,0,dEdx_vals[0]);
    xesp = interp1d(e_vals,dEdx_vals*rho,kind='cubic');

    # Compute the integral.
    tval = itg.quad(lambda x: (sqrt((x+0.511)**2-0.511**2)/(x+0.511)*xesp(x))**(-1),0.001,KEinit,limit=500)[0];
    tval = tval/(100*pc_clight);
    
    return tval;

# Calculate the average track time if it has not been set.
#if(Ttrack < 0):
Ttrack = calc_track_time();
print " *** Average time to create track is {0}".format(Ttrack);

# Calculate the cyclotron frequency.
wcyc = -1.0*(pc_eC/pc_me)*Bfield;

# Record the values of the mean curvature and fraction of positive curvature.
l_scurv_mean = []; l_pcurv_frac = []; l_scurv_vertex = []; l_npseg = []; l_nnseg = [];
l_time = [];
for trk_num in range(num_tracks):

#    if(rev_trk):
#        trk_file = "/Users/jrenner/IFIC/pylosk/tracks/{0}/rev/{1}_{2}.dat".format(trk_name,trk_name,trk_num);
#    else:
    #trk_file = "/Users/jrenner/IFIC/pylosk/tracks/{0}/{1}_{2}.dat".format(trk_name,trk_name,trk_num);
    trk_file = "{0}/{1}/toyMC/{2}/{3}_{4}.dat".format(trk_outdir,trk_name,trk_name,trk_name,trk_num);

    if(not os.path.isfile(trk_file)):
        print "Note: did not find file {0}".format(trk_file);
        continue;

    # Read in the track.
    trktbl = np.loadtxt(trk_file);
    # xM yM ziM zfM ux uy uz E deltaE deltaX
    #trktbl = np.loadtxt("data/magse07H_lowerr_b30_f0.dat");
    #trktbl = np.loadtxt("data/mag2e00H_f0.dat");
    trk_xM = trktbl[:,0]; #trktbl[:,1];
    trk_yM = trktbl[:,1]; #trktbl[:,2];
    trk_zM = trktbl[:,2]; #trktbl[:,3];
    trk_ux = trktbl[:,4];
    trk_uy = trktbl[:,5];
    trk_uz = trktbl[:,6];
    trk_eM = trktbl[:,7];
    trk_dxM = trktbl[:,9];
    
    print "Processing track {0} with {1} samples".format(trk_num,len(trk_xM));

    # Calculate the time elapsed over the creation of the track.
    telapsed = 0.;
    for ee,dx in zip(trk_eM,trk_dxM):
        if(ee > 0.):
            vel = sqrt((ee+0.511)**2-0.511**2)/(ee+0.511)*pc_clight*1000;
            telapsed += dx/vel;
    print "Actual time elapsed = {0}".format(telapsed);
    l_time.append(telapsed);
    if(abs(telapsed-Ttrack) > Ttrack):
        print "*** WARNING: Elapsed time not as expected... using computed average time.\n\n";
        telapsed = Ttrack;
    
    # Design the Butterworth filter.
    NN = len(trk_xM);
    fsamp = NN/telapsed;
    fcbar = abs(wcyc)/(2*pi*fsamp);
    #[od, wn] = signal.buttord(1.0*fcbar, 1.5*fcbar, 3, 20);
    [od, wn] = signal.buttord(1.5*fcbar, 2.0*fcbar, 3, 20);
    print "fcbar = {0}".format(fcbar);

    [bc, ac] = signal.butter(od,wn,'lowpass');
    #print "Filter coefficients are: b = ";
    #print bc;
    #print "And a = ";
    #print ac;
    
    print "Designing filter.";
    
    # -------------------------------------------------------------------------
    # FIR filter from: http://wiki.scipy.org/Cookbook/FIRFilter
    # The Nyquist rate of the signal.
    nyq_rate = 1.0 / 2.0;  # assuming sampling freq = 1
    
    # The desired width of the transition from pass to stop,
    # relative to the Nyquist rate.
    width = 0.1 / nyq_rate;
    
    # The desired attenuation in the stop band, in dB.
    ripple_db = 40.0;

    # Compute the order and Kaiser parameter for the FIR filter.
    N, beta = signal.kaiserord(ripple_db, width);

    # The cutoff frequency of the filter.
    cutoff_freq = 0.5*fcbar;

    # Use firwin with a Kaiser window to create a lowpass FIR filter.
    taps = signal.firwin(N, cutoff_freq/nyq_rate, window=('kaiser', beta), nyq=nyq_rate);
    #print taps

    # Use lfilter to filter x with the FIR filter.
    #filtered_x = lfilter(taps, 1.0, x);

    # Compute the filter delay.
    fdelay = int(N/2);

    # -------------------------------------------------------------------------

    print "Drawing filter.";
    
    # Draw the filter if the option is set to do so.
    fnfreqs = np.fft.fftfreq(len(trk_xM));
    wfreq, hfreq = signal.freqz(bc, ac, fnfreqs*pi);
    wfreq_norm = [];
    for ww in wfreq:
        wfreq_norm.append(ww/pi);
    if(plt_drawfilter):
    
        fig = plt.figure(1);
        fig.set_figheight(5.0);
        fig.set_figwidth(7.5);

        plt.plot(wfreq_norm, abs(hfreq)); #20*np.log10(abs(hfreq)));
        plt.xlabel("Frequency (half cycles/sample)");
        plt.ylabel("Frequency Response");
        plt.xscale('log');
    
        plt.savefig("{0}/filter_freq_response.pdf".format(plt_base), bbox_inches='tight');
        plt.close();
        
        fig = plt.figure(2);
        plt.clf();
        w, h = signal.freqz(taps, worN=8000);
        plt.plot(w/pi, np.absolute(h), linewidth=2);
        plt.xlabel('Frequency (half cycles/sample)');
        plt.ylabel('Gain');
        plt.title('Frequency Response');
        plt.ylim(-0.05, 1.05);
        plt.grid(True);
        
        plt.savefig("{0}/FIR_freq_resp.pdf".format(plt_base), bbox_inches='tight');
        plt.close();
    
    # Calculate the discrete cyclotron frequency.
    ux0 = trk_ux[0]; uy0 = trk_uy[0]; print "ux = {0}, uy = {1}".format(ux0,uy0);
    fsamp = len(trk_xM)/telapsed;
    fcbar = abs(wcyc)/(2*pi*fsamp);
    rcyc = sqrt((KEinit+0.511)**2-0.511**2)/(KEinit+0.511)*sqrt(ux0**2 + uy0**2)*pc_clight*100/abs(wcyc);
    print "Sampling frequency = {0} samp/s, cyclotron freq = {1} cyc/s, fcbar = {2}, max. rcyc = {3} cm; should see {4} cycles".format(fsamp,wcyc/(2*pi),fcbar,rcyc,abs(telapsed*wcyc/(2*pi)));
    
    # Mean subtract and save the original unfiltered track.
#    trk_xM_0 = []; trk_yM_0 = []; trk_zM_0 = [];
#    muX = np.mean(trk_xM); muY = np.mean(trk_yM); muZ = np.mean(trk_zM);
#    for xx in trk_xM: trk_xM_0.append(xx-muX);
#    for yy in trk_yM: trk_yM_0.append(yy-muY);
#    for zz in trk_zM: trk_zM_0.append(zz-muZ);

    # Reverse the hits if this is a reverse fit.
    if(rev_trk):
        trk_xM_0 = trk_xM[::-1];
        trk_yM_0 = trk_yM[::-1];
        trk_zM_0 = trk_zM[::-1];
    else:
        trk_xM_0 = trk_xM;
        trk_yM_0 = trk_yM;
        trk_zM_0 = trk_zM;
    
    # Apply the filters.
    if(apply_lpf):
        trk_xM = signal.lfilter(taps,1.0,trk_xM_0);
        trk_yM = signal.lfilter(taps,1.0,trk_yM_0);
        trk_zM = signal.lfilter(taps,1.0,trk_zM_0);
        #trk_xM = signal.lfilter(bc,ac,trk_xM_0);
        #trk_yM = signal.lfilter(bc,ac,trk_yM_0);
        #trk_zM = signal.lfilter(bc,ac,trk_zM_0);
        print "Applying filter with coefficients:";
        print taps;
        
        # Correct for the FIR delay.
        print "Applying delay of {0} samples".format(fdelay);
        trk_xM = np.roll(trk_xM,-1*fdelay);
        trk_yM = np.roll(trk_yM,-1*fdelay);
        trk_zM = np.roll(trk_zM,-1*fdelay);
    
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
    remove_outliers(dxdn);
    remove_outliers(dydn);
    remove_outliers(dzdn);
    
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
    d2xdz2 = []; d2ydz2 = []; rcurv = []; scurv = []; sscurv = [];
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
    
    # Determine the breakpoints (potential vertices) from the second derivative.
    brk_pts = find_brkpts(d2xdz2);
    #print "Found breakpoints:";
    #print brk_pts;
    #print d2xdz2;
    
    # Remove outliers from the second derivatives.
    #remove_outliers(d2xdz2);
    #remove_outliers(d2ydz2);
    
    for dzn,dxz,dyz,d2xz,d2yz in zip(dzdn,dxdz,dydz,d2xdz2,d2ydz2):
    
        # Radius of curvature and signed curvature.
        if(dxz*d2yz - dyz*d2xz == 0):
            sc = 0.;
            rc = 0.;
        else:
            sc = (dxz*d2yz - dyz*d2xz)/(dxz**2 + dyz**2)**1.5;
            if(dzn < 0): sc *= -1;  # correct for the direction of travel
            rc = abs(1.0/sc);
        rcurv.append(rc);
        scurv.append(sc);
        if(sc > 0): sscurv.append(1);
        else: sscurv.append(-1);
        
    # Remove outliers.
    remove_outliers(rcurv);
    remove_outliers(scurv);
    
    # Analyze the track curvature profile based on these breakpoints.
    # Want segments with decided sign, length, and location.
    # Single-electron condition is having beginning and end segments long enough.
    npseg,nnseg = analyze_se(brk_pts,sscurv);
    
    # Calculate the mean curvature.
    #scurv_mean = np.mean(scurv[0:len(scurv)/5]);
    #scurv_mean = np.mean(scurv[4*len(scurv)/5:]);
    #scurv_mean = np.mean(scurv[3*len(scurv)/5:5*len(scurv)/5]);
    scurv_mean = 1.0*(np.mean(sscurv[0:len(sscurv)/2]) - np.mean(sscurv[len(sscurv)/2:]))/len(sscurv);

    # Calculate the vertex.
    scurv_vertex = calc_vertex(sscurv);
    
    # Determine the fraction of total points with positive curvature.
    npos = 0;
    for cval in scurv:
        if(cval > 0):
            npos += 1;
    pfrac = 1.0*npos/len(sscurv);
    
    # Correct for tracks traveling in the -z direction.
#    if(trk_zM[-1] < trk_zM[0]): 
#        scurv_mean *= -1;    # flip the sign for the opposite z-direction
#        pfrac = 1 - pfrac;   # reverse the fraction
        
    print "Mean curvature is {0}; positive fraction is {1}; vertex is {2}".format(scurv_mean,pfrac,scurv_vertex);
    l_scurv_mean.append(scurv_mean);
    l_pcurv_frac.append(pfrac);
    l_scurv_vertex.append(scurv_vertex);
    if(npseg + nnseg > 0): l_npseg.append(1.0*npseg/(npseg+nnseg));
    else: l_npseg.append(0.);
    if(npseg + nnseg > 0): l_nnseg.append(1.0*nnseg/(npseg+nnseg));
    else: l_nnseg.append(0.);
    #l_scurv_vertex_seg.append(scurv_vertex_seg);
    
    # Compute FFTs
    ffxvals = np.fft.fft(trk_xM_0);
    #rffdxdz = np.real(ffdxdz);
    #print "Number of samples = {0}".format(len(dxdz));
    #print fnfreqs;
    
    # Make the plot.
    if(plt_drawtrk):
        
        # Create the lists of (x,y,z) for positive and negative curvature.
        pc_x = []; pc_y = []; pc_z = [];
        nc_x = []; nc_y = []; nc_z = [];
        for cv,xx,yy,zz in zip(sscurv, trk_xM, trk_yM, trk_zM):
            if(cv > 0):
                pc_x.append(xx);
                pc_y.append(yy);
                pc_z.append(zz);
            else:
                nc_x.append(xx);
                nc_y.append(yy);
                nc_z.append(zz);
        
        fig = plt.figure(2);
        fig.set_figheight(15.0);
        fig.set_figwidth(10.0);
        
        ax1 = fig.add_subplot(321);
        #ax1.plot(nn,dxdz,'.-',color='red');
        ax1.plot(nn,trk_xM_0[0:-1],'.-',color='red');
        ax1.plot(nn,trk_xM[0:-1],'.-',color='red');
        ax1.set_xlabel("hit number n");
        ax1.set_ylabel("dx/dz");
        
        ax2 = fig.add_subplot(322);
#        ax2.plot(nn,dydz,'.',color='green')
#        ax2.plot(nn,f_dydz,'-',color='red');
        ax2.plot(nn,d2xdz2,'.-',color='green');
        ax2.set_ylim([-15.,15.]);
        ax2.set_xlabel("hit number n");
        ax2.set_ylabel("d$^2$x/dz$^2$");
        
        ax3 = fig.add_subplot(323); #fig = plt.figure(1);
        scurvn, scurvbins, scurvpatches = ax3.hist(scurv, 40, normed=0, histtype='step',color='blue');
        ax3.set_xlabel("signed curvature");
        ax3.set_ylabel("Counts/bin");
        
        ax3 = fig.add_subplot(323);
        ax3.plot(fnfreqs,abs(ffxvals),'.',color='black',label='FFT');
        ax3.plot(wfreq_norm, 20*abs(hfreq),label='LP filter x20');
        lnd = plt.legend(loc=1,frameon=False,handletextpad=0,fontsize=8);
        ax3.set_ylim([0,100]);
        ax3.set_xlabel("frequency");
        ax3.set_ylabel("FFT(x)");
        #print fnfreqs;
        #print ffxvals;
        
        ax4 = fig.add_subplot(324);
        ax4.plot(nn,sscurv,'-',color='green')
        ax4.set_ylim([-1.5,1.5]);
        ax4.set_xlabel("hit number n");
        ax4.set_ylabel("sign of curvature");
        ax4.set_title("mean sgn(curvature) = {0}".format(round(np.mean(sscurv),4)));
        
        # Create the 3D track plot.
        ax5 = fig.add_subplot(325, projection='3d');
        ax5.plot(trk_xM,trk_yM,trk_zM,'-',color='black');
        ax5.plot(pc_x,pc_y,pc_z,'+',color='red');
        ax5.plot(nc_x,nc_y,nc_z,'.',color='blue');
        ax5.set_xlabel("x (mm)");
        ax5.set_ylabel("y (mm)");
        ax5.set_zlabel("z (mm)");
        
        lb_x = ax5.get_xticklabels();
        lb_y = ax5.get_yticklabels();
        lb_z = ax5.get_zticklabels();
        for lb in (lb_x + lb_y + lb_z):
            lb.set_fontsize(8);  
        
        # Create the x-y projection.
        ax6 = fig.add_subplot(326);
        ax6.plot(trk_xM,trk_yM,'-',color='black');
        ax6.plot(pc_x,pc_y,'+',color='red');
        ax6.plot(nc_x,nc_y,'.',color='blue');        
        ax6.set_xlabel("x (mm)");
        ax6.set_ylabel("y (mm)");
        
        plt.savefig("{0}/plt_signals_{1}.pdf".format(plt_base,trk_num), bbox_inches='tight');
        
        plt.close();
        #plt.show();
        
    print "\n";

# Output the list of mean curvature values and fraction of positive curvature values.
if(output_means):
    fm = open("{0}/scurv_means.dat".format(plt_base),"w");
    fm.write("# (trk) (scurv_avg) (pos_frac) (vertex frac.) (npseg) (nnseg) (time)\n");
    ntrk = 0;
    for scurv,nfpos,vert,npseg,nnseg,tval in zip(l_scurv_mean,l_pcurv_frac,l_scurv_vertex,l_npseg,l_nnseg,l_time):
        fm.write("{0} {1} {2} {3} {4} {5} {6}\n".format(ntrk,scurv,nfpos,vert,npseg,nnseg,tval));
        ntrk += 1;
    fm.close();