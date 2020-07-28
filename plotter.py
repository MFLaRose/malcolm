#!/usr/bin/env python

# import standard packages
import numpy as np
import scipy.optimize as optimize
import scipy.interpolate as interpolate
# import pypulse stuff
import pypulse as pp
import pypulse.archive as arch
import pypulse.utils as u
import pypulse.dynamicspectrum as DS
import pypulse.functionfit as ffit
# import plotting stuff
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# import admin stuff
import datetime as datetime
import os as os
import subprocess as subproc

firsttime = datetime.datetime.now()

threshold = 3.5

def zapbin(dspec, threshold):
    # We get the median and standard deviation of the data
    #threshold -= 0.5
    spec_med = np.median(dspec)
    print(spec_med, 'median dspec value')
    spec_std = np.std(dspec)
    bins = []
    # Now we go and zap bins
    row = 0
    for r in dspec:
        col = 0
        for b in r:
            if np.abs(b-spec_med) > threshold*spec_std:
                dspec[row, col] = 0.0
                bins.append([col, row])
            col += 1
        row += 1
    bins = np.array(bins)
    # And we return the spectra
    #print(bins)
    return dspec

def zap_spectra(array,threshold):
    # This zaps the bad subints
    badsubints = []
    y = [np.std(array[:,i]) for i in range(array.shape[1])]
    for i in range(array.shape[1]):
        if np.abs(y[i] - np.median(y)) > threshold*np.std(y):
            array[:,i] = 0
            badsubints.append(i)
    # This zaps the bad channels
    badfchans = []
    y = [np.std(array[i,:]) for i in range(array.shape[0])]
    for i in range(array.shape[0]):
        if np.abs(y[i] - np.median(y)) > threshold*np.std(y):
            array[i,:] = 0
            badfchans.append(i)
    # Option to save the new dynamic spectrum
    #newdspec = dspecfile.split("b32")[0]+"b32.zapped.txt"
    #np.savetxt(newdspec, array)
    # Now we return the zapped dynamic spectrum
    #print(badsubints,badfchans)
    return array

rawdatapath = 'RawData/'
txtresultspath = 'IntermediateData/'
plotpath = 'Plots/'
parpath = ''

DStxtfiles = []
zfitsfiles = []


for r, d, f in os.walk(txtresultspath):
    for file in f:
        if '.txt' in file:
            DStxtfiles.append(file)
            
for r, d, f in os.walk(rawdatapath):
    for file in f:
        if 'par' in file:
            parpath = rawdatapath+file
            
for r, d, f in os.walk(txtresultspath):
    for file in f:
        if '.zfits' in file:
            zfitsfiles.append(file)
            
for i in range(len(zfitsfiles)):
    
    #dspecdata = np.loadtxt(txtresultspath + zfitsfiles[i])
    zfitsfile = zfitsfiles[i]
    #filename = zfitsfile.replace('.zfits','')
    #filename = filename + 'sFactor'
    sFactor = 1
    while sFactor <= 12:
        filename = zfitsfile.replace('.zfits','')
        filename = filename + 'sFactor'
        filename = filename + str(sFactor) + 'DynSpecTxt.txt'
        dspecdata = np.loadtxt(txtresultspath + filename)
        print(zfitsfile,sFactor)
        

        outstring = 'psrstat -Qq -c '
        dispmeasure = subproc.check_output(outstring + 'dm ' + txtresultspath + zfitsfile, stderr=subproc.STDOUT,shell=True)
        obslen = float(subproc.check_output(outstring + 'length ' + txtresultspath + zfitsfile, stderr=subproc.STDOUT,shell=True))
        numbins = subproc.check_output(outstring + 'nbin ' + txtresultspath + zfitsfile, stderr=subproc.STDOUT,shell=True)
        nfbins = subproc.check_output(outstring + 'nchan ' + txtresultspath + zfitsfile, stderr=subproc.STDOUT,shell=True)
        bandwidth = float(subproc.check_output(outstring + 'bw ' + txtresultspath + zfitsfile, stderr=subproc.STDOUT,shell=True))
        centerfreq = float(subproc.check_output(outstring + 'freq ' + txtresultspath + zfitsfile, stderr=subproc.STDOUT,shell=True))

        minfreq = centerfreq - ((-0.5)*(bandwidth))
        maxfreq = centerfreq + ((-0.5)*(bandwidth))

        obslen_hr = float(obslen)/3600.0
        obslen_min = float(obslen)/60.0

        dspecdata = zapbin(dspecdata,threshold)
        dspecdata = zap_spectra(dspecdata,threshold)


        ds = DS.DynamicSpectrum(dspecdata, \
                                F = np.linspace(minfreq,maxfreq,np.shape(dspecdata)[0]), Funit='MHz',\
                                T = np.linspace(0.0,obslen,np.shape(dspecdata)[1]), Tunit='sec')

        DynSpec = ds.getData()


        fig = plt.figure(figsize=(6,6))
        ax_dspec = fig.add_subplot(111)
        # Plot the Dynamic Spectrum
        DynSpecPlot = ax_dspec.imshow(DynSpec, cmap = cm.hot, aspect = 'auto', interpolation = 'nearest', extent = [0.0, obslen, maxfreq, minfreq], origin = 'lower')
        # Set axis labels
        ax_dspec.set_xlabel('Time [s]')
        ax_dspec.set_ylabel('Frequency [MHz]')

        # Create Colorbar axis
        ax_dspec.axis(xmin=0.0, xmax=obslen, ymin = minfreq, ymax = maxfreq)

        # Create colorbar
        ax_dspec_cbar = plt.colorbar(DynSpecPlot, ax = ax_dspec, fraction = 0.05, pad = 0.03)
        ax_dspec_cbar.set_label(r"Normalized Power")

        # Correct the inverted y-axis
        plt.gca().invert_yaxis()

        # Save the dynamic spectrum with a new extension
        #print(plotpath+filename+'dynSpec.png')
        plt.title('Dynamic Spectrum')
        plt.savefig(plotpath+filename+str(sFactor)+'dynSpec.png',format='png')
        #plt.show()
        plt.close()

     ### Produce secondary spectrum
        secondspec = ds.secondary_spectrum(log=True)

        bw = ds.getBandwidth()
        chunkDuration = ds.dT
        nchans = len(ds.F)
        conj_time_max = np.multiply(np.divide(1.0, np.multiply(2.0, chunkDuration)), \
                                    1000.0)
        conj_time_min = -conj_time_max
        conj_freq_max = np.divide(nchans, np.multiply(2.0, bw))
        conj_freq_min = -conj_freq_max


        fig = plt.figure(figsize=(6,6))
        ax_secondspec = fig.add_subplot(111)
        secspecplot = ax_secondspec.imshow(secondspec/np.max(secondspec), \
                                           cmap = cm.viridis, aspect = 'auto', \
                                           extent=[conj_time_min, conj_time_max, conj_freq_min, conj_freq_max])
        ax_secondspec.set_ylabel(r"Differential Delay (us)")
        ax_secondspec.set_xlabel(r"Differential Doppler Frequency (mHz)")

        ax_secondspec.axis(xmin=conj_time_min, xmax=conj_time_max, ymin=0.0, \
                           ymax=conj_freq_max)

        secspec_cbar = plt.colorbar(secspecplot, ax = ax_secondspec, \
                                    fraction = 0.05, pad = 0.03)
        secspec_cbar.set_label(r"Normalized Power")
        plt.title('Secondary Spectrum')
        plt.tight_layout()
        plt.savefig(plotpath+filename+str(sFactor)+"Secondary.png",format='png')
        plt.close()

        ### Produce ACF
        fig = plt.figure(figsize=(6,6))

        acf2d = ds.acf2d()
        ACF_2D = ds.acf
        plotacf = ds.getACF()
        ax_acf = fig.add_subplot(111) # Trying this...
        acfPlot = ax_acf.imshow(plotacf, cmap=cm.hot, origin = 'lower', aspect = 'auto', \
                                interpolation = 'nearest')
        plt.tight_layout()
        ax_acf.set_xlabel("Time [s]")
        ax_acf.set_ylabel("Frequency [MHz]")
        plt.title('Autocorrelation')
        plt.savefig(plotpath+filename+str(sFactor)+'Autocorrelation.png',format='png')
        #plt.show()
        plt.close()

        # CAN CHANGE DIVISOR FOR BETTER FIT
        print(filename)
        maxr = np.shape(ds.getData())[0]/4 # frequency shape
        maxc = np.shape(ds.getData())[1]/4 # time shape
        #print('Hmm1',np.shape(ds.getData())[1]/1,np.shape(ds.getData())[0]/1)
        delta_t_d, err_t_d, delta_nu_d, err_nu_d \
        = ds.scintillation_parameters(simple=True, full_output=True, show = False, \
                                      maxr=maxr, maxc=maxc, cmap=cm.hot)

        acfparamlist = []
        acfparamlist.append("delta_t_d: "+str(delta_t_d)+"\n")
        acfparamlist.append("err_t_d: "+str(err_t_d)+"\n")
        acfparamlist.append("delta_nu_d: "+str(delta_nu_d)+"\n")
        acfparamlist.append("err_nu_d: "+str(err_nu_d)+"\n")
        #acfparamlist.append("rotation: "+str(rotation)+"\n")
        #acfparamlist.append("err_rot: "+str(err_rot))

        #acfparamarray = np.array(acfparamlist)
        #print(acfparamarray)
        acffile = plotpath+filename+str(sFactor)+'acfparam.txt'

        try:
            subprocess.call('touch '+acffile,shell=True)

        except:
            print("Something's wrong...")

        acffileopen = open(plotpath+filename+str(sFactor)+'acfparam.txt','w')


        acffileopen.writelines(acfparamlist)
        acffileopen.close()

        # Brent's Method
        # Define the initial figure
        fig = plt.figure(figsize = (6,6)) # NEW EDIT HERE
        ax_1Dfit_up = fig.add_subplot(211)
        ax_1Dfit_down = fig.add_subplot(212)
        # Define a function that returns a gaussian, will be use dlater
        def Gfit(x, p):
            return p[0] * np.exp(-((x-p[1])/(np.sqrt(2)*p[2]))**2) + p[3]

        # We need to get a few values from the spectra
        dT = ds.dT # size of time bin, units as input in DS.dynamicspectrum call (here minutes)
        dF = ds.dF # same as dT but the size of the frequency channels (here MHz)
        acfshape = np.shape(ds.acf) # array dimensions of the 2-D ACF
        eta = 0.2 # parameter used to determine scintillation parameters
        NF = len(ds.F) # Number of frequency channels
        Faxis = (np.arange(-(NF-1),NF,dtype=np.float)*np.abs(dF)) # frequency value of each channel
        NT = len(ds.T) # number of time bins
        Taxis = (np.arange(-(NT-1),NT,dtype=np.float)*np.abs(dT))[1:-1] # time value of each time bin
        plotbound = 1.0

        centerrind = np.where(Faxis == 0)[0] # index where frequency lag is 0, or center index in frequency lag
        centercind = np.where(Taxis == 0)[0] # as above but for time lag

        # Now we sum across the frequency axis
        Faxis_summed = np.sum(ds.acf, axis=0) # for timescale
        Faxis_summed = Faxis_summed[1:-1] # cut edges, they're usually bad values
        #print('Lookie!',Faxis_summed)
        # Now sum across the time axis
        Taxis_summed = np.sum(ds.acf, axis=1) # for frequency scale
        FnormFac = np.max(Faxis_summed) # get normalizing factor for fitting later
        TnormFac = np.max(Taxis_summed) # as above
        #print('Lookie2!',Taxis_summed)

        # We change the fitting regions to be the same number of bins as in the 2D fit
        # This can be varied a bit but helps better fit things and give higher consistancy
        low_t_bound = int(centercind-plotbound*maxc+1)
        high_t_bound = int(centercind+plotbound*maxc+1)
        low_f_bound = int(centerrind-plotbound*maxr+1)
        high_f_bound = int(centerrind+plotbound*maxr)
        # Make sure we don't exceed bounds of ACF -> probably not important here, but won't slow things down
        if low_t_bound < 0:
            low_t_bound = 0
        if high_t_bound > len(Taxis):
            high_t_bound = len(Taxis)
        if low_f_bound < 0:
            low_f_bound = 0
        if high_f_bound > len(Faxis):
            high_f_bound = len(Faxis)
        # regardless of bounds, we set the edges of the 1-D ACF to be the minimum value in that range as this
        # helps to better fit a gaussian. This can be played with as well to obtain better fits
        quarter_t_bins = int(np.floor(len(Taxis[low_t_bound:high_t_bound])/8))
        # This is to get a separate array for fitting
        Faxis_cor = []
        for v in Faxis_summed[low_t_bound:high_t_bound]:
            Faxis_cor.append(v)
        Faxis_cor = np.array(Faxis_cor)
            # We replace the values with the appropriate things

        print(low_t_bound,high_t_bound)
        print(Faxis_summed)
        print(Faxis_summed[low_t_bound:high_t_bound]) ### NEW PRINT STATEMENT
        f_rep_val = np.min(Faxis_summed[low_t_bound:high_t_bound])
        Faxis_cor[:quarter_t_bins] = f_rep_val
        Faxis_cor[len(Faxis_summed[low_t_bound:high_t_bound])-quarter_t_bins:] = f_rep_val
        # fit for scintillation timescale along axis better for fitting
        pout, errs = ffit.gaussianfit(Taxis[low_t_bound:high_t_bound],Faxis_cor,baseline=True)
        f = interpolate.interp1d(Taxis,ffit.funcgaussian(pout,Taxis,baseline=True)-(pout[3]+pout[0]/np.e))
        # plot stuff
        # for timescale
        ax_1Dfit_up.plot(Taxis[low_t_bound:high_t_bound], Faxis_summed[low_t_bound:high_t_bound]/FnormFac, c = 'k')
        ax_1Dfit_up.plot(Taxis[low_t_bound:high_t_bound], Gfit(Taxis[low_t_bound:high_t_bound], pout)/FnormFac, c= 'r')
        #ax_1Dfit_up.plot(Taxis[low_t_bound:high_t_bound], Faxis_cor/FnormFac, c = 'g') # this is the line that is fit to
        # This actually finds what the scintillation timescale is, as opposed to just fitting the functions used to find this
        try:
            delta_t_d = optimize.brentq(f,0,Taxis[-1])
            # If this fails we assume the timescale cannot be resolved and set it to the length of the observation
        except:
            print("Error with delta_t_d, setting to max value")
            delta_t_d = dT*NT
        # more plot labeling
        ax_1Dfit_up.set_xlabel(r"Time Lag (s)")
        ax_1Dfit_up.set_ylabel(r"Normalized Power")
        ax_1Dfit_up.set_xlim([Taxis[low_t_bound], Taxis[high_t_bound-1]])
        ax_1Dfit_up.set_ylim([np.min(Faxis_summed[low_t_bound:high_t_bound])/FnormFac-0.05, 1.0])

        # for scintillation bandwidth

        # We replace the edge values with zero to better fit things
        quarter_f_bins = int(np.floor(len(Faxis[low_f_bound:high_f_bound])/8))
        Taxis_cor = []
        for v in Taxis_summed[low_f_bound:high_f_bound]:
            Taxis_cor.append(v)
        Taxis_cor = np.array(Taxis_cor)

        t_rep_val = np.min(Taxis_summed[low_f_bound:high_f_bound])
        Taxis_cor[:quarter_f_bins] = t_rep_val
        Taxis_cor[len(Taxis_summed[low_f_bound:high_f_bound])-quarter_f_bins:] = t_rep_val


        # And do some plotting for this
        ax_1Dfit_down.plot(Faxis[low_f_bound:high_f_bound], Taxis_summed[low_f_bound:high_f_bound]/TnormFac, c = 'k')
        # fit for scintilattion bandwidth
        pout, errs = ffit.gaussianfit(Faxis[low_f_bound:high_f_bound],Taxis_cor,baseline=True)
        f = interpolate.interp1d(Faxis,ffit.funcgaussian(pout,Faxis,baseline=True)-(pout[3]+pout[0]/2))

        ax_1Dfit_down.plot(Faxis[low_f_bound:high_f_bound], Gfit(Faxis[low_f_bound:high_f_bound], pout)/TnormFac, c= 'r')
        #ax_1Dfit_down.plot(Faxis[low_f_bound:high_f_bound], Taxis_cor/TnormFac, c= 'g') # This shows what we used to fit the Gaussian
        # Now actually try to get the scintillation bandwidth
        try:
            delta_nu_d = optimize.brentq(f,0,Faxis[-1])
            # And if it fails assume it's again, unresolvable
        except:
            print "Error with delta_nu_d, setting to max value"
            delta_nu_d = dF*NF
        # more plotting things...
        ax_1Dfit_down.set_xlabel(r"Frequency Lag (MHz)")
        ax_1Dfit_down.set_ylabel(r"Normalized Power")
        ax_1Dfit_down.set_xlim([Faxis[low_f_bound], Faxis[high_f_bound-1]])
        ax_1Dfit_down.set_ylim([np.min(Taxis_summed[low_f_bound:high_f_bound]/TnormFac)-0.05, 1.0])
        # We estimate the errors on the scintillation parameters
        bw = ds.getBandwidth()
        T = ds.getTspan()
        if delta_t_d == 0.0:
            N_d = (1+eta * bw/delta_nu_d)
        elif delta_nu_d == 0.0:
            N_d = (1+eta*T/delta_t_d)
        else:
            N_d = (1+eta * bw/delta_nu_d) * (1+eta*T/delta_t_d)
        fse_nu_d = delta_nu_d/(2*np.log(2)*np.sqrt(N_d)) #log because of FWHM?
        fse_t_d = delta_t_d/(2*np.sqrt(N_d))

        err_nu_d = fse_nu_d
        err_t_d = fse_t_d 

        # And we add these values to the plots to see them easily    
        try:
            scinttxt = r"$\Delta t_{d}$  = %.3f $\pm$ %.3f seconds" % (delta_t_d, err_t_d)
            ax_1Dfit_up.text(0.98, 0.95, scinttxt, ha='right', va='top', color = 'k', \
                             transform=ax_1Dfit_up.transAxes, fontsize=12)
        except:
            print("Error with scintillation timescale")
        try:
            scinttxt = r"$\Delta \nu_{d}$  = %.3f $\pm$ %.3f MHz" % (delta_nu_d, err_nu_d)
            ax_1Dfit_down.text(0.95, 0.95, scinttxt, ha='right', va='top', color = 'k', \
                               transform=ax_1Dfit_down.transAxes, fontsize=12)
        except:
            print("Error with scintillation bandwidth")


        plt.savefig(plotpath+filename+str(sFactor)+'Values.png',format='png')
        plt.close()    
        if sFactor == 1:
            sFactor += 1
        else:
            sFactor += 2