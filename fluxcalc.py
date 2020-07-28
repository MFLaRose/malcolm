#!/usr/bin/env python

import numpy as np
import os as os
import subprocess as subproc

### TO-DO:
# Add function that detects reciever and chooses the right gain and system temp
# Add more comments

### Telescope parameters --> These are specific to Rcvr_342 (300-400 MHz) at the GBT
gain = 2.0 # K/Jy
Tsys = 70.0 # K

### Pulsar parameters --> These are specific to J0038-25
period = 1 # Normalized to one, so pwidth is a fraction between zero and one
pwidth = 0.1 # Fraction of phase with pulse


# This program should be ran from the directory containing the intdatapath
intdatapath = 'IntermediateData/' # Where your folded files are located
zfitsfiles = []
tempfiles = []

# Finds zapped files and temp files if they've already been created
for r, d, f in os.walk(intdatapath):
    for file in f:
        if '.zfits' in file:
            zfitsfiles.append(file)
        if '.temp' in file:
            tempfiles.append(file)
# Skips creaing temp files if they already exist
if len(tempfiles) == 0:
    for i in range(len(zfitsfiles)):
        pamstring = 'pam -e temp -FTp ' + intdatapath + str(zfitsfiles[i]) # Using PSRCHIVE to create temp files
        print(pamstring) 
        try:
            subproc.check_output(pamstring,stderr=subproc.STDOUT,shell=True)
        except subproc.CalledProcessError as e:
            raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd,e.returncode,e.output))
# Skips finding temp files if they already exist
if len(tempfiles) == 0:
    for r, d, f in os.walk(intdatapath):
        for file in f:
            if '.temp' in file:
                tempfiles.append(file)

# Grabs relevant parameters and calculates flux w/ radiometer equation
for i in range(len(tempfiles)):
    statstring = 'psrstat -Qq -c '
    SNR = float(subproc.check_output(statstring + 'snr ' + intdatapath + tempfiles[i], stderr=subproc.STDOUT,shell=True))
    bandwidth = (-float(subproc.check_output(statstring + 'bw ' + intdatapath + tempfiles[i], stderr=subproc.STDOUT,shell=True))-4.0)
    bandwidth = bandwidth*(10**6) # Converts MHz to Hz
    obslen = float(subproc.check_output(statstring + 'length ' + intdatapath + tempfiles[i], stderr=subproc.STDOUT,shell=True))
    npol = float(subproc.check_output(statstring + 'npol ' + intdatapath + tempfiles[i], stderr=subproc.STDOUT,shell=True))
    nbin = float(subproc.check_output(statstring + 'nbin ' + intdatapath + tempfiles[i], stderr=subproc.STDOUT,shell=True))
    onbins = float(subproc.check_output(statstring + 'on:count ' + intdatapath + tempfiles[i], stderr=subproc.STDOUT,shell=True))
    
    pwidth = onbins/nbin # Fraction of period on pulse
    #print(tempfiles[i],SNR,bandwidth,obslen,npol,nbin,onbins,pwidth)
    
    pfactor = np.sqrt(pwidth/(period-pwidth)) # Part of the radiometer equation, done seperately for ease of understanding
    tfactor = (SNR*Tsys)/(gain*np.sqrt(npol*obslen*bandwidth)) # Same as above
    flux = pfactor*tfactor # Gives flux in Jy
    flux = flux*1000.0 # Converts Jy to mJy
    flux = str("{:.4f}".format(flux)) # Converts to a string and gives 4 decimal places of precision
    print(tempfiles[i],flux+ ' mJy') # Displays flux with the file it was calculated from

print('Success!')
        
