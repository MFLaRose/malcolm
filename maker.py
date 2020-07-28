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
import datetime as DT
import os as os
import subprocess as subproc

# Define some useful functions
def Remove(duplicate):
    final_list = []
    for string in duplicate:
        if string not in final_list:
            final_list.append(string)
    return final_list

# Get time at program start
firsttime = DT.datetime.now()

# Directory paths
rawdatapath = 'RawData/'
intdatapath = 'IntermediateData/'
parpath = ''

# List initialization
rawfitsfiles = []
PSRNames = []
MJDs = []
Obsnums = []
Obslist = []
SnippedObsList = []
GUPPIFiles = [] # List of folded files
zapfiles = [] # List of zapped files <-- This is necessary since we rename the files after zapping
FZfiles = [] # List of folded and zapped files used for creating the dynamic spectra



sFactor = 1 # Loop through increasing sFactors for scrunching

for r, d, f in os.walk(rawdatapath):
    for file in f:
        if '.fits' in file:
            rawfitsfiles.append(file)
        if '.par' in file:
            parpath = rawdatapath+file
    
    
for f in rawfitsfiles:
    PSRNames.append(f.split("_")[2]) # Gets PSR name
    MJDs.append(f.split("_")[1]) # Gets obs MJD
    Obsnums.append(f.split("_")[3]) # Gets individual obs number
    f = f[:-10] # This removes the last 10 characters from the filename, which should be the trailing numbers and file extension
    Obslist.append(f) # So we should just have a list of all the files with the ends cut off
    f = f[5:] # Snips off 'guppi_' which is necessary for writing the correct filenames later
    SnippedObsList.append(f) # This should be a list of files looking like MJDXX_JRAXX_DECX_ONUM
    
#print(Obslist)    
PSRNames = Remove(PSRNames)
Obslist = Remove(Obslist) # This takes the list of files and removes duplicates, so we can loop through this correctly with the folding software
SnippedObsList = Remove(SnippedObsList) # And does the same to the above

print('Folding data...')
### BEWARE THAT THIS STEP TAKES AN EXTREMELY LONG TIME
# Changing NumFreqBins or SubIntTime or their while loop conditions can dramatically increase or reduce the amount of files created
# The current values will produce four files for each observation
for i in range(len(Obslist)): # The range(len()) part is necessary to correctly loop through the list here   
    
    NumFreqBins = 512 # This needs to be reinitialized every time, or else we wont loop through numfreqbins
    
    # This loops through number of frequency bins to use while folding
    while NumFreqBins >= 512:
        SubIntTime = 15 # This also needs to be reinitialized every time
        # This loops through the duration of time bins to use while folding
        while SubIntTime >= 15:
            
            foldstring = 'fold_psrfits -o ' + intdatapath + 'GUPPI' + SnippedObsList[i] + '_' + str(NumFreqBins) + 'f' + str(SubIntTime) + 't' + ' -b ' + str(NumFreqBins) + ' -t ' + str(SubIntTime) + ' -S 100 -P ' + parpath + ' ' + rawdatapath+Obslist[i] + '*'
            print(foldstring)
            
            SubIntTime = SubIntTime/2

            # This block of code ouputs the string to the console and raises an error if something goes wrong
            try:
                subproc.check_output(foldstring,stderr=subproc.STDOUT,shell=True)
            except subproc.CalledProcessError as e:
                raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        NumFreqBins = NumFreqBins/2
        
        
# Now that we have folded the data, we need to remove some RFI, or else we won't see the scintillation in the dynamic spectrum
print('Zapping RFI...')

for r, d, f in os.walk(intdatapath):
    for file in f:
        if 'GUPPI' in file:
            GUPPIFiles.append(file)

for z in range(len(GUPPIFiles)):
    zapstring1 = 'paz -e zap -E 2.0 ' + intdatapath + GUPPIFiles[z]
    
    try:
        subproc.check_output(zapstring1,stderr=subproc.STDOUT,shell=True)
    except subproc.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output)) 

# This detects the names of our newly zapped files and adds them to a list
for r, d, f in os.walk(intdatapath):
    for file in f:
        if '.zap' in file:
            zapfiles.append(file)
            
for z in range(len(zapfiles)):
    zapstring2 = "paz -v -m -j 'zap median exp={$off:max-$off:min}, zap median' " + intdatapath + zapfiles[z]
    
    zapstring3 = 'paz -v -m -F "360 380" ' + intdatapath + zapfiles[z]
    zapstring4 = 'paz -v -m -F "794.6 798.6" -F "814.1 820.7" ' + intdatapath + zapfiles[z]
    
    try:
        subproc.check_output(zapstring2,stderr=subproc.STDOUT,shell=True)
        subproc.check_output(zapstring3,stderr=subproc.STDOUT,shell=True)
        subproc.check_output(zapstring4,stderr=subproc.STDOUT,shell=True)
    except subproc.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        
# And scrunch polarization bins to total intensity
for p in range(len(zapfiles)):
    # This string does the scrunching and saves everything as a new file with the .zfits extension
    pamstring = 'pam -e .zfits -p ' + intdatapath + zapfiles[p]
    
    try:
        subproc.check_output(pamstring,stderr=subproc.STDOUT,shell=True)
    except subproc.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

# We need a list of folded and zapped files        
for r, d, f, in os.walk(intdatapath):
    for file in f:
        if '.zfits' in file:
            FZfiles.append(file)
            
# Now we can create the actual dynamic spectrum and scrunch if need be  
print('Creating dynamic spectra...')
for i in range(len(FZfiles)):
    fitsfile = intdatapath+FZfiles[i]
    filename = fitsfile.replace('.zfits','') # This is for giving each created DS a unique filename
    print('Processing ' + fitsfile)
    sFactor = 1
     ### SCRUNCH HERE
    while sFactor <= 12:
        ar = arch.Archive(fitsfile, lowmem = True, baseline_removal = False) # This reinitializes the data we're using so we don't scrunch already scrunched data
        ar = ar.fscrunch(factor=sFactor)
        obslen = ar.getDuration()
        numbins = ar.getNbin()
        nfbins = ar.getNchan()
        #ar = ar.fscrunch(factor=sFactor) Moved to before
        ds = ar.getDynamicSpectrum(windowsize=int(numbins/8),maketemplate=True)
        DynSpec = ds.getData()
        np.savetxt(filename+'sFactor'+str(sFactor)+'DynSpecTxt.txt',DynSpec)
        if sFactor == 1:
            sFactor += 1
        else:
            sFactor += 2

secondtime = DT.datetime.now()
timediff = secondtime - firsttime
print('Runtime: ' + str(timediff))
print('Success!')      
     
    
    
    
    
        




