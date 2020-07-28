#!/usr/bin/env python

import numpy as np
import os as os
import subprocess as subproc

intdatapath = 'IntermediateData/'
pavplots = 'PavPlots/'
zfitsfiles = []
cwd = subproc.check_output('pwd',stderr=subproc.STDOUT,shell=True)
print(cwd)

# This program uses psrplot to create a plot similar to pav -FYp
# It saves each file as pgplot.png, which will have to be renamed with every loop

def filefinder(path,ext): # Produces a list of filenames for a given path and extension (eg. 'Data/','.fits')
    outlist = []
    for r, d, f in os.walk(path):
        for file in f:
            if ext in file:
                outlist.append(path+file)
    return outlist

datafiles = filefinder(intdatapath,'.zfits') # Finds the already folded data

for f in datafiles:
    # first we need to get a unique identifier for each obs.
    # the original filename is good for this
    ident = f[:-11] # the original filename minus the last 11 characters
    print(ident) # For error checking
    
    # psrplot can preprocess and plot in one move
    plotstring = 'psrplot -p Y -j Fp ' + str(f)
    print(plotstring) # For error checking
    
    # execute plotting command
    p = subproc.Popen(plotstring,shell=True,stdin=subproc.PIPE,stdout=subproc.PIPE, universal_newlines=True)
    command = '/PNG'
    p.communicate(command)
    
    # rename pgplot.png
    os.rename('pgplot.png',ident+ 'FvPhi' + '.png')
    
# these plots were dumped in the intdatapath, we should move them to their own directory
mvstring = 'mv ' + intdatapath + '*.png ' + pavplots
print(mvstring)
subproc.call(mvstring,shell=True)
    
# indicate completion
print('Success!')
        