#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 19:40:43 2016

@author: philipp
"""

# PinAPL-Py Master File
# =======================================================================

import os
import yaml
import itertools
import time
import re
import pandas
from PrintStatus import *
from ReadDataSheet import *

configFile = open('configuration.yaml','r')
config = yaml.load(configFile)
configFile.close()    
ScreenType = config['ScreenType']
WorkingDir = config['WorkingDir']
LibDir = config['LibDir']
IndexDir = config['IndexDir']
DataDir = config['DataDir']
ScriptsDir = config['ScriptsDir']
AlignDir = config['AlignDir']
AlnQCDir = config['AlnQCDir']
script00 = config['script00']
script01 = config['script01']
script02 = config['script02']
script03 = config['script03']
script04 = config['script04']
script05 = config['script05']
script06 = config['script06']
script07 = config['script07']
script08 = config['script08']
script09 = config['script09']
script10 = config['script10']
script11 = config['script11']

# Print Header
os.system('python -u PrintStatus.py Header blank 2>&1 | tee PinAPL-Py.log')
start = time.time()

# Generate index if not present
if not os.path.exists(IndexDir):
    StatMsg = 'Building library index ...'
    os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+script00+'.py 2>&1 | tee -a PinAPL-Py.log')
    DoneMsg = 'Library index completed.'
    os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')    

# Read Samples
os.chdir(ScriptsDir)
StatMsg = 'Reading data ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+script01+'.py 2>&1 | tee -a PinAPL-Py.log')
SampleNames, Treatments, Replicates = GetSamples()
os.chdir(ScriptsDir)
os.system('python -u '+script02+'.py 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = '\nData input completed.'    
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')
    
# Align Reads
os.chdir(ScriptsDir)
StatMsg = 'Aligning reads ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    if os.path.exists(AlignDir+sample) and os.path.exists(AlnQCDir+sample):
        os.system('python -u PrintStatus.py SkipSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    else:        
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        os.system('python -u '+script03+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        
DoneMsg = 'Read alignments completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Normalize Counts
os.chdir(ScriptsDir)
StatMsg = 'Normalizing read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+script04+'.py'+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Normalization completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Analyze Counts
os.chdir(ScriptsDir)
StatMsg = 'Analyzing read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+script05+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read count analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Analyze Controls
os.chdir(ScriptsDir)
StatMsg = 'Analyzing controls ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+script06+'.py'+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Control analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Find sgRNA hits
os.chdir(ScriptsDir)
StatMsg = 'sgRNA '+ScreenType+' analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames: 
    if 'Control' not in sample:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        os.system('python -u '+script07+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'sgRNA '+ScreenType+' analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Rank genes
os.chdir(ScriptsDir)
StatMsg = 'Gene ranking analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    if 'Control' not in sample:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        os.system('python -u '+script08+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Gene ranking analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare sample scatterplots
os.chdir(ScriptsDir)
StatMsg = 'Read count scatterplots ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    if 'Control' not in sample:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')        
        os.system('python -u '+script09+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read count scatterplots completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare replicate scatterplots
os.chdir(ScriptsDir)
StatMsg = 'Replicate correlation ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for treatment in Treatments:
    os.system('python -u PrintStatus.py ProcessSample '+treatment+' 2>&1 | tee -a PinAPL-Py.log')
    if len(Replicates[treatment]) >= 2: 
        PairIt = itertools.combinations(Replicates[treatment],2)
        for pair in PairIt:
            os.system('python -u '+script10+'.py '+pair[0]+' '+pair[1]+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = '\nReplicate correlation analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare heatmap
os.chdir(ScriptsDir)
StatMsg = 'Clustering analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+script11+'.py'+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Clustering analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Final Time Stamp
end = time.time()
os.system('python -u PrintStatus.py AllDone blank 2>&1 | tee -a PinAPL-Py.log') 
sec_elapsed = end - start
if sec_elapsed < 60:
    time_elapsed = sec_elapsed
    StatMsg = 'Time elapsed [secs]: ' + '%.3f' % time_elapsed +'\n'
    os.system('python -u PrintStatus.py TimeStamp "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
elif sec_elapsed < 3600:
    time_elapsed = sec_elapsed/60
    StatMsg = 'Time elapsed [mins]: ' + '%.3f' % time_elapsed +'\n'
    os.system('python -u PrintStatus.py TimeStamp "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
else:
    time_elapsed = sec_elapsed/3600
    StatMsg = 'Time elapsed [hours]: ' + '%.3f' % time_elapsed +'\n'
    os.system('python -u PrintStatus.py TimeStamp "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
    
# Move Log File  
os.system('mv PinAPL-Py.log '+WorkingDir)