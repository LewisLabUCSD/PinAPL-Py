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
SeqQCDir = config['SeqQCDir']
LogFileDir = config['LogFileDir']
IndexScript = config['IndexScript']
LoaderScript = config['LoaderScript']
ReadDepthScript = config['ReadDepthScript']
SeqQCScript = config['SeqQCScript']
AlignScript = config['AlignScript']
NormalizeScript = config['NormalizeScript']
StatsScript = config['StatsScript']
ControlScript = config['ControlScript']
sgRNARankScript = config['sgRNARankScript']
GeneRankScript = config['GeneRankScript']
ScatterScript = config['ScatterScript']
ReplicateScript = config['ReplicateScript']
ClusterScript = config['ClusterScript']

# Make config file accessible
os.system('cp configuration.yaml '+ScriptsDir)
os.chdir(ScriptsDir)

# Print Header
os.system('python -u PrintStatus.py Header blank 2>&1 | tee PinAPL-Py.log')
start = time.time()

# Generate index if not present
if not os.path.exists(IndexDir):
    StatMsg = 'Building library index ...'
    os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+IndexScript+'.py 2>&1 | tee -a PinAPL-Py.log')
    DoneMsg = 'Library index completed.'
    os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')    

# Read Samples
StatMsg = 'Reading data ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+LoaderScript+'.py 2>&1 | tee -a PinAPL-Py.log')
SampleNames, Treatments, Replicates = GetSamples()
DoneMsg = '\nData input completed.'    
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Run Sequence Quality Control
StatMsg = 'Running sequence quality control ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
if os.path.exists(SeqQCDir):
    os.system('python -u PrintStatus.py SkipSeqQC blank 2>&1 | tee -a PinAPL-Py.log')
else: 
    os.system('python -u '+SeqQCScript+'.py 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Sequence quality check completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')    
    
# Align Reads
StatMsg = 'Aligning reads ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    if os.path.exists(AlignDir+sample) and os.path.exists(AlnQCDir+sample):
        os.system('python -u PrintStatus.py SkipSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    else:        
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        os.system('python -u '+AlignScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')        
os.system('python -u '+ReadDepthScript+'.py 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read alignments completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Normalize Counts
StatMsg = 'Normalizing read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+NormalizeScript+'.py'+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Normalization completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Analyze Counts
StatMsg = 'Analyzing read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+StatsScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read count analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Analyze Controls
StatMsg = 'Analyzing controls ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+ControlScript+'.py'+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Control analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Find sgRNA hits
StatMsg = 'sgRNA '+ScreenType+' analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames: 
    if 'Control' not in sample:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        os.system('python -u '+sgRNARankScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'sgRNA '+ScreenType+' analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare sample scatterplots
StatMsg = 'Read count scatterplots ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    if 'Control' not in sample:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')        
        os.system('python -u '+ScatterScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read count scatterplots completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare replicate scatterplots
StatMsg = 'Replicate correlation ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for treatment in Treatments:
    os.system('python -u PrintStatus.py ProcessSample '+treatment+' 2>&1 | tee -a PinAPL-Py.log')
    if len(Replicates[treatment]) >= 2: 
        PairIt = itertools.combinations(Replicates[treatment],2)
        for pair in PairIt:
            os.system('python -u '+ReplicateScript+'.py '+pair[0]+' '+pair[1]+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = '\nReplicate correlation analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare heatmap
StatMsg = 'Clustering analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+ClusterScript+'.py'+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Clustering analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Rank genes
StatMsg = 'Gene ranking analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    if 'Control' not in sample:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
        os.system('python -u '+GeneRankScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Gene ranking analysis completed.'
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
    
# Move Log Files 
if not os.path.exists(LogFileDir):
	os.makedirs(LogFileDir) 
os.system('cp configuration.yaml '+LogFileDir)
os.system('cp PinAPL-Py.log '+LogFileDir)
os.chdir(WorkingDir)
os.system('cp DataSheet.xlsx '+LogFileDir)