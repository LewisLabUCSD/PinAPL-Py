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
import glob
import pandas
from PrintStatus import *
from ReadDataSheet import *
from ExtractTop10Genes import *

# Open configuration file
configFile = open('configuration.yaml','r')
config = yaml.load(configFile)
configFile.close()  
# Get parameters 
ScreenType = config['ScreenType']
IncludeDensityPlots = config['IncludeDensityPlots']
IncludeGeneRankCombination = config['IncludeGeneRankCombination']
AutoHighlight = config['AutoHighlight']
# Get directories
WorkingDir = config['WorkingDir']
LibDir = config['LibDir']
IndexDir = config['IndexDir']
DataDir = config['DataDir']
ScriptsDir = config['ScriptsDir']
AlignDir = config['AlignDir']
AlnQCDir = config['AlnQCDir']
SeqQCDir = config['SeqQCDir']
LogFileDir = config['LogFileDir']
sgRNAReadCountDir = config['sgRNAReadCountDir']
GeneDir = config['GeneDir']
# Get scripts
SanityScript = config['SanityScript']
IndexScript = config['IndexScript']
LoaderScript = config['LoaderScript']
ReadDepthScript = config['ReadDepthScript']
SeqQCScript = config['SeqQCScript']
TrimScript = config['TrimScript']
AlignScript = config['AlignScript']
ClassifyScript = config['ClassifyScript']
CutoffScript = config['CutoffScript']
CleanUpScript = config['CleanUpScript']
NormalizeScript = config['NormalizeScript']
StatsScript = config['StatsScript']
ControlScript = config['ControlScript']
sgRNARankScript = config['sgRNARankScript']
zFCScript = config['zFCScript']
vFCScript = config['vFCScript']
AverageCountsScript = config['AverageCountsScript']
GeneRankScript = config['GeneRankScript']
GenePlotScript = config['GenePlotScript']
CombineScript = config['CombineScript']
ScatterScript = config['ScatterScript']
DensityScript = config['DensityScript']
ReplicateScript = config['ReplicateScript']
ClusterScript = config['ClusterScript']
ExtractTop10Script = config['ExtractTop10Script']

# Check for error file
if os.path.exists('ErrorFound.txt'):
    os.remove('ErrorFound.txt')

# Make config file accessible
os.system('cp configuration.yaml '+ScriptsDir)
os.chdir(ScriptsDir)

# Print Header
os.system('python -u PrintStatus.py Header blank 2>&1 | tee PinAPL-Py.log')
start = time.time()

# Character sanity check
StatMsg = 'Running character sanity check ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+SanityScript+'.py 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Character sanity check completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

## Generate index if not present
#if not os.path.exists(IndexDir):
#    StatMsg = 'Building library index ...'
#    os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
#    os.system('python -u '+IndexScript+'.py 2>&1 | tee -a PinAPL-Py.log')
#    DoneMsg = 'Library index completed.'
#    os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')    

# Read Samples
StatMsg = 'Reading sample definition input ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+LoaderScript+'.py 2>&1 | tee -a PinAPL-Py.log')
SampleNames, Treatments, TreatmentSamples, Replicates = GetSamples()
DoneMsg = 'Sample definition completed.'    
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

## Run Sequence Quality Control
#StatMsg = 'Running sequence quality control ...'
#os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
#if os.path.exists(SeqQCDir):
#    os.system('python -u PrintStatus.py SkipSeqQC blank 2>&1 | tee -a PinAPL-Py.log')
#else: 
#    os.system('python -u '+SeqQCScript+'.py 2>&1 | tee -a PinAPL-Py.log')
#DoneMsg = 'Sequence quality control completed.'
#os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')    
#  
## Trim Adapters
#StatMsg = 'Trimming reads ...'
#os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
#if os.path.exists(AlignDir):
#    os.system('python -u PrintStatus.py SkipTrim blank 2>&1 | tee -a PinAPL-Py.log')
#else:
#    os.system('python -u '+TrimScript+'.py 2>&1 | tee -a PinAPL-Py.log')
#DoneMsg = 'Read trimming completed.'
#os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')
#  
## Align Reads
#StatMsg = 'Mapping reads to library ...'
#os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
#NewAlignmentDone = False
#for sample in SampleNames:
#    if os.path.exists(AlignDir+sample):
#        os.system('python -u PrintStatus.py SkipAlignment '+sample+' 2>&1 | tee -a PinAPL-Py.log')
#    else:        
#        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
#        os.system('python -u '+AlignScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')     
#        NewAlignmentDone = True
#DoneMsg = 'Read alignments completed.'
#os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')        

# Get sgRNA and Gene Read Counts
StatMsg = 'Collecting read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    # Analyze alignment output and get sgRNA counts
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    #os.system('python -u '+ClassifyScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    # Apply sgRNA count cutoff and get gene counts
    os.system('python -u '+CutoffScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read count acquisition completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')     

## Delete temporary alignment output
#if NewAlignmentDone:
#    StatMsg = 'Clearing temporary alignment files ...'
#    os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
#    os.system('python -u '+CleanUpScript+'.py 2>&1 | tee -a PinAPL-Py.log')
#    DoneMsg = 'Alignment output processing completed.'
#    os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

## Show read depth
#os.system('python -u '+ReadDepthScript+'.py 2>&1 | tee -a PinAPL-Py.log')

# Average counts over replicates
StatMsg = 'Averaging read counts over replicates ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for treatment in Treatments:
    os.system('python -u '+AverageCountsScript+'.py '+treatment+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Read count averaging completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Normalize Counts
StatMsg = 'Normalizing read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+NormalizeScript+'.py'+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Read count normalization completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Update sample names (to include the averaged samples)
os.chdir(sgRNAReadCountDir)
Filenames = glob.glob('*GuideCounts.txt')
SampleNames = [filename[0:-16] for filename in Filenames]
TreatmentSamples = [samplename for samplename in SampleNames if 'Control' not in samplename]
os.chdir(ScriptsDir)

# Analyze Counts
StatMsg = 'Analyzing sgRNA read count distribution ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in SampleNames:
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+StatsScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Read count distribution analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Analyze Controls
StatMsg = 'Analyzing control samples ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+ControlScript+'.py'+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Control sample analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Rank sgRNA
StatMsg = 'sgRNA '+ScreenType+' analysis ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in TreatmentSamples: 
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+sgRNARankScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'sgRNA '+ScreenType+' analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare sgRNA scatterplots
StatMsg = 'Plotting sgRNA read counts ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in TreatmentSamples:
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')        
    os.system('python -u '+ScatterScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+zFCScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')    
    os.system('python -u '+vFCScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')    
DoneMsg = 'Read count scatterplots completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare density plots
if IncludeDensityPlots:
    StatMsg = 'Plotting sgRNA read counts densities ...'
    os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
    for sample in TreatmentSamples:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')        
        os.system('python -u '+DensityScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    DoneMsg = 'Density plots completed.'
    os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare replicate scatterplots
StatMsg = 'Analyzing replicate correlation ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for treatment in Treatments:
    os.system('python -u PrintStatus.py ProcessSample '+treatment+' 2>&1 | tee -a PinAPL-Py.log')
    if len(Replicates[treatment]) >= 2: 
        PairIt = itertools.combinations(Replicates[treatment],2)
        for pair in PairIt:
            os.system('python -u '+ReplicateScript+'.py '+pair[0]+' '+pair[1]+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Replicate correlation analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare heatmap
StatMsg = 'Analyzing sample clusters ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
os.system('python -u '+ClusterScript+'.py'+' 2>&1 | tee -a PinAPL-Py.log' )
DoneMsg = 'Clustering analysis completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Rank genes 
StatMsg = 'Gene ranking ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in TreatmentSamples:
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+GeneRankScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')
DoneMsg = 'Gene ranking completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Prepare gene score scatterplots
StatMsg = 'Plotting gene ranks ...'
os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
for sample in TreatmentSamples:
    os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')
    os.system('python -u '+GenePlotScript+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')    
DoneMsg = 'Gene rank scatterplots completed.'
os.system('python -u PrintStatus.py Done "'+DoneMsg+'" 2>&1 | tee -a PinAPL-Py.log')

# Auto-highlight top 10 genes
if AutoHighlight:
    StatMsg = 'Extracting top 10 ranked genes ...'
    os.system('python -u PrintStatus.py SubHeader "'+StatMsg+'" 2>&1 | tee -a PinAPL-Py.log')
    for sample in TreatmentSamples:
        os.system('python -u PrintStatus.py ProcessSample '+sample+' 2>&1 | tee -a PinAPL-Py.log')  
        os.system('python -u '+ExtractTop10Script+'.py '+sample+' 2>&1 | tee -a PinAPL-Py.log')         
        Top10Genes = GetTop10Genes(sample)
        for gene in Top10Genes:
            print('--------------------------------')
            print('Plotting scores: '+gene)
            print('--------------------------------')
            os.system('python -u '+GenePlotScript+'.py '+sample+' '+gene)
            os.system('python -u '+zFCScript+'.py '+sample+' '+gene)
            os.system('python -u '+vFCScript+'.py '+sample+' '+gene)
            os.system('python -u '+ScatterScript+'.py '+sample+' '+gene)    

# Combine gene ranks across replicates
if IncludeGeneRankCombination:
    os.system('python -u PrintStatus.py CombineReplicates '+sample+' 2>&1 | tee -a PinAPL-Py.log') 
    for treatment in Treatments:  
        if 'Control' not in treatment:            
            os.system('python -u '+CombineScript+'.py '+treatment+' 2>&1 | tee -a PinAPL-Py.log')        
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
    
# Copy Log Files 
if not os.path.exists(LogFileDir):
	os.makedirs(LogFileDir) 
os.system('cp configuration.yaml '+LogFileDir)
os.system('mv PinAPL-Py.log '+LogFileDir)
os.chdir(WorkingDir)
os.system('cp DataSheet.xlsx '+LogFileDir)

# Check for errors in log file
os.chdir(LogFileDir)
LogFile = open('PinAPL-Py.log','r')
LogFileText = LogFile.read()
if 'Traceback' in LogFileText:
    print('### ERROR(S) DETECTED IN PROGRAM LOG! ###')
    print('### Please check the log file (PinAPL-Py.log) ###')
    print('### You can use the "Submit a Bug" button on the website to contact the developers team! ###')    
    ErrorFile = open('ErrorFound.txt','w')
    ErrorFile.write('Python runtime error detected. Please check program log')
    ErrorFile.close()    
else:
    print('*** RUN FINISHED SUCCESSFULLY ***')
    print('LOADING RESULTS PAGE. PLEASE REFRESH PERIODICALLY...')