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

os.chdir('/workingdir/')
configFile = open('configuration.yaml','r')
config = yaml.load(configFile)
configFile.close()    
LibDir = config['LibDir']
IndexDir = config['IndexDir']
DataDir = config['DataDir']
ScriptsDir = config['ScriptsDir']
CtrlPrefix = config['CtrlPrefix']
script00 = config['script00']
script01 = config['script01']
script02 = config['script02']
script03 = config['script03']
script04 = config['script04']
script05 = config['script05']
script06 = config['script06']
script07 = config['script07']
script08 = config['script08']


print('*******************************************')
print('Launching PinAPL-Py ...')
print('P. Spahn et al., UC San Diego (10/2016)')
print('*******************************************')
start = time.time()

# Generate index if not present
if not os.path.exists(IndexDir):
    os.chdir(ScriptsDir)
    print('++++++++++++++++++++++++++++++++++++++++')
    print('Building Library Index ...')
    print('++++++++++++++++++++++++++++++++++++++++')
    command = 'python '+script00+'.py'
    os.system(command)  
    print('Library Index completed.')    
    print('\n\n')

# Read Samples
print('++++++++++++++++++++++++++++++++++++++++')
print('Reading Samples ...')
print('++++++++++++++++++++++++++++++++++++++++')
os.chdir(DataDir)
Filenames = [filename for filename in os.listdir(DataDir) if '.fastq.gz' in filename]
SampleNames = list()
for filename in Filenames:
    regexp = re.search('[a-zA-Z0-9/-]*_R[0-9a-z]*',filename)
    SampleNames.append(regexp.group())
    print('Found sample '+SampleNames[-1])
ToxSamples = [sample for sample in SampleNames if CtrlPrefix not in sample]
Treatments = list()
for tox_sample in ToxSamples:
    U = [pos for pos, char in enumerate(tox_sample) if char == '_']
    Treatments.append(tox_sample[0:U[0]])
Treatments = list(set(Treatments))
Replicates = dict()
for treatment in Treatments:
    replist = [sample for sample in SampleNames if treatment in sample]
    Replicates[treatment] = replist
print('\n\n')
    
# Align Reads
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Aligning Reads ...')
print('++++++++++++++++++++++++++++++++++++++++')
for sample in SampleNames:
    print('\nProcessing sample '+sample+' ... ')
    command = 'python '+script01+'.py '+sample
    os.system(command)
print('Read Alignment completed.')
print('\n\n')
        
# Analyze Counts
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Analyzing Read Counts ...')
print('++++++++++++++++++++++++++++++++++++++++')
for sample in SampleNames:
    print('\nProcessing sample '+sample+' ... ')
    command = 'python '+script02+'.py '+sample
    os.system(command)
print('Read Count Analysis completed.')
print('\n\n')

# Analyze Controls
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Analyzing Control Samples ...')
print('++++++++++++++++++++++++++++++++++++++++')
command = 'python '+script03+'.py'
os.system(command)
print('Control Analysis completed.')
print('\n\n')

# List Candidate sgRNA
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Preparing Candidate sgRNA Lists ...')
print('++++++++++++++++++++++++++++++++++++++++')
for tox_sample in ToxSamples:
    command = 'python '+script04+'.py '+tox_sample
    os.system(command)
print('Candidate sgRNA Lists completed.')
print('\n\n')

# List Candidate genes
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Preparing Candidate Gene Lists ...')
print('++++++++++++++++++++++++++++++++++++++++')
for tox_sample in ToxSamples:
    command = 'python '+script05+'.py '+tox_sample
    os.system(command)
print('Candidate Gene Lists completed.')
print('\n\n')

# Prepare sample scatterplots
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Preparing Sample Scatterplots ...')
print('++++++++++++++++++++++++++++++++++++++++')
for tox_sample in ToxSamples:
    print('\nProcessing sample '+tox_sample+' ... ')
    command = 'python '+script06+'.py '+tox_sample
    os.system(command)
print('Sample Scatterplots completed.')
print('\n\n')

# Prepare replicate scatterplots
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Preparing Replicate Scatterplots ...')
print('++++++++++++++++++++++++++++++++++++++++')
for treatment in Treatments:
    print('\nProcessing '+treatment+' ... ')
    if len(Replicates[treatment]) >= 2: 
        PairIt = itertools.combinations(Replicates[treatment],2)
        for pair in PairIt:
            command = 'python '+script07+'.py '+pair[0]+' '+pair[1]
            os.system(command)
print('Replicate Scatterplots completed.')
print('\n\n')

# Prepare heatmap
os.chdir(ScriptsDir)
print('++++++++++++++++++++++++++++++++++++++++')
print('Preparing Heatmap ...')
print('++++++++++++++++++++++++++++++++++++++++')
command = 'python '+script08+'.py' 
os.system(command)
print('Heatmap completed.')
print('\n\n')

# Final Time Stamp
end = time.time()
print('------------------------------------------------')
print('------------------------------------------------')
print('PinAPL-Py completed.')    
sec_elapsed = end - start
if sec_elapsed < 60:
    time_elapsed = sec_elapsed
    print('Time elapsed [secs]: ' + '%.3f' % time_elapsed +'\n')
elif sec_elapsed < 3600:
    time_elapsed = sec_elapsed/60
    print('Time elapsed [mins]: ' + '%.3f' % time_elapsed +'\n')
else:
    time_elapsed = sec_elapsed/3600
    print('Time elapsed [hours]: ' + '%.3f' % time_elapsed +'\n') 
