#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 19:11:28 2016

@author: philipp
"""

# Perform Bowtie2 alignment
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import pandas
from Bowtie2 import RunBowtie2
import yaml
import os
import glob
import time
import sys


def ReadAlignment(sample):    
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++++++++++++')  
    start_total = time.time()   

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    WorkingDir = config['WorkingDir']
    TempDataDir = config['TempDataDir']
    bw2Dir = config['bw2Dir']
    IndexDir = config['IndexDir']
    LibDir = config['LibDir']
    AlnStemDir = config['AlignDir']
    AlnDir = AlnStemDir+sample+'/'
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    if LibFormat == 'tsv':
        libsep = '\t'
    elif LibFormat == 'csv':
        libsep = ','
    L_bw = config['L_bw']
    N_bw = config['N_bw']
    i_bw = config['i_bw']

    
    # ------------------------------------------------
    # Read library
    # ------------------------------------------------  
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pandas.read_table(LibFilename, sep = libsep, skiprows = 1, names = LibCols)
    LibFile = LibFile.sort_values(['gene','ID'])    
    sgIDs = list(LibFile['ID'])
    global L
    L = len(sgIDs)
    global geneIDs
    geneIDs = list(LibFile['gene'])
    G = len(set(geneIDs))
    
    # ------------------------------------------------
    # Get sample read file
    # ------------------------------------------------      
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    n = len(FileNames)
    Samples = list(DataSheet['SAMPLE NAME'].values)
    ReadsFilename = [FileNames[j] for j in range(n) if Samples[j] == sample][0] 
    ReadsFilename0 = 'Trim_'+ReadsFilename    

    # ----------------------------------------------
    # Run alignment
    # ----------------------------------------------                  
    start = time.time()  
    print('Aligning reads to library ...')        
    RunBowtie2(ReadsFilename0,TempDataDir,AlnDir,bw2Dir,IndexDir,L_bw,N_bw,i_bw)
    print('Alignment completed.')
    end = time.time()
    # Time stamp
    aln_time = end-start
    if aln_time < 60: 
        time_elapsed = aln_time
        print('Time elapsed (Alignment) [secs]: ' + '%.3f' % time_elapsed)
    elif aln_time < 3600:
        time_elapsed = aln_time/60
        print('Time elapsed (Alignment) [mins]: ' + '%.3f' % time_elapsed)
    else:
        time_elapsed = aln_time/3600
        print('Time elapsed (Alignment) [hours]: ' + '%.3f' % time_elapsed)


    # --------------------------------------
    # Final time stamp
    # --------------------------------------        
    os.chdir(ScriptsDir)    
    end_total = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')    
    sec_elapsed = end_total - start_total
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Total) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Total) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Total) [hours]: ' + '%.3f' % time_elapsed +'\n')            
      
if __name__ == "__main__":
    input1 = sys.argv[1]
    ReadAlignment(input1)
